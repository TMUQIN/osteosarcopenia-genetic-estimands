"""Portable I/O adapter around exact recovered R2 functions; no GWAS/LDSC/SEM calls."""
import argparse,csv,hashlib,json,math,platform,subprocess
from pathlib import Path
import numpy as np
import scipy
from scipy.stats import norm
import r2_closure_core as core
import r2_wald_core as wald

KEEP=np.array([1,2,3,4,5,9,10,11,12,16,17,18,22,23,27])-1
WEIGHTS=[0,.25,.5,.75,1]
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def save(path,rows):
    cols=list(dict.fromkeys(k for r in rows for k in r))
    with Path(path).open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols,lineterminator='\n');w.writeheader();w.writerows(rows)
def load_rds(path,rscript):
    expr="x<-readRDS(commandArgs(TRUE)[1]); stopifnot(identical(colnames(x$S),c('ALM','ArmFFM','LegFFM','ArmFM','LegFM','eBMD','FNBMD','LSBMD')));cat(jsonlite::toJSON(lapply(x[c('S','V')],function(a)list(dim=dim(a),values=sprintf('%.17g',as.vector(a)))),auto_unbox=TRUE))"
    run=subprocess.run([rscript,'-e',expr,str(path)],capture_output=True,text=True,check=True)
    obj=json.loads(run.stdout)
    return {k:np.array(v['values'],dtype=float).reshape(v['dim'],order='F') for k,v in obj.items()}
def validate(s,v):
    if s.shape!=(8,8) or v.shape!=(36,36):raise ValueError('Expected S8x8 and V36x36')
    for a in (s,v):
        if not np.isfinite(a).all() or not np.allclose(a,a.T,rtol=0,atol=1e-14):raise ValueError('Nonfinite/asymmetric matrix')
        if np.linalg.eigvalsh(a).min()<=0:raise ValueError('Matrix must be positive definite; no repair applied')
def five_trait(s,v):
    validate(s,v)
    return s[:5,:5],v[np.ix_(KEEP,KEEP)]
def evaluate(s,v,side):
    validate(s,v);five_trait(s,v);x=core.vech_from_matrix(s)
    vals=wald.composite(x);keys=list(vals)
    jac,_,_=wald.gradients(lambda z:np.array(list(wald.composite(z).values())),x)
    cov=jac@v@jac.T
    if np.any(np.diag(cov)<=0):raise ValueError('Nonpositive derived sampling variance')
    comps=[dict(side=side,metric=k,**wald.single(np.real(vals[k]),math.sqrt(cov[i,i]))) for i,k in enumerate(keys)]
    inds=[keys.index(k+'_delta_rg') for k in ('eBMD','FNBMD','LSBMD')]
    shifts=np.array([np.real(vals[keys[k]]) for k in inds]);scov=cov[np.ix_(inds,inds)]
    eq=np.array([[1.,-1.,0.],[1.,0.,-1.]])
    globals_=[dict(side=side,test='W3_joint_zero',**wald.wald(shifts,scov)),dict(side=side,test='W2_equal_shifts',**wald.wald(eq@shifts,eq@scov@eq.T))]
    pairs=[]
    for i,j in ((0,1),(0,2),(1,2)):
        est=shifts[i]-shifts[j];se=math.sqrt(scov[i,i]+scov[j,j]-2*scov[i,j])
        pairs.append(dict(side=side,contrast=['eBMD','FNBMD','LSBMD'][i]+'_minus_'+['eBMD','FNBMD','LSBMD'][j],**wald.single(est,se)))
    grid=[];std=[];partial=[];checks=[];restrictions=[]
    for w in WEIGHTS:
        for m in WEIGHTS:
            fn=lambda z:np.array(list(core.composite(z,w,m).values()))
            values=core.composite(x,w,m);j=core.gradient(fn,x);cv=j@v@j.T
            checks.append(dict(side=side,**core.gradient_crosscheck(fn,x,f'grid_{w}_{m}')))
            if np.any(np.diag(cv)<=0):raise ValueError('Invalid grid sampling variance')
            for i,k in enumerate(values):grid.append(dict(side=side,w=w,v=m,metric=k,**core.estimate(values[k],math.sqrt(cv[i,i]))))
    for tag,fn,rows in [('genetic_SD',lambda z:core.composite(z,.5,.5,True),std),('partial_rg',core.partial_rg,partial)]:
        values=fn(x);j=core.gradient(lambda z:np.array(list(fn(z).values())),x);cv=j@v@j.T
        checks.append(dict(side=side,**core.gradient_crosscheck(lambda z:np.array(list(fn(z).values())),x,tag)))
        for i,k in enumerate(values):rows.append(dict(side=side,metric=k,**core.estimate(values[k],math.sqrt(cv[i,i]))))
    fv=core.free_restrictions(x);fj=core.gradient(core.free_restrictions,x);fc=fj@v@fj.T
    for tag,ix in [('ALM_only_joint_3df',[0,1,5]),('ALM_external_joint_2df',[1,5]),('cross_block_tetrad_1df',[0])]:
        restrictions.append(dict(side=side,test=tag,**core.joint_wald(fv[ix],fc[np.ix_(ix,ix)])))
    return dict(composite=comps,bone_global=globals_,bone_pairwise=pairs,weight_grid=grid,genetic_sd=std,partial_rg=partial,alm_restrictions=restrictions,gradient_checks=checks)
def main():
    p=argparse.ArgumentParser();p.add_argument('--input-dir',required=True,type=Path);p.add_argument('--output-dir',required=True,type=Path)
    p.add_argument('--format',choices=['rds','synthetic-json'],default='rds');p.add_argument('--rscript',default='Rscript');args=p.parse_args()
    args.output_dir.mkdir(parents=True,exist_ok=True);allrows={};inputs=[]
    for side in ('right','left'):
        path=args.input_dir/(side+'_joint_8trait_ldsc.rds' if args.format=='rds' else side+'.json')
        obj=load_rds(path,args.rscript) if args.format=='rds' else {k:np.array(v,dtype=float) for k,v in json.loads(path.read_text())['matrices'].items()}
        for k,rows in evaluate(obj['S'],obj['V'],side).items():allrows.setdefault(k,[]).extend(rows)
        inputs.append(dict(side=side,input_name=path.name,sha256=sha(path),data_kind='REAL_FROZEN_AUTHOR_INPUT_NOT_DISTRIBUTED' if args.format=='rds' else 'SYNTHETIC_NOT_STUDY_DATA'))
    for key,rows in allrows.items():save(args.output_dir/(key+'.csv'),rows)
    save(args.output_dir/'input_identities.csv',inputs)
    env=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,mode=args.format,model_calls=0)
    (args.output_dir/'environment.json').write_text(json.dumps(env,indent=2)+'\n')
    print(json.dumps(dict(environment=env,row_counts={k:len(v) for k,v in allrows.items()},inputs=inputs),indent=2))
if __name__=='__main__':main()
