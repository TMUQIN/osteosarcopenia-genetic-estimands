"""Explicit keyed validation of fresh outputs against retained full-precision R2 references."""
import csv,json,math
from pathlib import Path
import numpy as np
import argparse
p=argparse.ArgumentParser();p.add_argument('--project-root',type=Path,required=True);p.add_argument('--recomputed-dir',type=Path,required=True);p.add_argument('--independent-r-dir',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);args=p.parse_args()
ROOT=args.project_root;OUT=args.output_dir;OUT.mkdir(parents=True,exist_ok=True);PUB=OUT;(PUB/'logs').mkdir(exist_ok=True)
def read(p,sep='\t'):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter=sep))
checks=[]
def compare(label,actual,expected,keys,fields):
    def keyed(rows):
        d={}
        for r in rows:
            k=tuple(str(r[x]) for x in keys);assert k not in d,(label,'duplicate',k);d[k]=r
        return d
    a,b=keyed(actual),keyed(expected);assert a.keys()==b.keys(),(label,'keys',a.keys()-b.keys(),b.keys()-a.keys())
    for k in a:
        for ak,bk in fields:
            av,bv=float(a[k][ak]),float(b[k][bk])
            ok=math.isclose(av,bv,rel_tol=1e-10,abs_tol=1e-12)
            checks.append(dict(comparison=label,key='|'.join(k),actual_field=ak,reference_field=bk,actual=av,reference=bv,absolute_difference=abs(av-bv),status='PASS' if ok else 'FAIL'))
    assert all(x['status']=='PASS' for x in checks),label
recomputed=args.recomputed_dir
old=ROOT/'R2_INDEPENDENT_WALD_VERIFICATION_20261110'
fields=[(x,x) for x in ('estimate','SE','z','P','CI_low','CI_high','P_bonf_m3','P_bonf_m6')]
compare('equal_composite',read(recomputed/'composite.csv',','),read(old/'R2_MANIFEST_COMPOSITE_RESULTS.tsv'),['side','metric'],fields)
compare('bone_global_W2_W3',read(recomputed/'bone_global.csv',','),read(old/'R2_MANIFEST_BONE_GLOBAL_TESTS.tsv'),['side','test'],[(x,x) for x in ('Q','df','P','min_eigen','max_eigen','condition')])
compare('bone_pairwise',read(recomputed/'bone_pairwise.csv',','),read(old/'R2_MANIFEST_BONE_PAIRWISE.tsv'),['side','contrast'],fields)
close=ROOT/'JCSM_R2_FINAL_EVIDENCE_CLOSURE'
for actual,reference in [('genetic_sd','R2_GENETIC_SD_SENSITIVITY.tsv'),('partial_rg','R2_PARTIAL_RG_SENSITIVITY.tsv')]:
    expected=[dict(x,metric=x['estimand']) for x in read(close/reference)]
    compare(actual,read(recomputed/(actual+'.csv'),','),expected,['side','metric'],[(x,x) for x in ('estimate','SE','CI95_low','CI95_high','P_nominal','log10_P_nominal')])
expected=[]
for x in read(close/'R2_WEIGHT_GRID_25.tsv'):
    for k in ('gamma','Var_g_F','Var_g_U','rg_F_ALM','rg_U_ALM','Delta_rg','Delta_Cov','variance_retention','covariance_retention'):
        expected.append(dict(side=x['side'],w=x['weight_w_ArmFFM'],v=x['weight_v_ArmFM'],metric=k,**{stat:x[k+'_'+stat] for stat in ('estimate','SE','CI95_low','CI95_high','P_nominal','log10_P_nominal')}))
actual=read(recomputed/'weight_grid.csv',',')
for rows in (actual,expected):
    for x in rows:
        x['w']=str(float(x['w']));x['v']=str(float(x['v']))
compare('all_25x2_grid',actual,expected,['side','w','v','metric'],[(x,x) for x in ('estimate','SE','CI95_low','CI95_high','P_nominal','log10_P_nominal')])
compare('ALM_subset_Wald',read(recomputed/'alm_restrictions.csv',','),read(close/'R2_ALM_ONLY_FREE_LOADING_WALD.tsv'),['side','test'],[(x,x) for x in ('Q','df','P_nominal')])
with (OUT/'FROZEN_RESULTS_NUMERIC_CHECKS.tsv').open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(checks[0]),delimiter='\t');w.writeheader();w.writerows(checks)
summary=dict(status='PASS_LOCAL_FROZEN_REFERENCE_CHECK_ONLY',numeric_comparisons=len(checks),failed=sum(x['status']!='PASS' for x in checks),
    tolerances=dict(rtol=1e-10,atol=1e-12,basis='Existing R2 closure comparison criterion; no target adjustment'),submitted_15_sheet_excel_verified=False,
    families={k:sum(x['comparison']==k for x in checks) for k in dict.fromkeys(x['comparison'] for x in checks)})
# The independent R adapter is run separately; compare all of its 470 rows.
r_rows=read(args.independent_r_dir/'R2_INDEPENDENT_R_CHECK.tsv');r_checks=[]
for r in r_rows:
    kind={'grid':'weight_grid','genetic_SD':'genetic_sd','partial_rg':'partial_rg'}[r['analysis']]
    source=read(recomputed/(kind+'.csv'),',')
    matched=[x for x in source if x['side']==r['side'] and x['metric']==r['metric'] and (kind!='weight_grid' or (float(x['w'])==float(r['w']) and float(x['v'])==float(r['v'])))]
    assert len(matched)==1
    a=matched[0]
    for field in ('estimate','SE'):
        ok=math.isclose(float(a[field]),float(r[field]),rel_tol=1e-8,abs_tol=1e-11)
        r_checks.append(dict(side=r['side'],analysis=r['analysis'],w=r['w'],v=r['v'],metric=r['metric'],field=field,status='PASS' if ok else 'FAIL',absolute_difference=abs(float(a[field])-float(r[field]))))
assert len(r_rows)==470 and all(x['status']=='PASS' for x in r_checks)
summary['independent_R_rows']=len(r_rows);summary['independent_R_numeric_comparisons']=len(r_checks)
summary['independent_R_tolerance']=dict(rtol=1e-8,atol=1e-11,basis='Independent Richardson finite-difference check, distinct from exact source-function comparison')
with (OUT/'INDEPENDENT_R_NUMERIC_CHECKS.tsv').open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(r_checks[0]),delimiter='\t');w.writeheader();w.writerows(r_checks)
(PUB/'logs/FROZEN_NUMERIC_VALIDATION.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
