"""Synthetic numerical tests only; none of these matrices are study estimates."""
import importlib.util,json,sys,tempfile,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code/r2'))
import run_r2 as run
import r2_closure_core as core
import r2_wald_core as wald

def synthetic(seed=274):
    rng=np.random.default_rng(seed)
    a=rng.normal(size=(8,8));s=a@a.T/40+np.eye(8)*.3
    a=rng.normal(size=(36,36));v=(a@a.T+np.eye(36)*10)*1e-7
    return s,v
class TestR2(unittest.TestCase):
    def setUp(self):self.s,self.v=synthetic()
    def test_round_trip_lower_triangle(self):
        sent=np.arange(64).reshape(8,8);sent=sent+sent.T
        np.testing.assert_array_equal(core.vech_matrix(core.vech_from_matrix(sent)),sent)
    def test_five_trait_both_axes_and_reject_top_left(self):
        s5,v5=run.five_trait(self.s,self.v)
        self.assertEqual(v5.shape,(15,15))
        self.assertEqual((run.KEEP+1).tolist(),[1,2,3,4,5,9,10,11,12,16,17,18,22,23,27])
        self.assertFalse(np.array_equal(v5,self.v[:15,:15]))
        for i,a in enumerate(run.KEEP):
            for j,b in enumerate(run.KEEP):self.assertEqual(v5[i,j],self.v[a,b])
    def test_orthogonality_and_retention(self):
        s=self.s;f=np.array([0,.5,.5,0,0,0,0,0]);m=np.array([0,0,0,.5,.5,0,0,0])
        gamma=(f@s@m)/(m@s@m);u=f-gamma*m
        self.assertAlmostEqual(u@s@m,0,places=14)
        ans=core.composite(core.vech_from_matrix(s),.5,.5)
        self.assertAlmostEqual(ans['gamma'],gamma)
        self.assertAlmostEqual(ans['Var_g_U'],f@s@f-(f@s@m)**2/(m@s@m))
    def test_paired_delta_uses_cross_covariance(self):
        x=core.vech_from_matrix(self.s)
        fn=lambda z:np.array([core.composite(z,.5,.5)[k] for k in ('rg_F_ALM','rg_U_ALM')])
        j=core.gradient(fn,x);cv=j@self.v@j.T
        delta=core.gradient(lambda z:np.diff(fn(z)),x)
        self.assertAlmostEqual(float((delta@self.v@delta.T)[0,0]),cv[0,0]+cv[1,1]-2*cv[0,1],places=15)
        self.assertNotAlmostEqual(cv[0,0]+cv[1,1],cv[0,0]+cv[1,1]-2*cv[0,1],places=10)
    def test_full_v_changes_se(self):
        x=core.vech_from_matrix(self.s);j=core.gradient(lambda z:core.composite(z,.5,.5)['Delta_rg'],x)
        self.assertGreater(abs(float((j@self.v@j.T-j@np.diag(np.diag(self.v))@j.T)[0,0])),1e-10)
    def test_sd_denominators_remain_in_gradient(self):
        x=core.vech_from_matrix(self.s)
        fn=lambda z:core.composite(z,.5,.5,True)['Delta_rg']
        j=core.gradient(fn,x)
        self.assertEqual(core.gradient_crosscheck(fn,x,'synthetic_SD')['status'],'PASS')
        self.assertGreater(np.max(abs(j[0,[core.PAIR_INDEX[(i,i)] for i in range(1,5)]])),1e-6)
    def test_all_grid_combinations_and_outputs(self):
        result=run.evaluate(self.s,self.v,'synthetic')
        self.assertEqual(len(result['weight_grid']),450//2)
        self.assertEqual(len({(x['w'],x['v']) for x in result['weight_grid']}),25)
        self.assertEqual({x['test'] for x in result['bone_global']},{'W2_equal_shifts','W3_joint_zero'})
        self.assertEqual(len(result['bone_pairwise']),3)
    def test_bad_input_is_rejected(self):
        with self.assertRaises(ValueError):run.validate(self.s,self.v[:15,:15])
        bad=self.v.copy();bad[0,0]=-1
        with self.assertRaises(ValueError):run.validate(self.s,bad)
    def test_synthetic_fixture_label(self):
        for side in ('right','left'):
            obj=json.loads((ROOT/'data/synthetic_r2'/f'{side}.json').read_text())
            self.assertEqual(obj['data_kind'],'SYNTHETIC_NOT_STUDY_DATA')
if __name__=='__main__':unittest.main(verbosity=2)
