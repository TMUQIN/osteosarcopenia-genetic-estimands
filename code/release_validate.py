#!/usr/bin/env python3
import csv, pathlib, math
root=pathlib.Path(__file__).resolve().parents[1]
data=root/'data'/'derived_public_candidate'
def read(name):
 with (data/name).open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
for name in ['MODEL_PARAMETERS_ALL_SCENARIOS.tsv','MODEL_FIT_ALL_SCENARIOS.tsv']:
 r=read(name);assert r and set(x['scenario'] for x in r)=={'D','R75','R50'},name
 assert not any('V201' in str(x) for x in r),name
 assert all('source_model_rds' not in x for x in r),name
expected_labels={'D':'Primary specification','R75':'75% assumed-N sensitivity','R50':'50% assumed-N sensitivity'}
for name in ['MODEL_PARAMETERS_ALL_SCENARIOS.tsv','MODEL_FIT_ALL_SCENARIOS.tsv','FIGURE2_DATA.tsv','FIGURE_S1_DATA.tsv']:
 r=read(name)
 for row in r:
  assert row['scenario_label']==expected_labels[row['scenario']],(name,row['scenario'],row['scenario_label'])
for r in read('MODEL_FIT_ALL_SCENARIOS.tsv'):
 assert 'AIC' not in r and 'p_chisq' not in r
 assert float(r['df'])==11
for r in read('BONE_OUTCOME_SHIFTS_DATA.tsv'):
 assert r['scenario']=='D'
 assert not any(k in r for k in ('expected_direction','moves_toward_ALM','direction_matches_prediction','distance_component_total','distance_component_specific'))
assert len(read('BONE_OUTCOME_SHIFTS_DATA.tsv'))==6
assert len(read('BONE_ASSUMED_N_SENSITIVITY.tsv'))==12
source={r['side']:r for r in read('FIGURE2_DATA.tsv') if r['scenario']=='D'}
summary={r['side']:r for r in read('DELTACOV_VERIFICATION_SUMMARY.tsv')}
assert set(summary)=={'right','left'}
for side in ('right','left'):
 a=source[side];b=summary[side]
 assert abs(float(a['paired_covariance_difference'])-float(b['DeltaCov_g']))<1e-12
 est,se,lo,hi=[float(b[k]) for k in ('DeltaCov_g','SE','CI95_low','CI95_high')]
 assert abs((lo+hi)/2-est)<1e-10
 assert abs((hi-lo)/(2*1.96)-se)<1e-8
 assert 0<float(b['nominal_P'])<1
 assert a['paired_covariance_difference_uncertainty'].startswith('VERIFIED_INDEPENDENT')
plot=root/'code'/'figures'/'submitted_figure3'/'FIGURE3_PLOT_SOURCE.tsv'
assert plot.exists() and len(list(csv.DictReader(plot.open(),delimiter='\t')))==12
print('CANDIDATE_RELEASE_CONTENT_VALIDATION=PASS')
print('primary_sides=2; frozen_scenarios=D,R75,R50; bone_primary=6; bone_assumed_N=12; cov_robust_inferences=2')
print('NO_PRIMARY_GWAS_LDSC_GENOMICSEM_RERUN')
