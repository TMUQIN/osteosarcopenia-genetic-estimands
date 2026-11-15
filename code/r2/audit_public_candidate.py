"""Inspect candidate file content without sending data or printing secret matches."""
import argparse,csv,hashlib,json,re,sys
from pathlib import Path
FORBIDDEN_SUFFIXES={'.rds','.RDS','.vcf','.bim','.bed','.fam','.bgz','.parquet','.sqlite','.gz','.zip','.xlsx','.docx','.pdf'}
PATTERNS={
 'github_token':r'gh[pousr]_[A-Za-z0-9]{20,}',
 'aws_access_key':r'AKIA[0-9A-Z]{16}',
 'private_key':r'-----BEGIN (?:RSA |OPENSSH |EC )?PRIVATE KEY-----',
 'credential_assignment':r'(?i)(?:access_token|api_key|password|client_secret)\s*[:=]\s*[\x22\x27][^\x22\x27\s]{12,}[\x22\x27]',
 'personal_absolute_path':r'(?i)(?:C|D):[\\/](?:Users|Research project)[\\/]',
}
def audit(root):
    issues=[];files=[]
    baseline_path=root/'r2/BASELINE_FILE_IDENTITIES.tsv'
    baseline={}
    if baseline_path.exists():
        with baseline_path.open(encoding='utf-8',newline='') as f:
            baseline={r['relative_path']:r['sha256'] for r in csv.DictReader(f,delimiter='\t')}
    retained=[]
    for p in sorted(root.rglob('*')):
        if not p.is_file() or '.git' in p.relative_to(root).parts or '__pycache__' in p.parts or p.suffix=='.pyc':continue
        if any(part.startswith('local_') for part in p.relative_to(root).parts):continue
        rel=p.relative_to(root).as_posix()
        data=p.read_bytes()
        same_as_baseline=baseline.get(rel)==hashlib.sha256(data).hexdigest()
        inherited_artwork=(same_as_baseline and rel.startswith(('code/figures/submitted_figure3/','figures/submitted_artwork/'))
                           and p.suffix in {'.pdf','.png','.eps','.svg'})
        if same_as_baseline:retained.append(rel)
        if p.suffix in FORBIDDEN_SUFFIXES and not inherited_artwork:issues.append(dict(file=rel,reason='FORBIDDEN_ARTIFACT_TYPE'))
        if p.stat().st_size>2_000_000 and not inherited_artwork:issues.append(dict(file=rel,reason='UNEXPECTED_LARGE_PUBLIC_FILE'))
        try:txt=data.decode('utf-8-sig')
        except UnicodeDecodeError:
            if not inherited_artwork:issues.append(dict(file=rel,reason='UNEXPECTED_BINARY_FILE'))
            files.append(rel);continue
        for label,pattern in PATTERNS.items():
            if re.search(pattern,txt):issues.append(dict(file=rel,reason=label))
        if p.suffix=='.json':
            obj=json.loads(txt)
            if isinstance(obj,dict) and 'matrices' in obj and obj.get('data_kind')!='SYNTHETIC_NOT_STUDY_DATA':
                issues.append(dict(file=rel,reason='UNLABELLED_OR_REAL_MATRIX'))
        if p.suffix in ('.csv','.tsv'):
            header=txt.splitlines()[0] if txt else ''
            if any(x in header.split('\t') for x in ('S_moment','V_diagonal','N_pairwise','Raw_BETA','Raw_P','variant_id')):
                issues.append(dict(file=rel,reason='POTENTIAL_RESTRICTED_MATRIX_OR_VARIANT_TABLE'))
        files.append(rel)
    return dict(status='PASS_CONTENT_EXCLUSION_SCAN' if not issues else 'FAIL',files_checked=len(files),issues=issues,
        scope='All payload files; unchanged licensed baseline artwork allowed only by exact SHA256; new restricted data excluded',
        baseline_unchanged_files=len(retained),license_clearance=False,remote_v1_baseline_compared=bool(baseline),
        limitations='Pattern and artifact-type scan complements the explicit file inventory; it does not grant rights. The verified current candidate tree is checked, not every historical Git commit.')
def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2]);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    obj=audit(a.root);a.output.write_text(json.dumps(obj,indent=2)+'\n');print(json.dumps(obj,indent=2));return 0 if not obj['issues'] else 1
if __name__=='__main__':sys.exit(main())
