"""Fail-closed comparison of explicitly mapped Excel cells to derived CSV values.

Uses standard-library OOXML reading; does not modify workbooks or infer table identity
from sheet count. Every target cell must be explicitly mapped, including its row key.
"""
import argparse,csv,hashlib,json,math,posixpath,sys,zipfile
from pathlib import Path
from xml.etree import ElementTree as ET
NS={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main','r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
REQUIRED={'S2','S3','S3b','S3c','S4','S5'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def workbook(path):
    result={}
    with zipfile.ZipFile(path) as z:
        shared=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            shared=[''.join(e.itertext()) for e in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('m:si',NS)]
        rels={x.attrib['Id']:x.attrib['Target'] for x in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
        sheets=ET.fromstring(z.read('xl/workbook.xml')).findall('m:sheets/m:sheet',NS)
        for s in sheets:
            dest=rels[s.attrib['{'+NS['r']+'}id']]
            dest=dest.lstrip('/') if dest.startswith('/') else posixpath.normpath(posixpath.join('xl',dest))
            cells={}
            for c in ET.fromstring(z.read(dest)).findall('.//m:sheetData/m:row/m:c',NS):
                value=c.find('m:v',NS);kind=c.attrib.get('t');f=c.find('m:f',NS)
                text=value.text if value is not None else None
                if kind=='s' and text is not None:text=shared[int(text)]
                elif kind=='inlineStr':text=''.join(c.find('m:is',NS).itertext())
                cells[c.attrib['r']]=dict(value=text,formula=f is not None,type=kind)
            result[s.attrib['name']]=cells
    return result
def validate(excel,mapping,derived):
    book=workbook(excel)
    if len(book)!=15:raise ValueError(f'Expected author-designated 15-sheet workbook; observed {len(book)}')
    with mapping.open(encoding='utf-8-sig',newline='') as f:spec=list(csv.DictReader(f,delimiter='\t'))
    if not spec or not REQUIRED.issubset({r['table'] for r in spec}):raise ValueError('Explicit coverage for S2/S3/S3b/S3c/S4/S5 is missing')
    checks=[];cache={};seen=set()
    for r in spec:
        token=(r['sheet'],r['cell'])
        if token in seen:raise ValueError('Duplicate mapped Excel cell')
        seen.add(token);cell=book[r['sheet']][r['cell']]
        if cell['type']=='e' or cell['value'] is None:raise ValueError('Missing/error/cached-formula value: '+str(token))
        path=(derived/r['csv_file']).resolve()
        if not path.is_relative_to(derived.resolve()):raise ValueError('CSV path escapes derived directory')
        if path not in cache:
            with path.open(encoding='utf-8-sig',newline='') as f:cache[path]=list(csv.DictReader(f))
        key=json.loads(r['key_json']);matches=[x for x in cache[path] if all(x.get(k)==str(v) for k,v in key.items())]
        if len(matches)!=1:raise ValueError('CSV key is not unique: '+str(key))
        a=float(matches[0][r['field']]);target=float(cell['value']);atol=float(r['absolute_tolerance']);rtol=float(r['relative_tolerance'])
        if not math.isfinite(a) or not math.isfinite(target) or atol<0 or rtol<0:raise ValueError('Invalid numeric value or tolerance')
        if not r['precision_basis'].strip():raise ValueError('Tolerance precision basis is required')
        ok=math.isclose(a,target,rel_tol=rtol,abs_tol=atol)
        checks.append(dict(table=r['table'],sheet=r['sheet'],cell=r['cell'],csv_file=r['csv_file'],key=key,field=r['field'],computed=a,excel=target,status='PASS' if ok else 'FAIL',precision_basis=r['precision_basis']))
    # Coverage must be asserted by author mapping, not inferred from finding matching numbers.
    return dict(status='PASS_EXPLICITLY_MAPPED_CELLS' if all(r['status']=='PASS' for r in checks) else 'FAIL',workbook_sha256=sha(excel),mapping_sha256=sha(mapping),sheet_count=len(book),mapped_numeric_cells=len(checks),checks=checks,
                limitation='Explicit cell comparison only; completeness of main-result mapping must be separately reviewed.')
def main():
    p=argparse.ArgumentParser();p.add_argument('--excel',type=Path);p.add_argument('--mapping',type=Path);p.add_argument('--derived-dir',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if not a.excel or not a.mapping or not a.excel.is_file() or not a.mapping.is_file():
        obj=dict(status='NOT_VERIFIED',reason='Author-designated 15-sheet R2 workbook and explicit main-result cell map are required; no substitute workbook used.');rc=2
    else:
        try:obj=validate(a.excel,a.mapping,a.derived_dir);rc=0 if obj['status']=='PASS_EXPLICITLY_MAPPED_CELLS' else 1
        except Exception as e:obj=dict(status='FAIL',reason=str(e));rc=1
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');print(json.dumps(obj,ensure_ascii=False));return rc
if __name__=='__main__':sys.exit(main())
