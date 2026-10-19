#!/usr/bin/env python3
"""Ensure that the exact checker rejects three deliberately corrupted certificates."""
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT=Path(__file__).resolve().parent
FILES=['verify.py','sos_certificate.json','certificate_basis0.json','certificate_basis1.json',
       'certificate_basis2.json','relation_generators.json','equality_certificate.json']

def run_test(name, mutate, expected):
    with tempfile.TemporaryDirectory(prefix='chollet-mutation-') as td:
        root=Path(td)
        for f in FILES:shutil.copy2(ROOT/f,root/f)
        mutate(root)
        result=subprocess.run([sys.executable,str(root/'verify.py')],capture_output=True,text=True,timeout=30)
        text=result.stdout+result.stderr
        assert result.returncode!=0 and expected in text, (name,result.returncode,text)
        return {'test':name,'result':'PASS: corruption rejected','expected_failure':expected}

def alter_identity(root):
    path=root/'sos_certificate.json';d=json.loads(path.read_text());qd=int(d['q_denominator'])
    d['blocks'][0]['Q_numerators'][0][0]=str(int(d['blocks'][0]['Q_numerators'][0][0])+qd)
    path.write_text(json.dumps(d,separators=(',',':')))

def alter_positivity(root):
    path=root/'sos_certificate.json';d=json.loads(path.read_text());qd=int(d['q_denominator'])
    d['blocks'][0]['Q_numerators'][0][0]=str(int(d['blocks'][0]['Q_numerators'][0][0])-1000000*qd)
    path.write_text(json.dumps(d,separators=(',',':')))

def alter_equality(root):
    path=root/'equality_certificate.json';d=json.loads(path.read_text())
    for a in d['identities'][0]['coefficients']:
        if a[0]!='0':
            a[0]=str(int(a[0])+int(a[1]));break
    path.write_text(json.dumps(d,separators=(',',':')))

if __name__=='__main__':
    results=[run_test('positivity',alter_positivity,'Q >= (1/5)I not certified'),
             run_test('global polynomial identity',alter_identity,'Polynomial identity failed'),
             run_test('equality implications',alter_equality,'Equality implication 0 failed')]
    for r in results:print(r['test']+': '+r['result'])
    (ROOT/'mutation_test_report.json').write_text(json.dumps({'status':'PASS','tests':results},indent=2)+'\n')
