"""Negative controls for the new exact arithmetic implementation."""
from pathlib import Path
import json
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[1]
source=(ROOT/'research/check_freshman_walkthrough.py').read_text()
mutations=[
    ('wrong_upward_balance', 'rate*q[j]/q[i]', 'rate*q[i]/q[j]', 'balances'),
    ('wrong_slow_eigenvalue', '-F(4966671,250000)', '-F(4966672,250000)', 'left eigenpair'),
    ('wrong_overlap_sign', "'-1/250'", "'1/250'", 'bath-orthogonal mode'),
    ('invalid_fast_root_interval', '(-60,-59)', '(-59,-58)', 'one fast root'),
    ('reversed_temperature_preparations', 'hot,Dh,ph=preparation(1000); cold,Dc,pc=preparation(100)',
     'hot,Dh,ph=preparation(100); cold,Dc,pc=preparation(1000)', 'printed hot bracket'),
]
records=[]
for optimized in [False,True]:
    for name,old,new,reason in mutations:
        if source.count(old)!=1: raise ValueError('Ambiguous mutation '+name)
        mutated=source.replace(old,new)
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'bad_check.py';path.write_text(mutated)
            command=[sys.executable,'-I','-S']+(['-O'] if optimized else [])+[str(path)]
            run=subprocess.run(command,capture_output=True,text=True,timeout=30)
            if run.returncode==0 or reason not in run.stderr:
                raise ValueError((name,optimized,run.returncode,run.stderr))
        records.append({'mutation':name,'optimized':optimized,'rejected':True,'reason':reason})
result={'passed':True,'scope':'Five deliberate corruptions must fail in both ordinary and optimized Python.',
        'negative_controls':records}
(ROOT/'research/detailed_audit/freshman_negative_controls.json').write_text(json.dumps(result,indent=2)+'\n')
print('PASS: ten negative controls rejected at the intended check.')
