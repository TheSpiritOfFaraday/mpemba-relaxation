"""Negative controls for exact decimal binding in the physical proof runner."""
from decimal import Decimal
import json
from pathlib import Path
from run_robust_finite_time_certificate import verify_embedded_parameters,ROOT,SOURCES
checks=[]
def check(name,condition):
 if not condition:raise AssertionError(name)
 checks.append(name)
def source(values):return 'const char*params[]={'+','.join('"'+str(v)+'"' for v in values)+'};'
expected=[Decimal('1.00000000000000001')]+[Decimal(str(i)) for i in range(2,9)]
check('Exact parameters accepted',verify_embedded_parameters(source(expected),expected)==expected)
near=expected.copy();near[0]=Decimal('1.00000000000000002')
check('Regression values collapse to the same binary64 value',float(near[0])==float(expected[0]))
for name,text in [('Different decimals with equal floats rejected',source(near)),('Missing parameter rejected',source(expected[:-1])),('Nonfinite parameter rejected',source([Decimal('Infinity')]+expected[1:])),('Missing array rejected','int main(){}')]:
 rejected=False
 try:verify_embedded_parameters(text,expected)
 except ValueError:rejected=True
 check(name,rejected)
live=json.loads((ROOT/'research/physical_rate_refined.json').read_text(),parse_float=Decimal,parse_int=Decimal)['result']['parameters']
for name,p in SOURCES.items():check('Exact current input matches '+name,verify_embedded_parameters(p.read_text(),live)==live)
out={'passed':True,'count':len(checks),'checks':checks,'finding':'Runner float-only equality could hide a decimal input mismatch. Current inputs already match exactly; completion audit had a separate exact Decimal check. Runner now rejects that mismatch before compilation.'}
(ROOT/'research/detailed_audit/physical_runner_input_regressions.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
