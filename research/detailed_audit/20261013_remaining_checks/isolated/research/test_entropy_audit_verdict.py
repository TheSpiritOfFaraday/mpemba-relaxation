"""Exercise the audit's actual final verdict code with injected result counters.
Avoid importing the sampler, whose top level runs the full numerical experiment.
"""
import ast,copy,hashlib,json
from pathlib import Path
source=Path(__file__).with_name('finite_time_entropy_independent_audit.py')
tree=ast.parse(source.read_text());start=next(i for i,n in enumerate(tree.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='failures' for t in n.targets))
stop=next(i for i in range(start,len(tree.body)) if isinstance(tree.body[i],ast.Expr))
raise_node=tree.body[-1]
assert isinstance(raise_node,ast.If) and isinstance(raise_node.body[0],ast.Raise)
program=compile(ast.Module(body=tree.body[start:stop]+[raise_node],type_ignores=[]),str(source),'exec')
base={k:{'KL':{'negative':0,'skipped_tiny':0},'ratio_order_failure_count':0} for k in ['general_G','tridiagonal']}
cases=[]
def run(label,data,expected):
 env={'results':copy.deepcopy(data)};rejected=False
 try:exec(program,env)
 except AssertionError:rejected=True
 assert rejected==expected and env['results']['passed']==(not expected),label
 cases.append(label)
run('Clean counters accepted',base,False)
d=copy.deepcopy(base);d['general_G']['KL']['negative']=1;run('One near-bath KL contradiction rejected',d,True)
d=copy.deepcopy(base);d['tridiagonal']['other_divergence']={'negative':1};run('Contradiction in any sampled convex divergence rejected',d,True)
d=copy.deepcopy(base);d['tridiagonal']['ratio_order_failure_count']=1;run('One tridiagonal ratio-order contradiction rejected',d,True)
d=copy.deepcopy(base);d['general_G']['ratio_order_failure_count']=100;run('General G-only ratio failures correctly permitted',d,False)
d=copy.deepcopy(base);d['tridiagonal']['quartic']={'negative':0,'skipped_tiny':320};run('Declared tiny-value exclusions retained without pretending they passed',d,False)
out={'passed':True,'checks':cases,'count':len(cases),'sampler_source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'scope':__doc__.strip()}
Path('research/detailed_audit/entropy_audit_verdict_regressions.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
