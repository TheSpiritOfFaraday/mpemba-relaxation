"""Build and run all four continuous-coverage obligations.

Requires a C++17 compiler and MPFR/GMP (Homebrew paths on this Mac).
An incomplete run exits nonzero and is never labeled certified.
"""
import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--radius',default='.05')
    parser.add_argument('--maxnodes',type=int,default=1000000)
    args=parser.parse_args()
    source=ROOT/'research/physical_box_spectral.cpp'
    binary=ROOT/'tmp/physical_box_certificate_runner'
    inputs=ROOT/'research/physical_rate_refined.json'
    x=json.loads(inputs.read_text())['result']['parameters']
    literals=re.search(r'const char\*params\[\]=\{([^}]+)\}',source.read_text()).group(1)
    encoded=re.findall(r'"([^"]+)"',literals)
    assert encoded==[str(v) for v in x], 'compiled decimal inputs differ from authoritative candidate'
    command=['clang++','-O2','-std=c++17','-ffp-contract=off','-I/opt/homebrew/include',str(source),
             '-L/opt/homebrew/lib','-lmpfr','-lgmp','-o',str(binary)]
    subprocess.run(command,check=True)
    name=f'spectral_box_{100*float(args.radius):g}pct'
    results=[]
    # Endpoint obligations precede the performance checks.
    for mode in [2,3,1,0]:
        log=ROOT/f'research/{name}_mode{mode}.log'
        with log.open('w') as stream:
            proc=subprocess.run([str(binary),'cover',args.radius,str(args.maxnodes),str(mode),'.1','2'],stdout=stream,stderr=subprocess.STDOUT)
        lines=log.read_text().splitlines()
        result=json.loads(lines[-1]) if lines and lines[-1].startswith('{') else dict(complete=False,error=lines[-1] if lines else 'missing output')
        result.update(exit_code=proc.returncode,log=str(log.relative_to(ROOT)),log_sha256=sha(log))
        results.append(result)
        print(json.dumps(result),flush=True)
    complete=all(r['complete'] and r['exit_code']==0 and r.get('pending')==0 for r in results)
    receipt=dict(complete=complete,relative_radius=args.radius,
                 parameters=x,source_sha256=sha(source),binary_sha256=sha(binary),inputs_sha256=sha(inputs),
                 compile_command=command,checks=results,
                 claims=dict(KL_distance_advantage_at_least=.10,normalized_distance_window=[.4,.8],cold_KL_floor=.001,
                             time_saving_at_least=.05,all_KL_thresholds=[.001,.01],threshold_crossing_window_tau=[.1,2]),
                 scope='four-state phenomenological Arrhenius rate model; no continuum or experimental claim',
                 mathematical_dependencies='See physical_robustness_certificate.tex; executable interval certificate, not Lean')
    (ROOT/f'research/{name}_certificate.json').write_text(json.dumps(receipt,indent=2)+'\n')
    raise SystemExit(0 if complete else 2)
