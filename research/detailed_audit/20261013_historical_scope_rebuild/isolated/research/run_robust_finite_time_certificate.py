"""Rebuild and replay all proof obligations for the complete +/-5% result.

No saved success flag is accepted in lieu of executing a coverage engine.
Requires clang++, MPFR and GMP (Homebrew paths on this Mac).
"""
import concurrent.futures
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import platform
import re
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'research/robust_finite_time_certificate'
SOURCES = {name: ROOT / f'research/{name}.cpp' for name in
           ['physical_temperature_reduction', 'physical_temperature_reduced_combined']}
OBLIGATIONS = {
  1: 'distance_window', 2: 'both_early_endpoints', 3: 'cold_late_endpoint',
  4: 'cold_temperature_concavity', 5: 'hot_temperature_monotonicity',
  6: 'hot_late_endpoint', 7: 'cold_low_hot_low_time', 8: 'cold_high_hot_low_time'}

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def verify_embedded_parameters(source_text, expected):
    """Compare exact decimal values; binary64 equality can hide different inputs."""
    match = re.search(r'const char\*params\[\]=\{([^}]+)\}', source_text)
    if match is None:
        raise ValueError('Missing embedded physical parameter array')
    actual = [Decimal(x) for x in re.findall(r'"([^"]+)"', match.group(1))]
    if len(actual) != 8 or any(not x.is_finite() for x in actual):
        raise ValueError('Expected eight finite embedded decimal parameters')
    if actual != expected:
        raise ValueError('Embedded decimal parameters differ from the input file')
    return actual

def main():
    OUT.mkdir(exist_ok=True)
    input_path = ROOT/'research/physical_rate_refined.json'
    input_before = digest(input_path)
    input_text = input_path.read_text()
    parameters = json.loads(input_text)['result']['parameters']
    exact_parameters = json.loads(input_text, parse_float=Decimal, parse_int=Decimal)['result']['parameters']
    binaries = {}
    builds = {}
    for name, source in SOURCES.items():
        source_before = digest(source)
        verify_embedded_parameters(source.read_text(), exact_parameters)
        binary = ROOT / f'tmp/final_{name}'
        command = ['clang++','-O2','-std=c++17','-ffp-contract=off',
                   '-I/opt/homebrew/include',str(source),'-L/opt/homebrew/lib',
                   '-lmpfr','-lgmp','-o',str(binary)]
        subprocess.run(command,check=True,cwd=ROOT)
        if digest(source) != source_before:
            raise RuntimeError('Source changed during compilation')
        binaries[name]=binary
        builds[name]={'source_sha256':digest(source),'binary_sha256':digest(binary),'command':command}
    manifest = {'complete':False,'parameters':parameters,'relative_radius':0.05,
                'parameter_coordinates':['g1','g2','g3','hA','hB','hX','Tc-1','Th-1'],
                'platform':platform.platform(),'compiler':subprocess.check_output(['clang++','--version'],text=True),
                'input_sha256':input_before,
                'parameter_decimal_strings':list(map(str, exact_parameters)),
                'runner_sha256':digest(Path(__file__)), 'builds':builds,'obligations':{}}
    receipt_path=OUT/'certificate.json'
    receipt_path.write_text(json.dumps(manifest,indent=2)+'\n')
    def run(mode):
        name='physical_temperature_reduction' if mode in (4,5) else 'physical_temperature_reduced_combined'
        command=[str(binaries[name]),'cover','.05','1000000',str(mode),'.1','2']
        logfile=OUT/f'{mode}_{OBLIGATIONS[mode]}.log'
        with logfile.open('w') as stream:
            result=subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT,cwd=ROOT)
        lines=[line for line in logfile.read_text().splitlines() if line.startswith('{')]
        record=json.loads(lines[-1]) if lines else {}
        assert result.returncode==0 and record.get('complete') is True, (mode,record)
        assert record['pending']==0 and record['mode']==mode and record['radius']>=.05
        assert record['start_u'][0]<=.1<=record['start_u'][1] and record['end_u']==[2,2]
        assert record['cells']==2*record['passed']-1
        assert digest(SOURCES[name])==builds[name]['source_sha256']
        assert digest(binaries[name])==builds[name]['binary_sha256']
        record.update(source=name,command=command,log_sha256=digest(logfile))
        print(f"Completed mode {mode}: {OBLIGATIONS[mode]}, {record['cells']} nodes",flush=True)
        return mode,record
    # Independent proof obligations share read-only sources and write separate logs.
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        for mode,record in pool.map(run,[4,5,7,8,1,2,3,6]):
            manifest['obligations'][str(mode)]=record
            receipt_path.write_text(json.dumps(manifest,indent=2)+'\n')
    assert set(manifest['obligations'])==set(map(str,OBLIGATIONS))
    if digest(input_path) != input_before:
        raise RuntimeError('Physical input changed during the proof run')
    manifest.update(complete=True,conclusion={
        'distance':'Dh <= 0.9 Dc and Dc >= 0.001 throughout u in [0.4,0.8]',
        'time':'th(delta) <= 0.95 tc(delta) for every delta in [0.001,0.01]',
        'domain':'entire simultaneous independent +/-5% box in the eight listed coordinates',
        'reduction':'Cold KL concavity makes the earliest cold threshold time occur at a temperature endpoint. Hot KL temperature derivative at threshold roots is nonpositive, making the lowest hot temperature worst. Modes 7 and 8 certify both remaining cases.',
        'scope':'Executable interval proof with analytic dependencies; not Lean or external peer review.'})
    receipt_path.write_text(json.dumps(manifest,indent=2)+'\n')
    print('All eight obligations completed:',receipt_path,flush=True)

if __name__=='__main__':
    main()
