"""Replay proof-supporting computations; a passing script is not an analytic proof.

Receipts preserve commands, source hashes, exit status, and complete logs. No
discovery optimizer is rerun. Invoke separately for core, extensions, physical.
The inventory deliberately keeps manuscript-proof review separate from replay.
"""
from pathlib import Path
import argparse
import datetime
import hashlib
import json
import os
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
GROUPS = {
    'core': [
        'standalone_primary_check.py', 'test_standalone_primary_check.py',
        'adversarial_certificate.py', 'adversarial_rational_mode_certificate.py',
        'independent_mode_audit.py', 'finite_time_certificate.py',
        'mechanism_certificate.py', 'full_triplet_amplitude_certificate.py',
        'extended_arithmetic_certificate.py', 'extended_mode_obstruction_certificate.py',
        'quarter_limit_certificate.py', 'quarter_limit_thermal_certificate.py',
    ],
    'extensions': [
        'ordered_conductance_certificate.py', 'energy_threshold_certificate.py',
        'classification_exact_certificate.py', 'degenerate_exact_certificate.py',
        'inverse_realizability_audit.py', 'test_exact_commensurate_classifier.py',
        'transient_only_chi_certificate.py', 'five_state_extension_certificate.py',
        'one_fast_mode_certificate.py', 'compact_endpoint_equality_certificate.py',
        'independent_rational_classifier_audit.py', 'independent_classifier_nearzero_audit.py',
        'high_dimension_amplitude_certificate.py', 'clone_rate_bound_certificate.py',
        'five_state_large_advantage_certificate.py', 'thermal_ratio_obstruction_certificate.py',
        'reset_trace_obstruction.py', 'reset_trace_independent_audit.py',
        'kl_contrast_gradient_obstruction.py',
    ],
    'physical': ['run_robust_finite_time_certificate.py'],
}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--group', choices=GROUPS, required=True)
    args = parser.parse_args()
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    out = ROOT / 'research/detailed_audit' / (stamp + '_' + args.group)
    out.mkdir(parents=True)
    receipt = {'group': args.group, 'started_utc': stamp, 'complete': False,
               'scope': 'Executable replay only; analytic argument review remains separate.',
               'python': sys.version, 'checks': []}
    path = out / 'receipt.json'
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    env.pop('PYTHONOPTIMIZE', None)
    for name in GROUPS[args.group]:
        script = ROOT / 'research' / name
        command = [sys.executable, str(script)]
        log = out / (script.stem + '.log')
        check = {'script': str(script.relative_to(ROOT)), 'command': command,
                 'source_sha256_before': sha(script), 'status': 'running',
                 'log': str(log.relative_to(ROOT))}
        receipt['checks'].append(check)
        path.write_text(json.dumps(receipt, indent=2) + '\n')
        print('START', name, flush=True)
        start = time.monotonic()
        with log.open('w') as stream:
            result = subprocess.run(command, cwd=ROOT, env=env, stdout=stream,
                                    stderr=subprocess.STDOUT)
        check.update(exit_code=result.returncode, seconds=time.monotonic()-start,
                     source_sha256_after=sha(script), log_sha256=sha(log))
        check['status'] = ('passed' if result.returncode == 0 and
                          check['source_sha256_before'] == check['source_sha256_after']
                          else 'failed')
        path.write_text(json.dumps(receipt, indent=2) + '\n')
        print(check['status'].upper(), name, flush=True)
    receipt.update(complete=True, all_passed=all(c['status']=='passed' for c in receipt['checks']))
    path.write_text(json.dumps(receipt, indent=2) + '\n')
    print(path, flush=True)
    return 0 if receipt['all_passed'] else 1

if __name__ == '__main__':
    sys.exit(main())
