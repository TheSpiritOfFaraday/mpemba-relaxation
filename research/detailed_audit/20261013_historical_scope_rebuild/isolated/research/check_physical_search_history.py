"""Independently recompute saved discovery points, not optimizer trajectories.

Uses direct scipy.linalg.expm in physical time, never the discovery spectral
propagator or objective. A passing result is numerical correspondence only.
Historical seeds, evaluation counts, and wall times cannot be authenticated
by this check; no global optimum or continuous-box guarantee follows.
"""
from pathlib import Path
import hashlib
import itertools
import json
import math
import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
TOL = 2e-8
CORNERS = np.array(list(itertools.product([-1, 1], repeat=8)))
errors = []


def require(condition, message):
    if not condition:
        raise ValueError(message)


def compare(actual, recorded, label):
    error = abs(float(actual) - float(recorded))
    require(math.isfinite(error) and error < TOL, (label, error))
    errors.append(error)


def perturb(x, signs):
    x = np.asarray(x, dtype=float)
    origin = np.array([0.] * 6 + [1., 1.])
    return origin + (x - origin) * (1 + .05 * np.asarray(signs))


class Direct:
    def __init__(self, x):
        self.x = np.asarray(x, dtype=float)
        E = np.r_[0., np.cumsum(self.x[:3])]
        self.pi = self.thermal(E, 1.)
        self.p = [self.thermal(E, t) for t in self.x[-2:]]
        self.M = np.zeros((4, 4))
        for i, j in itertools.combinations(range(4), 2):
            h = x[3] if (i, j) == (0, 3) else x[4] if (i, j) == (1, 2) else x[5]
            saddle = max(E[i], E[j]) + h
            self.M[i, j] = math.exp(E[j] - saddle)
            self.M[j, i] = math.exp(E[i] - saddle)
        for j in range(4):
            self.M[j, j] = -sum(self.M[i, j] for i in range(4) if i != j)

    @staticmethod
    def thermal(E, T):
        q = np.exp(-E / T)
        return q / sum(q)

    def distance(self, t, j):
        p = expm(self.M * t) @ self.p[j]
        z = p / self.pi - 1
        return float(sum(self.pi * ((1 + z) * np.log1p(z) - z)))

    def hitting(self, delta, j):
        hi = 1.
        require(self.distance(0, j) > delta, 'Target must start below distance')
        while self.distance(hi, j) > delta:
            hi *= 2
            require(hi < 1e10, 'Failed to bracket hitting time')
        return brentq(lambda t: self.distance(t, j) - delta, 0, hi,
                      xtol=1e-10, rtol=1e-13)

    def score(self):
        initial = min(self.distance(0, j) for j in (0, 1))
        if initial <= .011:
            return 10 + .011 - initial
        return max(self.distance(.95 * self.hitting(d, 0), 1) / d - 1
                   for d in [.001, .0031622776601683794, .01])


def objective(x, profile):
    x = np.asarray(x)
    if sum(x[:3]) > 4.76 or x[7] < x[6] + .3:
        return 10 + max(sum(x[:3]) - 4.76, x[6] + .3 - x[7])
    model = Direct(x)
    rates = model.M[~np.eye(4, dtype=bool)]
    if profile == 'nominal':
        ratio, limit = max(rates) / min(rates), 850
    else:
        # Evaluate every activation at its box endpoints, independently of
        # the closed-form discovery bound. Positivity makes these monotone.
        low, high = perturb(x, -np.ones(8)), perturb(x, np.ones(8))
        low_rates = Direct(low).M[~np.eye(4, dtype=bool)]
        high_rates = Direct(high).M[~np.eye(4, dtype=bool)]
        ratio, limit = max(low_rates) / min(high_rates), 995
    if ratio > limit:
        return 10 + math.log(ratio / limit)
    value = model.score()
    if profile == 'robust_screen':
        for j in range(8):
            for sign in (-1, 1):
                s = np.zeros(8); s[j] = sign
                value = max(value, Direct(perturb(x, s)).score())
    return value


def main():
    names = ['physical_rate_search', 'physical_rate_search_robust_screen',
             'physical_rate_refined', 'physical_rate_adversarial_robustness']
    records = [json.loads((ROOT / f'research/{n}.json').read_text()) for n in names]
    histories = []
    for record, profile in zip(records[:2], ['nominal', 'robust_screen']):
        history = record['history']
        require(all(a['evaluation'] < b['evaluation'] and a['value'] > b['value']
                    for a, b in zip(history, history[1:])), 'History must strictly improve')
        require(history[-1]['evaluation'] <= record['evaluations'], 'Evaluation range')
        for i, entry in enumerate(history):
            compare(objective(entry['parameters'], profile), entry['value'], (profile, i))
        compare(objective(record['result']['parameters'], profile), record['score'], profile)
        histories.append({'profile': profile, 'saved_improvements_checked': len(history)})
        print('Checked', profile, len(history), 'saved improvements', flush=True)
    x = np.array(records[1]['result']['parameters'])
    active = set()
    stages = []
    for stage in records[2]['stages']:
        before = np.array([Direct(perturb(x, s)).score() for s in CORNERS])
        active.update(np.argsort(before)[-12:].tolist())
        require(set(stage['active_corner_indices']) == active, 'Active corner selection')
        y = np.array(stage['parameters'])
        after = np.array([Direct(perturb(y, s)).score() for s in CORNERS])
        compare(max(after), stage['all_corner_score'], 'All-corner score')
        compare(max(after[list(active)]), stage['score'], 'Active-corner score')
        accepted = bool(max(after) < max(before))
        require(accepted == stage['accepted'], 'Stage acceptance')
        if accepted:
            x = y
        stages.append({'stage': stage['stage'], 'before_worst': float(max(before)),
                       'after_worst': float(max(after)), 'accepted': accepted})
        active.update(np.argsort(after)[-12:].tolist())
        print('Checked refinement stage', stage['stage'], flush=True)
    require(x.tolist() == records[2]['result']['parameters'], 'Final candidate identity')
    compare(max(after), records[2]['all_corner_score'], 'Final all-corner score')
    for run in records[3]['runs']:
        y = perturb(x, run['perturbation'])
        require(np.max(np.abs(y - run['parameters'])) < 1e-14, 'Adversarial parameters')
        require(np.max(np.abs(run['perturbation'])) <= 1 and .001 <= run['threshold'] <= .01,
                'Adversarial point in domain')
        model = Direct(y)
        saving = 1 - model.hitting(run['threshold'], 1) / model.hitting(run['threshold'], 0)
        compare(saving, run['saving'], 'Adversarial saving')
    compare(min(r['saving'] for r in records[3]['runs']), records[3]['smallest_saving_found'],
            'Smallest saved saving')
    result = {'passed': True, 'method': 'Independent direct matrix exponential in physical time',
              'histories': histories, 'stages': stages, 'adversarial_points_checked': 3,
              'comparisons': len(errors), 'absolute_tolerance': TOL, 'maximum_error': max(errors),
              'optimizer_success_flags': [r['optimizer_success'] for r in records[3]['runs']],
              'limitations': ['Does not replay or authenticate historical optimizer trajectories, evaluation counts or wall times',
                             'All three adversarial optimizer runs report success=false; their points remain valid numerical candidates',
                             'No exhaustive optimum or continuous robustness proof is inferred'],
              'input_sha256': {f'research/{n}.json': hashlib.sha256((ROOT / f'research/{n}.json').read_bytes()).hexdigest() for n in names},
              'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (ROOT / 'research/detailed_audit/physical_search_history_checks.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
