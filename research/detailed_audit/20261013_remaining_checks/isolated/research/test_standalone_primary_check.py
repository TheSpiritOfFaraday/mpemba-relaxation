"""Adversarial checks of the separate verifier, including isolated execution.

Run: python3 -I -S research/test_standalone_primary_check.py
These tests exercise rejection paths and small analytically known examples;
passing them is not a substitute for the mathematical argument in the report.
"""
from fractions import Fraction as F
from pathlib import Path
import json
import runpy
import shutil
import subprocess
import sys
import tempfile
import unittest

SOURCE = Path(__file__).with_name('standalone_primary_check.py')
C = runpy.run_path(str(SOURCE))


class AdversarialChecks(unittest.TestCase):
    def rejected(self, message, **changes):
        with self.assertRaisesRegex(C['VerificationError'], message):
            C['verify'](**changes)

    def test_isolated_no_dependencies_and_optimized_mode(self):
        with tempfile.TemporaryDirectory() as directory:
            copied = Path(directory)/'check.py'
            shutil.copyfile(SOURCE, copied)
            outputs = []
            for flags in [[], ['-O']]:
                result = subprocess.run([sys.executable, '-I', '-S', *flags, str(copied)],
                                        cwd=directory, capture_output=True, text=True, check=True)
                outputs.append(json.loads(result.stdout))
            self.assertEqual(outputs[0], outputs[1])
            self.assertTrue(outputs[0]['verified'])
            self.assertEqual(sorted(p.name for p in Path(directory).iterdir()), ['check.py'])

    def test_triplet_violation_is_rejected(self):
        rates = list(C['RATES'])
        rates[2] = F(1000)  # a14 large: the published 124 upper-branch test passes.
        self.rejected('published triplet predicate', rates=rates)

    def test_wrong_eigenvector_is_rejected(self):
        self.rejected('weighted eigenvector residual', prefix=(F(1), F(0), F(0)))

    def test_wrong_mode_interval_is_rejected(self):
        self.rejected('trial belongs to slow interval', trial=F('-59.5'))

    def test_reversed_temperature_order_is_rejected(self):
        self.rejected('finite cooling temperature ordering', hot_n=100, cold_n=1000)

    def test_zero_rate_is_rejected(self):
        rates = list(C['RATES'])
        rates[0] = F(0)
        self.rejected('strict rate positivity', rates=rates)

    def test_transposed_generator_is_rejected_by_weighted_symmetry(self):
        L = C['row_generator'](C['Q'], C['RATES'])
        transposed = [list(row) for row in zip(*L)]
        with self.assertRaisesRegex(C['VerificationError'], 'symmetric pencil'):
            C['positive_inertia'](transposed, C['Q'], F(-20))

    def test_inertia_on_known_nondiagonal_spectrum(self):
        # Two 2x2 symmetric blocks, each with eigenvalues -2 and 0.
        A = [[F(v) for v in row] for row in
             [[-1, 1, 0, 0], [1, -1, 0, 0], [0, 0, -1, 1], [0, 0, 1, -1]]]
        for shift, expected in [(F(-3), 4), (F('-0.5'), 2), (F(1), 0)]:
            count, _ = C['positive_inertia'](A, [F(1)]*4, shift)
            self.assertEqual(count, expected)

    def test_zero_pivot_fails_closed(self):
        with self.assertRaisesRegex(C['VerificationError'], 'zero LDL pivot'):
            C['positive_inertia']([[F(0)]*4 for _ in range(4)], [F(1)]*4, F(0))

    def test_known_root_and_signed_overlap(self):
        lo, hi = C['thermal_root'](F(1, 4), 2)
        self.assertLessEqual(lo, F(1, 2))
        self.assertLessEqual(F(1, 2), hi)
        weights = [F(1), F(1, 4), F(1, 9), F(1, 16)]
        vector = [F(-2), F(3), F(-4), F(5)]
        actual = (-2+F(3, 2)-F(4, 3)+F(5, 4))/(1+F(1, 2)+F(1, 3)+F(1, 4))
        lower, upper = C['approximate_overlap'](weights, vector, 2)
        self.assertLessEqual(lower, actual)
        self.assertLessEqual(actual, upper)


if __name__ == '__main__':
    unittest.main(verbosity=2)
