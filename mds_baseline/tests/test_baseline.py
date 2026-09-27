import copy
import itertools
import json
from pathlib import Path
import random
import tempfile
import unittest
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import baseline as b


def reference_multiply(a, c, polynomial):
    # Carryless product followed by polynomial long division, independent of b's shift reduction.
    product = 0
    while c:
        low = c & -c
        product ^= a << (low.bit_length() - 1)
        c ^= low
    while product.bit_length() >= polynomial.bit_length():
        product ^= polynomial << (product.bit_length() - polynomial.bit_length())
    return product


class FieldAndConstructionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.f = b.Field()

    def test_all_field_products_and_inverses(self):
        for a in range(256):
            for c in range(256):
                self.assertEqual(self.f.mul[a][c], reference_multiply(a, c, 0x187))
            if a:
                self.assertEqual(self.f.mul[a][self.f.inv[a]], 1)

    def test_all_coefficient_maps_on_all_inputs(self):
        for a in range(256):
            rows = self.f.block(a)
            for x in range(256):
                expected = self.f.decode[reference_multiply(a, self.f.encode[x], 0x187)]
                self.assertEqual(b.binary_apply(rows, x), expected)

    def test_paper_basis_and_original_fixture(self):
        fixtures = json.loads((Path(__file__).parent / 'original_fixtures.json').read_text())
        inverse = self.f.inv[2]
        inverse_cube = self.f.mul[self.f.mul[inverse][inverse]][inverse]
        self.assertEqual(b.row_strings(self.f.block(inverse)), fixtures[3]['rows'])
        self.assertEqual(b.row_strings(self.f.block(inverse_cube)), fixtures[4]['rows'])

    def test_bad_domains_rejected(self):
        with self.assertRaises(ValueError):
            b.Field(4, 0x15, 'polynomial')
        with self.assertRaises(ValueError):
            b.cauchy(self.f, [0, 1], [1, 2])
        with self.assertRaises(ValueError):
            b.Field(4, 0x13, [1, 2, 4, 4])

    def test_normalization_and_exact_minors(self):
        f = b.Field(4, 0x13, 'polynomial')
        for normalized in [False, True]:
            a = b.cauchy(f, [0, 1, 2, 3], [4, 5, 6, 7], normalized)
            if normalized:
                self.assertEqual(a[0], [1] * 4)
                self.assertEqual([r[0] for r in a], [1] * 4)
            checked = 0
            for k in range(1, 5):
                for rows in itertools.combinations(range(4), k):
                    for cols in itertools.combinations(range(4), k):
                        minor = [[a[i][j] for j in cols] for i in rows]
                        self.assertNotEqual(b.field_det(f, minor), 0)
                        binary = b.expand(f, minor)
                        self.assertEqual(b.binary_rank(binary, k * 4), k * 4)
                        checked += 1
            result = b.verify_all_minors(f, a)
            self.assertTrue(result['mds'])
            self.assertEqual(result['checked'], checked)
            corrupted = copy.deepcopy(a)
            corrupted[0] = corrupted[1][:]
            self.assertFalse(b.verify_all_minors(f, corrupted)['mds'])

    def test_negative_circuit_and_lift_checks(self):
        rows = [3, 2]
        circuit = {'size': 2, 'xor_count': 1, 'input_permutation': [0, 1], 'gates': [[0, 1]]}
        b.verify_circuit(rows, circuit)
        for x in range(4):
            self.assertEqual(b.circuit_apply(circuit, x), b.binary_apply(rows, x))
        bad = copy.deepcopy(circuit)
        bad['gates'][0] = [1, 0]
        with self.assertRaises(ValueError):
            b.verify_circuit(rows, bad)
        bad = copy.deepcopy(circuit)
        bad['xor_count'] = 2
        with self.assertRaises(ValueError):
            b.verify_circuit(rows, bad)
        for target in [4, 8]:
            lifted = b.lift_rows(rows, 2, target)
            c = b.lift_circuit(circuit, 2, target)
            b.verify_circuit(lifted, c)
            b.verify_lift_structure(rows, lifted, 2, target)
            for x in range(1 << target):
                self.assertEqual(b.circuit_apply(c, x), b.binary_apply(lifted, x))
            lifted[0] ^= 1
            with self.assertRaises(ValueError):
                b.verify_lift_structure(rows, lifted, 2, target)


class OptimizerTests(unittest.TestCase):
    def test_11_original_coefficient_examples(self):
        fixtures = json.loads((Path(__file__).parent / 'original_fixtures.json').read_text())
        with tempfile.TemporaryDirectory() as td:
            for fixture in fixtures:
                with self.subTest(choice=fixture['choice']):
                    rows = b.parse_rows(fixture['rows'])
                    c = b.optimize(rows, Path(td) / str(fixture['choice']), seed=12345,
                                   iterations=8, windows=8, timeout=60)
                    self.assertFalse(c['timed_out'])
                    self.assertEqual(c['xor_count'], fixture['reference_count'])
                    b.verify_circuit(rows, c)

    def test_random_invertible_and_permutation_matrices(self):
        rng = random.Random(4321)
        with tempfile.TemporaryDirectory() as td:
            for trial in range(20):
                rows = [1 << i for i in range(8)]
                rng.shuffle(rows)
                for _ in range(trial * 3):
                    d, s = rng.sample(range(8), 2)
                    rows[d] ^= rows[s]
                c = b.optimize(rows, Path(td) / str(trial), seed=trial,
                               iterations=3, windows=12, timeout=30)
                b.verify_circuit(rows, c)
                for x in range(256):
                    self.assertEqual(b.circuit_apply(c, x), b.binary_apply(rows, x))

    def test_reproducible_fixed_seed(self):
        rng = random.Random(900)
        rows = [1 << i for i in range(8)]
        for _ in range(70):
            d, s = rng.sample(range(8), 2)
            rows[d] ^= rows[s]
        with tempfile.TemporaryDirectory() as td:
            a = b.optimize(rows, Path(td) / 'a', seed=121, iterations=5, windows=20)
            c = b.optimize(rows, Path(td) / 'b', seed=121, iterations=5, windows=20)
            self.assertEqual(a['gates'], c['gates'])
            self.assertEqual(a['input_permutation'], c['input_permutation'])

    def test_singular_input_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            with self.assertRaises(ValueError):
                b.optimize([1] * 8, td)


if __name__ == '__main__':
    unittest.main(verbosity=2)
