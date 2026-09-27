#!/usr/bin/env python3
"""Reproducible Cauchy MDS baselines with the supplied OILL XOR optimizer.

Python 3.10+, C++17 compiler; no third-party Python packages are required.
Column zero is the least significant bit of an internal integer row.
External matrix text prints column zero first (leftmost).
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
OUR_COUNTS = {9: {16: 588, 32: 1062, 64: 2086},
              10: {16: 754, 32: 1473, 64: 2484},
              11: {16: 1028, 32: 1992, 64: 4302}}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def write_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    os.replace(temp, path)


def poly_mul_mod(a, b, polynomial, bits):
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & (1 << bits):
            a ^= polynomial
    return result


class Field:
    """Small binary extension field; field integers encode polynomial coefficients."""
    def __init__(self, bits=8, polynomial=0x187, basis='paper'):
        require(1 <= bits <= 8, 'Base field degree must be in 1..8')
        require(polynomial.bit_length() == bits + 1, 'Incorrect polynomial degree')
        self.bits, self.polynomial, self.q = bits, polynomial, 1 << bits
        self.mul = [[poly_mul_mod(a, b, polynomial, bits) for b in range(self.q)]
                    for a in range(self.q)]
        self.inv = [0] * self.q
        for a in range(1, self.q):
            require(1 in self.mul[a], 'Reducible defining polynomial')
            self.inv[a] = self.mul[a].index(1)
        if isinstance(basis, list):
            vectors = basis
        elif basis == 'polynomial':
            vectors = [1 << i for i in range(bits)]
        elif basis == 'paper':
            require(bits == 8 and polynomial == 0x187,
                    'The paper basis requires GF(2^8)/0x187')
            vectors = [1, 2] + [(1 << i) ^ (1 << (i - 1)) for i in range(2, 8)]
        else:
            raise ValueError('Unknown basis')
        require(len(vectors) == bits and all(0 < x < self.q for x in vectors), 'Bad basis')
        self.basis = vectors
        self.encode = []
        for coordinate in range(self.q):
            v = 0
            for j, basis_vector in enumerate(vectors):
                if coordinate >> j & 1:
                    v ^= basis_vector
            self.encode.append(v)
        require(len(set(self.encode)) == self.q, 'Dependent basis vectors')
        self.decode = {v: i for i, v in enumerate(self.encode)}
        self.blocks = {}

    def block(self, value):
        require(0 <= value < self.q, 'Field element out of range')
        if value not in self.blocks:
            rows = [0] * self.bits
            for j, v in enumerate(self.basis):
                out = self.decode[self.mul[value][v]]
                for i in range(self.bits):
                    rows[i] |= ((out >> i) & 1) << j
            self.blocks[value] = rows
        return self.blocks[value]


def cauchy(field, x, y, normalized=False):
    m = len(x)
    require(m and len(y) == m and len(set(x + y)) == 2 * m,
            'X and Y must have distinct, disjoint elements')
    require(all(0 <= v < field.q for v in x + y), 'Cauchy element out of range')
    a = [[field.inv[u ^ v] for v in y] for u in x]
    if normalized:
        a = [[field.mul[field.mul[a[i][j]][a[0][0]]]
              [field.inv[field.mul[a[i][0]][a[0][j]]]]
              for j in range(m)] for i in range(m)]
    return a


def expand(field, a):
    m = len(a)
    require(m and all(len(r) == m for r in a), 'Expected a square matrix')
    rows = [0] * (m * field.bits)
    for i in range(m):
        for j in range(m):
            for bit, row in enumerate(field.block(a[i][j])):
                rows[i * field.bits + bit] |= row << (j * field.bits)
    return rows


def row_strings(rows):
    n = len(rows)
    return [''.join(str(r >> j & 1) for j in range(n)) for r in rows]


def parse_rows(strings):
    n = len(strings)
    require(n > 0 and all(len(r) == n and set(r) <= {'0', '1'} for r in strings),
            'Malformed binary rows')
    return [sum((c == '1') << j for j, c in enumerate(r)) for r in strings]


def matrix_text(rows):
    return str(len(rows)) + '\n' + '\n'.join(row_strings(rows)) + '\n'


def binary_apply(rows, x):
    return sum(((r & x).bit_count() & 1) << i for i, r in enumerate(rows))


def circuit_apply(circuit, x):
    registers = [(x >> j) & 1 for j in circuit['input_permutation']]
    for dst, src in circuit['gates']:
        registers[dst] ^= registers[src]
    return sum(v << i for i, v in enumerate(registers))


def verify_circuit(rows, circuit):
    """Symbolic verification on all basis inputs simultaneously, not random testing."""
    n = len(rows)
    require(circuit['size'] == n, 'Circuit dimension mismatch')
    p, gates = circuit['input_permutation'], circuit['gates']
    require(sorted(p) == list(range(n)), 'Invalid input permutation')
    require(circuit['xor_count'] == len(gates), 'XOR count does not equal gate count')
    registers = [1 << j for j in p]
    for gate in gates:
        require(len(gate) == 2, 'Malformed gate')
        d, s = gate
        require(0 <= d < n and 0 <= s < n and d != s, 'Invalid gate index')
        registers[d] ^= registers[s]
    require(registers == rows, 'Circuit does not implement the matrix')
    return True


def binary_rank(rows, columns):
    rows = list(rows)
    rank = 0
    for col in range(columns):
        pivot = next((r for r in range(rank, len(rows)) if rows[r] >> col & 1), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for r in range(rank + 1, len(rows)):
            if rows[r] >> col & 1:
                rows[r] ^= rows[rank]
        rank += 1
    return rank


def field_det(field, a):
    a = [list(r) for r in a]
    k, determinant = len(a), 1
    for col in range(k):
        p = next((r for r in range(col, k) if a[r][col]), None)
        if p is None:
            return 0
        a[p], a[col] = a[col], a[p]
        pivot = a[col][col]
        determinant = field.mul[determinant][pivot]
        for r in range(col + 1, k):
            factor = field.mul[a[r][col]][field.inv[pivot]]
            for j in range(col + 1, k):
                a[r][j] ^= field.mul[factor][a[col][j]]
    return determinant


def lift_index(index, base_bits, target_bits, lane):
    return (index // base_bits) * target_bits + lane * base_bits + index % base_bits


def lift_rows(rows, base_bits, target_bits):
    require(target_bits >= base_bits and target_bits % base_bits == 0, 'Invalid lifted word size')
    require(len(rows) % base_bits == 0, 'Invalid base dimensions')
    result = [0] * (len(rows) // base_bits * target_bits)
    for i, row in enumerate(rows):
        for lane in range(target_bits // base_bits):
            out = 0
            remaining = row
            while remaining:
                low = remaining & -remaining
                j = low.bit_length() - 1
                out |= 1 << lift_index(j, base_bits, target_bits, lane)
                remaining ^= low
            result[lift_index(i, base_bits, target_bits, lane)] = out
    return result


def lift_circuit(circuit, base_bits, target_bits):
    require(target_bits >= base_bits and target_bits % base_bits == 0, 'Invalid lifted word size')
    n = circuit['size'] // base_bits * target_bits
    p, gates = [0] * n, []
    for lane in range(target_bits // base_bits):
        for i, j in enumerate(circuit['input_permutation']):
            p[lift_index(i, base_bits, target_bits, lane)] = lift_index(j, base_bits, target_bits, lane)
        for d, s in circuit['gates']:
            gates.append([lift_index(d, base_bits, target_bits, lane),
                          lift_index(s, base_bits, target_bits, lane)])
    return {'size': n, 'xor_count': len(gates), 'input_permutation': p, 'gates': gates,
            'method': 'independent copies of verified base circuit'}


def write_circuit_text(path, circuit):
    lines = ['# All indices are zero-based. Each ^= costs one two-input XOR gate.',
             '# Initialize ALL registers simultaneously before running the gates.',
             'z = [x[j] for j in ' + repr(circuit['input_permutation']) + ']']
    lines += [f'z[{d}] ^= z[{s}]' for d, s in circuit['gates']]
    lines += ['y = z.copy()', f'# XOR gates: {circuit["xor_count"]}']
    Path(path).write_text('\n'.join(lines) + '\n', encoding='utf-8')


def compile_tool(size=None, sanitize=False):
    compiler = shutil.which(os.environ.get('CXX', 'g++'))
    require(compiler is not None, 'Install a C++17 compiler (g++), or set CXX')
    build = ROOT / 'build'
    build.mkdir(exist_ok=True)
    if size is None:
        sources = [ROOT / 'oill/verify_minors.cpp']
        name, flags = 'verify_minors', []
    else:
        require(1 <= size <= 1024, 'Binary evaluator supports sizes 1..1024')
        sources = [ROOT / ('oill/' + f) for f in ['driver.cpp', 'strategy.cpp', 'reduce.cpp']]
        name, flags = f'evaluate_{size}', [f'-DSIZE={size}', '-pthread']
    flags = ['-std=c++17', '-O2'] + flags
    if sanitize:
        name += '_sanitized'
        flags += ['-fsanitize=address,undefined', '-fno-omit-frame-pointer', '-g']
    exe = build / name
    dependencies = sources + list((ROOT / 'oill').glob('*.h'))
    digest = hashlib.sha256((' '.join(flags) + compiler).encode())
    for source in dependencies:
        digest.update(source.read_bytes())
    stamp = exe.with_suffix('.sha256')
    if not exe.exists() or not stamp.exists() or stamp.read_text() != digest.hexdigest():
        subprocess.run([compiler, *flags, *map(str, sources), '-o', str(exe)], check=True)
        stamp.write_text(digest.hexdigest())
    return exe


def optimize(rows, folder, seed=1, iterations=4, windows=0, window_max=24, timeout=120,
             sanitize=False):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    matrix_file, output = folder / 'matrix.txt', folder / 'optimizer.json'
    matrix_file.write_text(matrix_text(rows))
    if output.exists():
        output.unlink()  # Never mistake an old result for this run's checkpoint.
    exe = compile_tool(len(rows), sanitize)
    start, timed_out = time.perf_counter(), False
    try:
        proc = subprocess.run([str(exe), str(matrix_file), str(output), str(seed),
                               str(iterations), str(windows), str(window_max)],
                              capture_output=True, text=True, timeout=timeout)
        require(proc.returncode == 0, f'OILL failed: {proc.stderr}')
    except subprocess.TimeoutExpired:
        timed_out = True
    require(output.exists(), 'Optimizer ended before producing a valid checkpoint')
    circuit = json.loads(output.read_text())
    verify_circuit(rows, circuit)
    circuit['timed_out'] = timed_out
    circuit['wall_seconds'] = time.perf_counter() - start
    circuit['completed_requested_iterations'] = not timed_out
    write_json(output, circuit)
    return circuit


def verify_all_minors(field, a):
    exe = compile_tool()
    payload = f'{len(a)} {field.bits} {field.polynomial}\n' + '\n'.join(' '.join(map(str, r)) for r in a)
    start = time.perf_counter()
    proc = subprocess.run([str(exe)], input=payload, capture_output=True, text=True, check=True)
    result = json.loads(proc.stdout)
    result['seconds'] = time.perf_counter() - start
    if result['mds']:
        require(result['checked'] == math.comb(2 * len(a), len(a)) - 1, 'Incorrect minor coverage')
    return result


def verify_field_binary(field, a, rows, samples=32, seed=0):
    rng = random.Random(seed)
    for _ in range(samples):
        coordinates = [rng.randrange(field.q) for _ in a]
        values = [field.encode[x] for x in coordinates]
        output = []
        for row in a:
            v = 0
            for coefficient, x in zip(row, values):
                v ^= field.mul[coefficient][x]
            output.append(field.decode[v])
        x = sum(v << (i * field.bits) for i, v in enumerate(coordinates))
        y = sum(v << (i * field.bits) for i, v in enumerate(output))
        require(binary_apply(rows, x) == y, 'Field and binary implementations disagree')


def verify_lift_structure(base, lifted, base_bits, target_bits):
    """Compare each lifted coefficient block to I_k tensor the base block."""
    m, mask = len(base) // base_bits, (1 << base_bits) - 1
    for i in range(m):
        for lane in range(target_bits // base_bits):
            for bit in range(base_bits):
                row = lifted[i * target_bits + lane * base_bits + bit]
                for j in range(m):
                    block = (row >> (j * target_bits)) & ((1 << target_bits) - 1)
                    expected = ((base[i * base_bits + bit] >> (j * base_bits)) & mask) << (lane * base_bits)
                    require(block == expected, 'Incorrect lifted coefficient block')


def verify_artifact(path, exhaustive=True):
    path = Path(path)
    artifact = json.loads(path.read_text())
    require(artifact['schema'] == 'cauchy-oill-baseline-v1', 'Unknown artifact schema')
    metadata = artifact['field']
    field = Field(metadata['bits'], metadata['polynomial'], metadata['basis'])
    a = cauchy(field, artifact['x'], artifact['y'], artifact['normalized'])
    require(artifact['dimension'] == len(a), 'Incorrect matrix dimension metadata')
    require(a == artifact['field_matrix'], 'Cauchy construction does not match saved matrix')
    rows = expand(field, a)
    require(rows == parse_rows(artifact['binary_rows']), 'Binary expansion mismatch')
    verify_field_binary(field, a, rows)
    verify_circuit(rows, artifact['circuit'])
    result = {'base_circuit_exact': True, 'field_binary_checks': 32,
              'mds_certificate': 'Cauchy determinant formula with distinct disjoint X,Y; nonzero row/column scaling'}
    if exhaustive:
        result['base_minors'] = verify_all_minors(field, a)
        require(result['base_minors']['mds'], 'Matrix is not MDS')
    lifted_results = {}
    for n_str, entry in artifact['lifts'].items():
        n = int(n_str)
        saved = json.loads((path.parent / entry['file']).read_text())
        require(saved['word_bits'] == n and saved['dimension'] == len(a), 'Incorrect lift metadata')
        lifted = parse_rows(saved['binary_rows'])
        verify_lift_structure(rows, lifted, field.bits, n)
        verify_circuit(lifted, saved['circuit'])
        expected_count = artifact['circuit']['xor_count'] * (n // field.bits)
        require(saved['circuit']['xor_count'] == expected_count == entry['xor_count'], 'Lift cost mismatch')
        lifted_results[n_str] = {'xor_count': expected_count, 'circuit_exact': True,
                                'block_diagonal_structure_exact': True,
                                'mds': 'inherited from verified base; each minor is permutation-equivalent to copies of its base minor'}
    result['lifts'] = lifted_results
    return result


def run_search(args):
    require(args.samples > 0 and args.iterations >= 0 and args.windows >= 0, 'Invalid search budget')
    require(args.timeout > 0 and args.window_max >= 3, 'Invalid optimizer limits')
    dims = list(dict.fromkeys(args.dims))
    targets = list(dict.fromkeys(args.targets))
    require(all(2 <= m <= 12 for m in dims), 'Supported dimensions: 2..12')
    require(all(n >= 8 and n % 8 == 0 and n <= 64 for n in targets), 'Targets must be multiples of 8 up to 64')
    out = Path(args.out).resolve()
    require(not out.exists() or not any(out.iterdir()), 'Use a new or empty output directory')
    out.mkdir(parents=True, exist_ok=True)
    field = Field(8, args.polynomial, args.basis)
    config = vars(args).copy()
    config['command'] = ' '.join(sys.argv)
    config['python'] = sys.version
    config['compiler'] = subprocess.run([os.environ.get('CXX', 'g++'), '--version'], capture_output=True, text=True, check=True).stdout.splitlines()[0]
    config['normalization_candidates_per_draw'] = 2 if args.normalization == 'both' else 1
    config['metric'] = 'verified in-place two-input XOR circuit count; free bit permutations; upper bound, not minimum'
    write_json(out / 'config.json', config)
    summary = []
    for m in dims:
        started = time.perf_counter()
        rng = random.Random(args.seed + 1009 * m)
        folder = out / f'm{m}'
        folder.mkdir()
        records, best = [], None
        variants = [False, True] if args.normalization == 'both' else [args.normalization == 'normalized']
        for sample in range(args.samples):
            elements = rng.sample(range(256), 2 * m)
            x, y = elements[:m], elements[m:]
            for normalized in variants:
                a = cauchy(field, x, y, normalized)
                rows = expand(field, a)
                seed = (args.seed + m * 1000003 + sample * 2 + int(normalized)) % (2 ** 32)
                name = f'candidate_{sample:05d}_' + ('normalized' if normalized else 'raw')
                circuit = optimize(rows, folder / name, seed, args.iterations, args.windows, args.window_max, args.timeout)
                record = {'sample': sample, 'normalized': normalized, 'seed': seed,
                          'xor_count': circuit['xor_count'], 'naive_xor_count': sum(r.bit_count() - 1 for r in rows),
                          'timed_out': circuit['timed_out'], 'stage': circuit['stage'],
                          'wall_seconds': circuit['wall_seconds'], 'x': x, 'y': y}
                records.append(record)
                if best is None or circuit['xor_count'] < best['circuit']['xor_count']:
                    best = {'schema': 'cauchy-oill-baseline-v1', 'dimension': m,
                            'field': {'bits': field.bits, 'polynomial': field.polynomial, 'basis': field.basis},
                            'x': x, 'y': y, 'normalized': normalized, 'sample': sample,
                            'field_matrix': a, 'binary_rows': row_strings(rows), 'circuit': circuit, 'lifts': {}}
                    write_json(folder / 'best_checkpoint.json', best)
                write_json(folder / 'search_log.json', records)
                print(f'm={m} sample={sample+1}/{args.samples} normalized={normalized} XOR={circuit["xor_count"]} best={best["circuit"]["xor_count"]} timeout={circuit["timed_out"]}', flush=True)
        base = parse_rows(best['binary_rows'])
        for n in targets:
            rows = lift_rows(base, 8, n)
            circuit = lift_circuit(best['circuit'], 8, n)
            verify_circuit(rows, circuit)
            filename = f'best_n{n}.json'
            write_json(folder / filename, {'dimension': m, 'word_bits': n, 'binary_rows': row_strings(rows), 'circuit': circuit})
            write_circuit_text(folder / f'circuit_n{n}.txt', circuit)
            best['lifts'][str(n)] = {'file': filename, 'xor_count': circuit['xor_count']}
            ours = OUR_COUNTS.get(m, {}).get(n)
            summary.append({'m': m, 'n': n, 'base_xor': best['circuit']['xor_count'],
                            'baseline_xor': circuit['xor_count'], 'paper_xor': ours,
                            'reduction_percent': round(100 * (circuit['xor_count'] - ours) / circuit['xor_count'], 4) if ours is not None else None})
        write_json(folder / 'best.json', best)
        write_circuit_text(folder / 'circuit_n8.txt', best['circuit'])
        verification = verify_artifact(folder / 'best.json', exhaustive=True)
        verification['search_and_verification_seconds'] = time.perf_counter() - started
        verification['candidate_count'] = len(records)
        verification['timed_out_candidates'] = sum(r['timed_out'] for r in records)
        write_json(folder / 'verification.json', verification)
        print(f'm={m}: all {verification["base_minors"]["checked"]} base minors verified; lifted circuits verified', flush=True)
    with (out / 'comparison.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)
    lines = ['# Cauchy + OILL + subfield baseline', '',
             'Counts describe verified implementations found within the recorded budget, not minima.',
             'Paper counts are supplied reference values; this run does not independently re-evaluate the paper matrices.',
             '', '| m | n | baseline XOR | paper XOR | reduction (%) |', '|---:|---:|---:|---:|---:|']
    for r in summary:
        lines.append(f'| {r["m"]} | {r["n"]} | {r["baseline_xor"]} | {r["paper_xor"]} | {r["reduction_percent"]} |')
    (out / 'comparison.md').write_text('\n'.join(lines) + '\n')
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    p = sub.add_parser('run', help='Generate, optimize, lift, and fully verify Cauchy baselines')
    p.add_argument('--dims', nargs='+', type=int, default=[9, 10, 11])
    p.add_argument('--targets', nargs='+', type=int, default=[16, 32, 64])
    p.add_argument('--samples', type=int, default=32, help='Number of X,Y draws per dimension')
    p.add_argument('--normalization', choices=['raw', 'normalized', 'both'], default='both')
    p.add_argument('--seed', type=int, default=20260927)
    p.add_argument('--iterations', type=int, default=4)
    p.add_argument('--windows', type=int, default=0, help='Bounded equivalent-sequence attempts per iteration')
    p.add_argument('--window-max', type=int, default=24)
    p.add_argument('--timeout', type=float, default=120, help='Seconds per candidate; validated checkpoints survive timeout')
    p.add_argument('--polynomial', type=lambda x: int(x, 0), default=0x187)
    p.add_argument('--basis', choices=['paper', 'polynomial'], default='paper')
    p.add_argument('--out', default='results')
    p = sub.add_parser('verify', help='Independently reconstruct and verify a saved best.json')
    p.add_argument('artifact')
    p = sub.add_parser('evaluate', help='Evaluate any nonsingular square binary matrix')
    p.add_argument('matrix', help='Text file: dimension then 0/1 rows, column zero first')
    p.add_argument('--out', default='evaluation')
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--iterations', type=int, default=4)
    p.add_argument('--windows', type=int, default=0)
    p.add_argument('--window-max', type=int, default=24)
    p.add_argument('--timeout', type=float, default=120)
    args = parser.parse_args()
    if args.action == 'run':
        run_search(args)
    elif args.action == 'verify':
        print(json.dumps(verify_artifact(args.artifact), indent=2, ensure_ascii=False))
    elif args.action == 'evaluate':
        tokens = Path(args.matrix).read_text().split()
        require(tokens and int(tokens[0]) == len(tokens) - 1, 'Incorrect matrix file size')
        rows = parse_rows(tokens[1:])
        require(args.iterations >= 0 and args.windows >= 0 and args.window_max >= 3 and args.timeout > 0, 'Invalid optimizer limits')
        circuit = optimize(rows, Path(args.out).resolve(), args.seed, args.iterations, args.windows, args.window_max, args.timeout)
        write_circuit_text(Path(args.out) / 'circuit.txt', circuit)
        print(json.dumps({'xor_count': circuit['xor_count'], 'verified': True,
                          'timed_out': circuit['timed_out']}, indent=2))


if __name__ == '__main__':
    main()
