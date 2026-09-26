#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Verify 8x8 MDS paths for 8-bit and 16/32/64-bit configurations.

- 8-bit uses path_8bit and L matrix for F2^8 / 0x1c3 Exp. (10):
  [[8], [1,3], [1,2,3], [3], [4], [5], [6], [7]]

- 16/32/64-bit use path_higher and Exp. (11) GL(n, F2) matrices:
  n=16: [[2,16], [1], [2,5], [3], ..., [15]]
  n=32: [[11,32], [1], [2], ..., [31]]
  n=64: [[15,64], [1], [2], ..., [63]]

Run:
  sage -python verify_8x8_mds_8_16_32_64.py
"""

from itertools import combinations
from sage.all_cmdline import *  # noqa: F401,F403

P = PolynomialRing(GF(2), "L")
K = FractionField(P)
L = K.gen()
DIM = 8


def gen_binary_matrix(positions, ncols=None):
    if ncols is None:
        ncols = max((max(row) if row else 0) for row in positions)

    matrix = []
    for row in positions:
        vec = [0] * ncols
        for j in row:
            vec[j - 1] = 1
        matrix.append(vec)
    return matrix


def type3(i, j, coeff):
    E = identity_matrix(K, DIM)
    E[i, j] = coeff
    return E


def type2(i, coeff):
    E = identity_matrix(K, DIM)
    E[i, i] = coeff
    return E


def gen_matrix_from_path(path):
    M = identity_matrix(K, DIM)
    for item in path:
        if isinstance(item[0], tuple):
            i, j = item[0]
            M = type3(i, j, item[1]) * M
        else:
            i = item[0]
            M = type2(i, item[1]) * M
    return M


def get_minor_factors(M):
    """
    Extract all nonzero square-minor polynomial factors.
    If a symbolic minor is identically zero, the path is not MDS symbolically.
    """
    factors = set()
    factors.add(P.gen())  # L must be invertible because negative powers may appear.

    for k in range(1, DIM + 1):
        print(f"Calculating {k}x{k} minors...")
        checked = 0
        for rows in combinations(range(DIM), k):
            for cols in combinations(range(DIM), k):
                checked += 1
                minor = M.matrix_from_rows_and_columns(rows, cols).determinant()

                if minor == 0:
                    print(f"Zero symbolic minor found: rows={rows}, cols={cols}")
                    return "Not MDS"

                numerator = P(minor.numerator())
                denominator = P(minor.denominator())

                if numerator == 0:
                    print(f"Zero numerator found: rows={rows}, cols={cols}")
                    return "Not MDS"

                for fac, _ in numerator.factor():
                    if fac != 1:
                        factors.add(P(fac))

                for fac, _ in denominator.factor():
                    if fac != 1:
                        factors.add(P(fac))

        print(f"  checked {checked} minors.")

    return sorted(factors, key=lambda x: str(x))


def eval_poly_as_matrix(poly, L_matrix, bit_size):
    """Evaluate a polynomial p(L) as p(L_matrix) over GF(2)."""
    X = Matrix(GF(2), bit_size, bit_size, L_matrix)
    I = identity_matrix(GF(2), bit_size)
    Z = zero_matrix(GF(2), bit_size, bit_size)

    poly_obj = P(poly)
    result = Z
    for i, coeff in enumerate(poly_obj.list()):
        if coeff:
            result += I if i == 0 else X ** i
    return result


def check(poly_list, L_matrix, bit_size, label):
    print(f"\nChecking {label} ({bit_size}-bit)...")
    for idx, poly in enumerate(poly_list):
        M = eval_poly_as_matrix(poly, L_matrix, bit_size)
        if M.determinant() == 0:
            print(f"[{label}] failed at factor #{idx}: {poly}")
            return False
    print(f"{label} is MDS.")
    return True


if __name__ == "__main__":
    # 8-bit: F2^8 / 0x1c3 Exp. (10)
    binary_matrix_8 = gen_binary_matrix(
        [[8], [1, 3], [1, 2, 3], [3], [4], [5], [6], [7]],
        ncols=8,
    )

    # 16-bit: GL(16,F2) Exp. (11)
    binary_matrix_16 = gen_binary_matrix(
        [[2, 16], [1], [2, 5], [3], [4], [5], [6], [7],
         [8], [9], [10], [11], [12], [13], [14], [15]],
        ncols=16,
    )

    # 32-bit: GL(32,F2) Exp. (11)
    binary_matrix_32 = gen_binary_matrix(
        [[21, 32]] + [[i] for i in range(1, 32)],
        ncols=32,
    )

    # 64-bit: GL(64,F2) Exp. (11)
    binary_matrix_64 = gen_binary_matrix(
        [[15, 64]] + [[i] for i in range(1, 64)],
        ncols=64,
    )

    # Path for 8-bit result.
    path_8bit = [
        ((1, 0), 1),
        ((3, 2), 1),
        ((5, 4), 1),
        ((7, 6), 1),
        ((0, 3), L**3),
        (0, L),
        ((2, 5), L**3),
        (2, L),
        ((4, 7), L**3),
        (4, L),
        ((6, 1), L**3),
        (6, L),
        ((3, 4), 1),
        ((5, 6), 1),
        ((7, 0), 1),
        ((1, 2), 1),
        ((4, 1), 1),
        ((6, 3), 1),
        ((0, 5), 1),
        ((2, 7), 1),
        ((1, 6), L**-3),
        ((3, 0), L**-3),
        ((5, 2), L**-3),
        ((7, 4), L**-3),
        ((6, 7), 1),
        ((0, 1), 1),
        ((2, 3), 1),
        ((4, 5), 1),
        ((1, 2), L**-6),
        ((5, 6), L**-6),
    ]

    # Path for 16/32/64-bit results.
    path_higher = [
        ((1, 0), 1),
        ((3, 2), 1),
        ((5, 4), 1),
        ((7, 6), 1),
        ((0, 3), L),
        ((2, 5), L),
        ((4, 7), L),
        ((6, 1), L),
        ((3, 4), 1),
        ((5, 6), 1),
        ((7, 0), 1),
        ((1, 2), 1),
        ((4, 1), 1),
        ((6, 3), 1),
        ((0, 5), 1),
        ((2, 7), 1),
        (1, L**2),
        (6, L**-1),
        ((1, 6), 1),
        (3, L),
        ((3, 0), L**-2),
        (5, L**2),
        (2, L**-1),
        ((5, 2), 1),
        (7, L),
        ((7, 4), L**-2),
        ((6, 7), 1),
        ((0, 1), 1),
        ((2, 3), 1),
        ((4, 5), 1),
        ((1, 2), 1),
        ((5, 6), 1),
    ]

    print("=" * 80)
    print("Generating 8x8 symbolic matrix for 8-bit path...")
    print("=" * 80)
    M_8 = gen_matrix_from_path(path_8bit)
    poly_list_8 = get_minor_factors(M_8)
    assert poly_list_8 != "Not MDS", "8-bit symbolic path is not MDS."
    print(f"8-bit path unique factors: {len(poly_list_8)}")

    print("\n" + "=" * 80)
    print("Generating 8x8 symbolic matrix for 16/32/64-bit path...")
    print("=" * 80)
    M_higher = gen_matrix_from_path(path_higher)
    poly_list_higher = get_minor_factors(M_higher)
    assert poly_list_higher != "Not MDS", "16/32/64-bit symbolic path is not MDS."
    print(f"16/32/64-bit path unique factors: {len(poly_list_higher)}")

    ok_8 = check(poly_list_8, binary_matrix_8, 8, "8-bit Exp.(10)")
    ok_16 = check(poly_list_higher, binary_matrix_16, 16, "16-bit Exp.(11)")
    ok_32 = check(poly_list_higher, binary_matrix_32, 32, "32-bit Exp.(11)")
    ok_64 = check(poly_list_higher, binary_matrix_64, 64, "64-bit Exp.(11)")

    assert ok_8, "8-bit matrix is not MDS"
    assert ok_16, "16-bit matrix is not MDS"
    assert ok_32, "32-bit matrix is not MDS"
    assert ok_64, "64-bit matrix is not MDS"

    print("\n" + "=" * 80)
    print("All requested 8x8 matrices are MDS.")
    print("=" * 80)
