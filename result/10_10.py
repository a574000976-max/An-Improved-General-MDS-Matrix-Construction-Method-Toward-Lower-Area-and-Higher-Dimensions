from itertools import combinations
import argparse
from sage.all_cmdline import *

# ============================================================
# Symbolic polynomial ring for the path
# ============================================================

P = PolynomialRing(GF(2), "L")
K = FractionField(P)
L = K.gen()
DIM = 10


# ============================================================
# 10x10 path from the verified construction
# ============================================================

PATH = [
    ((0, 1), 1),
    ((2, 3), 1),
    ((4, 5), 1),
    ((6, 7), 1),
    (8, L**-1),
    ((8, 9), 1),
    (4, L),
    ((3, 4), 1),
    ((1, 2), L),
    ((5, 0), L),
    ((7, 8), 1),
    ((9, 6), L**-1),
    ((0, 7), L**-1),
    ((8, 3), L**-2),
    ((2, 9), L**-2),
    ((6, 5), L**-2),
    ((4, 1), L**-1),
    ((3, 6), L**-2),
    ((9, 4), L**-4),
    ((7, 2), L**-3),
    ((1, 0), L**-3),
    (5, L),
    ((5, 8), L**-1),
    ((8, 1), L**-4),
    ((2, 8), L**-2),
    ((4, 7), L**-2),
    ((1, 3), L**-3),
    ((3, 9), 1),
    ((7, 3), 1),
    ((6, 2), 1),
    ((0, 6), 1),
    ((2, 5), L**2),
    ((5, 7), L**4),
    ((3, 0), L**3),
    (0, L**-2),
    ((0, 5), 1),
    (9, L**-1),
    ((9, 0), 1),
    ((5, 4), L**-1),
    ((6, 9), 1),
    ((4, 6), 1),
    ((8, 4), 1),
    ((0, 8), L),
]


# ============================================================
# Fixed binary matrices for L
# ============================================================

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


def get_binary_matrices():
    # ------------------------------------------------------------
    # n = 16
    # Field/Ring:
    #   GF(2^16) / (x^16 + x^5 + x^3 + x^2 + 1)
    # L:
    #   fixed binary matrix exported from the Sage generator.
    # ------------------------------------------------------------
    binary_matrix_16 = gen_binary_matrix(
        [
            [16],
            [1],
            [2, 16],
            [3, 16],
            [4],
            [5, 16],
            [6],
            [7],
            [8],
            [9],
            [10],
            [11],
            [12],
            [13],
            [14],
            [15],
        ],
        ncols=16,
    )

    # ------------------------------------------------------------
    # n = 32
    # Field/Ring:
    #   GF(2^32) / (x^32 + x^15 + x^9 + x^7 + x^4 + x^3 + 1)
    # L:
    #   fixed binary matrix exported from the Sage generator.
    # ------------------------------------------------------------
    binary_matrix_32 = gen_binary_matrix(
        [
            [32],
            [1],
            [2],
            [3, 32],
            [4, 32],
            [5],
            [6],
            [7, 32],
            [8],
            [9, 32],
            [10],
            [11],
            [12],
            [13],
            [14],
            [15, 32],
            [16],
            [17],
            [18],
            [19],
            [20],
            [21],
            [22],
            [23],
            [24],
            [25],
            [26],
            [27],
            [28],
            [29],
            [30],
            [31],
        ],
        ncols=32,
    )

    # ------------------------------------------------------------
    # n = 64
    # Use the same 64-bit L matrix as the previous 9-dimensional
    # verification script.
    #
    # Field/Ring:
    #   GL(64, GF(2))
    # L:
    #   [[15,64], [1], [2], ..., [63]]
    # ------------------------------------------------------------
    binary_matrix_64 = gen_binary_matrix(
        [
            [15, 64],
            [1],
            [2],
            [3],
            [4],
            [5],
            [6],
            [7],
            [8],
            [9],
            [10],
            [11],
            [12],
            [13],
            [14],
            [15],
            [16],
            [17],
            [18],
            [19],
            [20],
            [21],
            [22],
            [23],
            [24],
            [25],
            [26],
            [27],
            [28],
            [29],
            [30],
            [31],
            [32],
            [33],
            [34],
            [35],
            [36],
            [37],
            [38],
            [39],
            [40],
            [41],
            [42],
            [43],
            [44],
            [45],
            [46],
            [47],
            [48],
            [49],
            [50],
            [51],
            [52],
            [53],
            [54],
            [55],
            [56],
            [57],
            [58],
            [59],
            [60],
            [61],
            [62],
            [63],
        ],
        ncols=64,
    )

    return {
        16: binary_matrix_16,
        32: binary_matrix_32,
        64: binary_matrix_64,
    }


FIELD_DESCRIPTIONS = {
    16: "GF(2^16) / (x^16 + x^5 + x^3 + x^2 + 1)",
    32: "GF(2^32) / (x^32 + x^15 + x^9 + x^7 + x^4 + x^3 + 1)",
    64: "GL(64, GF(2)), with L = [[15,64], [1], [2], ..., [63]]",
}


# ============================================================
# Type-2 / Type-3 matrices over symbolic field
# ============================================================

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
            coeff = item[1]
            M = type3(i, j, coeff) * M
        else:
            i = item[0]
            coeff = item[1]
            M = type2(i, coeff) * M

    return M


# ============================================================
# Minor factor extraction
# ============================================================

def get_minor_factors(M):
    factors = set()

    # The path contains negative powers of L.
    # Hence L itself must be invertible in each tested representation.
    factors.add(P.gen())

    for k in range(1, DIM + 1):
        print(f"Calculating {k}x{k} minors...")
        count = 0

        for rows in combinations(range(DIM), k):
            for cols in combinations(range(DIM), k):
                count += 1
                minor = M.matrix_from_rows_and_columns(rows, cols).determinant()

                if minor == 0:
                    print(f"Zero minor found. rows={rows}, cols={cols}")
                    return "Not MDS"

                numerator = P(minor.numerator())
                denominator = P(minor.denominator())

                if numerator == 0:
                    print(f"Zero numerator found. rows={rows}, cols={cols}")
                    return "Not MDS"

                for fac, _ in numerator.factor():
                    if fac != 1:
                        factors.add(P(fac))

                for fac, _ in denominator.factor():
                    if fac != 1:
                        factors.add(P(fac))

        print(f"  checked {count} minors.")

    return sorted(factors, key=lambda x: str(x))


# ============================================================
# Fixed binary matrix verification
# ============================================================

def eval_poly_as_matrix(poly, L_matrix, bit_size):
    """
    Substitute the fixed binary matrix L_matrix into a polynomial p(L).
    """
    X = Matrix(GF(2), bit_size, bit_size, L_matrix)
    I = identity_matrix(GF(2), bit_size)
    Z = zero_matrix(GF(2), bit_size, bit_size)

    result = Z
    poly_obj = P(poly)

    for i, coeff in enumerate(poly_obj.list()):
        if coeff:
            result += I if i == 0 else X ** i

    return result


def check_over_binary_matrix(poly_list, L_matrix, bit_size):
    """
    Check that every polynomial factor p(L_matrix) is nonsingular.
    """
    for idx, poly in enumerate(poly_list):
        M = eval_poly_as_matrix(poly, L_matrix, bit_size)

        if M.determinant() == 0:
            print(f"[{bit_size}-bit] failed at factor #{idx}: {poly}")
            return False

    return True


# ============================================================
# Main
# ============================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bits",
        type=str,
        default="16,32,64",
        help="Comma-separated bit sizes, for example: 16,32,64",
    )
    args = parser.parse_args()

    bit_sizes = [int(x.strip()) for x in args.bits.split(",") if x.strip()]

    binary_matrices = get_binary_matrices()

    for n in bit_sizes:
        if n not in binary_matrices:
            raise ValueError("Only 16, 32, and 64 are supported in this script.")

    print("=" * 80)
    print("Generating 10x10 symbolic matrix from path...")
    print("=" * 80)

    M = gen_matrix_from_path(PATH)

    print("Extracting polynomial factors from all square minors...")
    poly_list = get_minor_factors(M)

    assert poly_list != "Not MDS", "10x10 symbolic matrix is not MDS."

    print(f"\nTotal unique polynomial factors: {len(poly_list)}")

    print("\n" + "=" * 80)
    print("Verification over fixed binary L matrices")
    print("=" * 80)

    passed = []
    failed = []

    for n in bit_sizes:
        print("\n" + "-" * 80)
        print(f"n = {n}")
        print("-" * 80)
        print(f"Field/Ring: {FIELD_DESCRIPTIONS[n]}")

        ok = check_over_binary_matrix(
            poly_list,
            binary_matrices[n],
            n,
        )

        if ok:
            print(f"10x10 {n}-bit matrix is MDS.")
            passed.append(n)
        else:
            print(f"10x10 {n}-bit matrix is NOT MDS.")
            failed.append(n)

    print("\n" + "=" * 80)
    print("Summary for 10x10 path")
    print("=" * 80)
    print(f"Passed: {passed}")
    print(f"Failed: {failed}")