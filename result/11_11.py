from itertools import combinations
from sage.all_cmdline import *

P = PolynomialRing(GF(2), "L")
K = FractionField(P)
L = K.gen()

DIM = 11
BIT_SIZES = [16, 32, 64]


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
    # Field/Ring:
    #   GF(2^64) /
    #   (x^64 + x^33 + x^30 + x^26 + x^25 + x^24 + x^23
    #    + x^22 + x^21 + x^20 + x^18 + x^13 + x^12
    #    + x^11 + x^10 + x^7 + x^5 + x^4 + x^2 + x + 1)
    # ------------------------------------------------------------
    binary_matrix_64 = gen_binary_matrix(
        [
            [64],
            [1, 64],
            [2, 64],
            [3],
            [4, 64],
            [5, 64],
            [6],
            [7, 64],
            [8],
            [9],
            [10, 64],
            [11, 64],
            [12, 64],
            [13, 64],
            [14],
            [15],
            [16],
            [17],
            [18, 64],
            [19],
            [20, 64],
            [21, 64],
            [22, 64],
            [23, 64],
            [24, 64],
            [25, 64],
            [26, 64],
            [27],
            [28],
            [29],
            [30, 64],
            [31],
            [32],
            [33, 64],
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
    64: (
        "GF(2^64) / "
        "(x^64 + x^33 + x^30 + x^26 + x^25 + x^24 + x^23 "
        "+ x^22 + x^21 + x^20 + x^18 + x^13 + x^12 "
        "+ x^11 + x^10 + x^7 + x^5 + x^4 + x^2 + x + 1)"
    ),
}


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


def get_minor_factors(M):
    """
    Extract all unique polynomial factors from all square minors.
    If one symbolic minor is identically zero, return "Not MDS".
    """
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


def eval_poly_as_matrix(poly, L_matrix, bit_size):
    """
    Substitute the binary matrix L_matrix into a polynomial p(L).
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


if __name__ == "__main__":
    path = [
        ((10, 2), L**3),
        ((2, 1), L**7),
        ((9, 10), L**6),
        ((1, 9), L**2),
        ((0, 1), L**8),
        ((3, 0), L**2),
        ((10, 3), L**4),
        ((1, 4), L),
        ((4, 10), L**2),
        ((3, 5), 1),
        ((5, 6), L**-1),
        ((6, 1), 1),
        ((1, 7), L**2),
        ((7, 3), L**4),
        ((8, 1), 1),
        (1, L**4),
        ((10, 8), L**-4),
        ((1, 2), L**-6),
        ((2, 5), L**-9),
        ((5, 10), L**-1),
        ((3, 1), L**-2),
        ((9, 5), L**-1),
        (5, L**-3),
        ((10, 0), 1),
        ((0, 3), L**2),
        ((1, 9), L**2),
        ((9, 4), 1),
        ((8, 0), L**-1),
        ((0, 9), L**-3),
        ((4, 8), L**-4),
        (4, L**-1),
        ((3, 2), L**-2),
        ((2, 4), L**-4),
        ((4, 6), L**-1),
        ((6, 1), L**-1),
        ((5, 4), 1),
        ((1, 3), L**-2),
        ((4, 1), L**-6),
        ((3, 5), L**-4),
        ((8, 4), L**-5),
        ((7, 8), L**-5),
        ((1, 10), L**-1),
        ((10, 7), L**-1),
        ((8, 3), 1),
        (3, L**-2),
        ((5, 10), L**2),
        ((7, 0), 1),
        (0, L**-2),
        ((9, 1), L**2),
    ]

    print("=" * 80)
    print("Generating 11x11 matrix from path...")
    print("=" * 80)

    M = gen_matrix_from_path(path)

    print("Extracting polynomial factors from all square minors...")
    poly_list = get_minor_factors(M)

    assert poly_list != "Not MDS", "11x11 symbolic matrix is not MDS."

    print(f"\nTotal unique polynomial factors: {len(poly_list)}")

    binary_matrices = get_binary_matrices()

    passed = []
    failed = []

    print("\n" + "=" * 80)
    print("Verification over fixed binary L matrices")
    print("=" * 80)

    for bit_size in BIT_SIZES:
        print("\n" + "-" * 80)
        print(f"n = {bit_size}")
        print("-" * 80)
        print(f"Field/Ring: {FIELD_DESCRIPTIONS[bit_size]}")

        ok = check_over_binary_matrix(
            poly_list,
            binary_matrices[bit_size],
            bit_size,
        )

        if ok:
            print(f"11x11 {bit_size}-bit matrix is MDS.")
            passed.append(bit_size)
        else:
            print(f"11x11 {bit_size}-bit matrix is NOT MDS.")
            failed.append(bit_size)

    print("\n" + "=" * 80)
    print("Summary for 11x11 path")
    print("=" * 80)
    print(f"Passed: {passed}")
    print(f"Failed: {failed}")