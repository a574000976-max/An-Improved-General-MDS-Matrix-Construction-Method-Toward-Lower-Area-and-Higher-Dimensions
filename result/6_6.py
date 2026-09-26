import numpy as np
from sage.all_cmdline import *

P = PolynomialRing(GF(2), "L")
K = FractionField(P)
L = K.gen()

DIM = 6


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


def get_minors(M):

    factors = set()

    factors.add(P.gen())

    for k in range(1, M.nrows() + 1):
        print(f"Calculating {k}x{k} minors...")
        minors = M.minors(k)

        for minor in minors:
            if minor == 0:
                return "Not MDS"

            numerator = P(minor.numerator())
            denominator = P(minor.denominator())

            if numerator == 0:
                return "Not MDS"

            for fac, _ in numerator.factor():
                if fac != 1:
                    factors.add(P(fac))

            for fac, _ in denominator.factor():
                if fac != 1:
                    factors.add(P(fac))

    return list(factors)


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


def eval_poly_as_matrix(poly, L_matrix, bit_size):
    X = Matrix(GF(2), bit_size, bit_size, L_matrix)
    I = identity_matrix(GF(2), bit_size)
    Z = zero_matrix(GF(2), bit_size, bit_size)

    poly_obj = P(poly)
    result = Z

    for i, coeff in enumerate(poly_obj.list()):
        if coeff:
            if i == 0:
                result += I
            else:
                result += X ** i

    return result


def check(poly_list, L_matrix, bit_size):

    for idx, poly in enumerate(poly_list):
        M = eval_poly_as_matrix(poly, L_matrix, bit_size)
        determinant = M.determinant()

        if determinant == 0:
            print(f"[{bit_size}-bit] Failed at factor #{idx}: {poly}")
            return False, L_matrix

    return True, L_matrix


if __name__ == "__main__":

    # 8-bit
    binary_matrix_8 = gen_binary_matrix(
        [
            [8],
            [1, 2],
            [2, 8],
            [3],
            [4],
            [5],
            [6],
            [7],
        ],
        ncols=8
    )

    # 16-bit
    binary_matrix_16 = gen_binary_matrix(
        [
            [2, 16],
            [1],
            [2, 5],
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
        ],
        ncols=16
    )

    # 32-bit
    binary_matrix_32 = gen_binary_matrix(
        [
            [11, 32],
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
        ],
        ncols=32
    )

    # 64-bit
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
        ncols=64
    )

    path = [
        ((1, 3), 1),
        ((4, 2), 1),
        ((0, 5), 1),
        ((2, 0), 1),
        ((5, 4), 1),

        (5, L**-1),
        (4, L),

        ((0, 1), L),
        ((3, 2), 1),
        ((1, 5), L**-2),
        ((4, 3), L**-1),
        ((3, 1), L**-2),
        ((2, 4), 1),
        ((5, 0), 1),
        ((0, 2), L),

        (2, L**-1),

        ((4, 5), 1),
        ((1, 0), 1),
        ((2, 3), 1),
    ]

    print("=" * 80)
    print("Generating 6x6 matrix from path...")
    print("=" * 80)

    M = gen_matrix_from_path(path)

    print("\nCalculating minors...")
    poly_list = get_minors(M)

    assert poly_list != "Not MDS", "Symbolic matrix is not MDS."

    print(f"\nTotal unique polynomial factors: {len(poly_list)}")
    print("First several factors:")
    for p in poly_list[:20]:
        print(f"  {p}")

    print("\n" + "=" * 80)
    print("Checking against 8-bit, 16-bit, 32-bit, and 64-bit configurations...")
    print("=" * 80)

    result_8bit = check(poly_list, binary_matrix_8, 8)
    result_16bit = check(poly_list, binary_matrix_16, 16)
    result_32bit = check(poly_list, binary_matrix_32, 32)
    result_64bit = check(poly_list, binary_matrix_64, 64)

    assert result_8bit[0], "8-bit matrix is not MDS"
    print("8-bit matrix is MDS.")

    assert result_16bit[0], "16-bit matrix is not MDS"
    print("16-bit matrix is MDS.")

    assert result_32bit[0], "32-bit matrix is not MDS"
    print("32-bit matrix is MDS.")

    assert result_64bit[0], "64-bit matrix is not MDS"
    print("64-bit matrix is MDS.")

    print("\n" + "=" * 80)
    print("All matrices are MDS.")
    print("=" * 80)