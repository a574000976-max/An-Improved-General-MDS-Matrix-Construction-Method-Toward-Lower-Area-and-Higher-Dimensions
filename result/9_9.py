from itertools import combinations
from sage.all_cmdline import *

P = PolynomialRing(GF(2), "L")
K = FractionField(P)
L = K.gen()
DIM = 9


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
    binary_matrix_16 = gen_binary_matrix(
        [[2, 16], [1], [2, 5], [3], [4], [5], [6], [7],
         [8], [9], [10], [11], [12], [13], [14], [15]],
        ncols=16,
    )

    binary_matrix_32 = gen_binary_matrix(
        [[11, 32], [1], [2], [3], [4], [5], [6], [7], [8], [9],
         [10], [11], [12], [13], [14], [15], [16], [17], [18],
         [19], [20], [21], [22], [23], [24], [25], [26], [27],
         [28], [29], [30], [31]],
        ncols=32,
    )

    binary_matrix_64 = gen_binary_matrix(
        [[15, 64], [1], [2], [3], [4], [5], [6], [7], [8], [9],
         [10], [11], [12], [13], [14], [15], [16], [17], [18],
         [19], [20], [21], [22], [23], [24], [25], [26], [27],
         [28], [29], [30], [31], [32], [33], [34], [35], [36],
         [37], [38], [39], [40], [41], [42], [43], [44], [45],
         [46], [47], [48], [49], [50], [51], [52], [53], [54],
         [55], [56], [57], [58], [59], [60], [61], [62], [63]],
        ncols=64,
    )

    return {
        16: binary_matrix_16,
        32: binary_matrix_32,
        64: binary_matrix_64,
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
            M = type3(i, j, item[1]) * M
        else:
            M = type2(item[0], item[1]) * M
    return M


def get_minor_factors(M):
    factors = set()
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
    X = Matrix(GF(2), bit_size, bit_size, L_matrix)
    I = identity_matrix(GF(2), bit_size)
    Z = zero_matrix(GF(2), bit_size, bit_size)

    result = Z
    poly_obj = P(poly)

    for i, coeff in enumerate(poly_obj.list()):
        if coeff:
            result += I if i == 0 else X ** i

    return result


def check(poly_list, L_matrix, bit_size):
    for idx, poly in enumerate(poly_list):
        M = eval_poly_as_matrix(poly, L_matrix, bit_size)
        if M.determinant() == 0:
            print(f"[{bit_size}-bit] failed at factor #{idx}: {poly}")
            return False
    return True


if __name__ == "__main__":
    path = [
        ((7, 6), 1),
        ((6, 3), 1),
        ((1, 7), 1),
        ((3, 1), 1),
        ((1, 2), L**-1),
        ((2, 0), L**-1),
        ((0, 6), 1),
        ((4, 1), 1),
        (6, L**-1),
        ((6, 4), L**2),
        ((5, 0), L**3),
        ((8, 5), L**3),
        (1, L**-1),
        ((1, 8), L**2),
        ((0, 1), 1),
        ((8, 6), 1),
        ((7, 0), L**-1),
        ((6, 7), L**-2),
        ((0, 3), L**-2),
        ((2, 5), L**-2),
        ((5, 6), L**-1),
        (7, L),
        ((7, 2), L**-1),
        (3, L**2),
        ((3, 5), 1),
        (5, L**2),
        ((5, 0), 1),
        ((4, 3), L),
        ((0, 8), L**3),
        ((2, 4), 1),
        ((3, 7), L),
        ((4, 0), L**2),
        ((8, 3), 1),
        ((1, 2), L**2),
        ((2, 5), L),
        ((0, 1), 1),
        ((7, 4), 1),
    ]

    print("=" * 80)
    print("Generating 9x9 matrix from path...")
    print("=" * 80)

    M = gen_matrix_from_path(path)

    print("Extracting polynomial factors from all square minors...")
    poly_list = get_minor_factors(M)
    assert poly_list != "Not MDS", "9x9 symbolic matrix is not MDS."

    print(f"Total unique polynomial factors: {len(poly_list)}")

    binary_matrices = get_binary_matrices()

    passed = []
    failed = []

    for bit_size in [16, 32, 64]:
        print(f"\nChecking 9x9 path against {bit_size}-bit configuration...")
        ok = check(poly_list, binary_matrices[bit_size], bit_size)

        if ok:
            print(f"9x9 {bit_size}-bit matrix is MDS.")
            passed.append(bit_size)
        else:
            print(f"9x9 {bit_size}-bit matrix is NOT MDS.")
            failed.append(bit_size)

    print("\n" + "=" * 80)
    print("Summary")
    print("=" * 80)
    print(f"Passed: {passed}")
    print(f"Failed: {failed}")

    if failed:
        print("Some bit-size configurations failed. This path is not universal for all tested L matrices.")
    else:
        print("All 9x9 checks passed.")