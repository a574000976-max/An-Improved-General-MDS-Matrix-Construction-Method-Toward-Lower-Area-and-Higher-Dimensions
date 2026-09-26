import numpy as np
from itertools import combinations
from sage.all_cmdline import *

P = PolynomialRing(GF(2), "L")
K = FractionField(P)
L = K.gen()


def det(matrix):
    row = matrix.shape[0]
    for i in range(row):
        non_zero_row = np.argmax(matrix[i:, i]) + i
        if matrix[non_zero_row, i] == 0:
            return 0
        matrix[[i, non_zero_row]] = matrix[[non_zero_row, i]]
        for j in range(i + 1, row):
            if matrix[j, i] == 1:
                matrix[j] = np.bitwise_xor(matrix[j], matrix[i])

    return 1


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
    R = [L]
    lst = []
    for i in range(len(M.rows()) + 1)[1:]:
        ll = M.minors(i)
        for l in ll:
            if l != 0:
                F = list(l.factor())
                for f in F:
                    R.append(f[0])
            else:
                return "Not MDS"
    lst = list(set(R))
    return lst


def type3(l, m, k):
    E = identity_matrix(K, 4)
    E[l, m] = k
    return E


def type2(l, k):
    E = identity_matrix(K, 4)
    E[l, l] = k
    return E


def gen_matrix_from_path(path):
    init_matrix = identity_matrix(K, 4)
    for i in path:
        if type(i[0]) == tuple:
            init_matrix = type3(i[0][0], i[0][1], i[1]) * init_matrix
        else:
            init_matrix = type2(i[0], i[1]) * init_matrix
    return init_matrix


def check(poly, L, bit_size):
    poly_list = poly
    X = Matrix(GF(2), bit_size, bit_size, L)
    I = identity_matrix(GF(2), bit_size)

    matrix_list = []
    for p in poly_list:
        poly = P(p)
        mat = sum(c * (X**i) if i > 0 else c * I for i, c in enumerate(poly.list()))
        matrix_list.append(mat)

    det_list = []
    for idx, M in enumerate(matrix_list):
        det = M.determinant()
        if det == 0:
            return False, L
        det_list.append(det)
    if all(det_list):
        return True, L


def verify_binary_path(path, positions, bit_size):
    """Check every square block submatrix independently of factor collection."""
    F = GF(2)
    X = Matrix(F, bit_size, bit_size, positions)
    I = identity_matrix(F, bit_size)
    Z = zero_matrix(F, bit_size)

    def evaluate_polynomial(poly):
        value = zero_matrix(F, bit_size)
        for coefficient in reversed(P(poly).list()):
            value = value * X + coefficient * I
        return value

    blocks = [[I.copy() if i == j else Z.copy() for j in range(4)]
              for i in range(4)]
    for row_position, coefficient in path:
        coefficient = K(coefficient)
        A = (evaluate_polynomial(coefficient.numerator()) *
             evaluate_polynomial(coefficient.denominator()).inverse())
        if isinstance(row_position, tuple):
            i, j = row_position
            blocks[i] = [blocks[i][k] + A * blocks[j][k] for k in range(4)]
        else:
            i = row_position
            blocks[i] = [A * blocks[i][k] for k in range(4)]

    checked = 0
    for size in range(1, 5):
        for rows in combinations(range(4), size):
            for cols in combinations(range(4), size):
                submatrix = block_matrix([[blocks[i][j] for j in cols] for i in rows])
                assert submatrix.rank() == size * bit_size, (rows, cols)
                checked += 1
    return checked


if __name__ == "__main__":
    binary_matrix_4 = gen_binary_matrix([[4], [1, 4], [2], [3]])

    binary_matrix_8 = gen_binary_matrix([[2, 8], [1], [2], [3], [4], [5], [6], [7]])

    binary_matrix_16 = gen_binary_matrix(
        [[1, 16], [1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11], [12], [13], [14], [15]]
    )

    binary_matrix_32 = gen_binary_matrix(
        [
            [2, 32],
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
        ]
    )
    
    binary_matrix_64 = gen_binary_matrix(
        [
            [4, 64],
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
        ]
    )

    path_4dim = [
        ((2, 1), 1),
        ((3, 0), 1),
        ((0, 2), L**-1),
        ((1, 3), 1),
        (1, L**-1),
        ((2, 1), 1),
        ((3, 0), L**-1),
        ((0, 2), 1),
        ((1, 3), 1),
    ]

    m_4dim = gen_matrix_from_path(path_4dim)
    poly_list_4dim = get_minors(m_4dim)

    assert poly_list_4dim != "Not MDS", "4-dim Matrix is not MDS"

    result_4bit = check(poly_list_4dim, binary_matrix_4, 4)
    result_8bit = check(poly_list_4dim, binary_matrix_8, 8)
    result_16bit = check(poly_list_4dim, binary_matrix_16, 16)
    result_32bit = check(poly_list_4dim, binary_matrix_32, 32)
    result_64bit = check(poly_list_4dim, binary_matrix_64, 64)

    assert result_4bit[0], "4-bit matrix representation is not MDS"
    assert result_8bit[0], "8-bit matrix representation is not MDS"
    assert result_16bit[0], "16-bit matrix representation is not MDS"
    assert result_32bit[0], "32-bit matrix representation is not MDS"
    assert result_64bit[0], "64-bit matrix representation is not MDS"
    
    # Verify the one-XOR implementation of the new 8-bit inverse.
    X8 = Matrix(GF(2), 8, 8, binary_matrix_8)
    inverse8 = Matrix(GF(2), 8, 8,
                      gen_binary_matrix([[2], [3], [4], [5], [6], [7], [8], [1, 3]], 8))
    assert X8 * inverse8 == identity_matrix(GF(2), 8)
    assert inverse8 * X8 == identity_matrix(GF(2), 8)
    xor_inverse8 = sum(sum(int(value) for value in row) - 1 for row in inverse8.rows())
    assert xor_inverse8 == 1
    type3_count = sum(isinstance(position, tuple) for position, _ in path_4dim)
    inverse_count = sum(K(coefficient) == L**-1 for _, coefficient in path_4dim)
    assert type3_count == 8 and inverse_count == 3
    assert all(K(coefficient) in (K(1), L**-1) for _, coefficient in path_4dim)
    xor_count8 = type3_count * 8 + inverse_count * xor_inverse8
    assert xor_count8 == 67
    checked8 = verify_binary_path(path_4dim, binary_matrix_8, 8)
    assert checked8 == 69
    print(f"8-bit construction: {checked8}/69 square block submatrices invertible; "
          f"L^-1 uses {xor_inverse8} XOR; total = {xor_count8} XOR gates.")

    print("All matrices (4, 8, 16, 32, 64-bit) are MDS. Verification Passed!")