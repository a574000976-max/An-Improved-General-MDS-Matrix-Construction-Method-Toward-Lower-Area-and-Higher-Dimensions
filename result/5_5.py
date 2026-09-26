
import numpy as np
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
    E = identity_matrix(K, 5)
    E[l, m] = k
    return E


def type2(l, k):
    E = identity_matrix(K, 5)
    E[l, l] = k
    return E


def gen_matrix_from_path(path):
    init_matrix = identity_matrix(K, 5)
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


if __name__ == "__main__":
    binary_matrix_4 = gen_binary_matrix([[4], [1, 4], [2], [3]])

    binary_matrix_8 = gen_binary_matrix([[8], [1, 2], [2, 8], [3], [4], [5], [6], [7]])

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

    path_4bit = [
        ((1, 3), L**-1),
      ((4, 0), L),
      ((0, 1), 1),
      (0, L**-1),
      ((3, 4), 1),
      (4, L),
      ((2, 3), 1),
      ((1, 2), L**3),
      (2, L**-1),
      ((4, 0), L**-1),
      ((3, 1), 1),
      ((0, 3), 1),
      ((2, 4), 1),
      ((3, 2), 1),
      ((4, 0), 1),
    ]
    path_higher = [
        ((0, 1), 1),
        ((2, 3), 1),
        ((3, 0), L**-1),
        ((1, 2), 1),
        ((4, 3), 1),
        ((2, 4), L**-3),
        (2, L**-1),
        ((3, 2), 1),
        (3, L**-1),
        ((0, 1), L**-1),
        ((4, 0), 1),
        ((1, 3), 1),
        ((3, 4), 1),
        ((0, 1), 1),
    ]
    m_4bit = gen_matrix_from_path(path_4bit)
    m_higher = gen_matrix_from_path(path_higher)
    
    poly_list_4bit = get_minors(m_4bit)
    poly_list_higher = get_minors(m_higher)


    assert poly_list_4bit != "Not MDS", "Matrix is not MDS"
    assert poly_list_higher != "Not MDS", "Matrix is not MDS"

    result_4bit = check(poly_list_4bit, binary_matrix_4, 4)
    result_8bit = check(poly_list_higher, binary_matrix_8, 8)
    result_16bit = check(poly_list_higher, binary_matrix_16, 16)
    result_32bit = check(poly_list_higher, binary_matrix_32, 32)
    result_64bit = check(poly_list_higher, binary_matrix_64, 64)

    assert result_4bit[0], "4-bit matrix is not MDS"
    assert result_8bit[0], "8-bit matrix is not MDS"
    assert result_16bit[0], "16-bit matrix is not MDS"
    assert result_32bit[0], "32-bit matrix is not MDS"
    assert result_64bit[0], "64-bit matrix is not MDS"
    print("All matrices are MDS.")
