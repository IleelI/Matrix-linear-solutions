def get_diagonal_matrix(matrix: [[]]):
    output = [[0 for _ in matrix] for _ in matrix]
    for (row_idx, row) in enumerate(matrix):
        output[row_idx][row_idx] = matrix[row_idx][row_idx]
    return output


def get_triangular_matrix(matrix: [[]], isUpper: bool):
    output = [[0 for _ in matrix] for _ in matrix]
    for (row_idx, row) in enumerate(matrix):
        for (col_idx, col) in enumerate(row):
            triangular_condition = col_idx > row_idx if isUpper else col_idx < row_idx
            output[row_idx][col_idx] = col if triangular_condition else 0
    return output


def add_matrices(matrix_1: [[]], matrix_2: [[]]):
    output = [[0 for _ in matrix_1] for _ in matrix_1]
    for (row_idx, row) in enumerate(matrix_1):
        for (col_idx, col) in enumerate(row):
            output[row_idx][col_idx] = matrix_1[row_idx][col_idx] + matrix_2[row_idx][col_idx]
    return output


def log_matrix(matrix):
    print("Start log")
    for row in matrix:
        print(row)
    print("End log")
