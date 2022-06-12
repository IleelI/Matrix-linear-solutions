import time
import constants
import matplotlib.pyplot as plt
from copy import deepcopy
from generators import *
from matrix_utils import *
from utils import *


def jacobi_method(A, b, max_convergence, logStats=True):
    x = [0.0] * len(A)
    prev_x = [1.0] * len(A)
    D = get_diagonal_matrix(A)
    UL = add_matrices(get_triangular_matrix(A, True), get_triangular_matrix(A, False))
    iterations = 0
    start_time = time.time()
    while True:
        iterations += 1
        for row_idx, rows in enumerate(UL):
            row_sum = sum(row_val * x_val for row_val, x_val in zip(rows, prev_x))
            x[row_idx] = (b[row_idx] - row_sum) / D[row_idx][row_idx]
        for idx, _ in enumerate(x):
            prev_x[idx] = x[idx]
        normalised_residual = normalise(get_residual(A, b, x))
        if normalised_residual <= max_convergence:
            residual = normalised_residual
            break
    total_time = time.time() - start_time
    if logStats:
        log_statistics('Jacobi', total_time, residual, iterations)
    return total_time


def gauss_seidel_method(A, b, max_convergence, logStats=True):
    x = [0.0] * len(A)
    prev_x = [1.0] * len(A)
    iterations = 0
    start_time = time.time()
    while True:
        iterations += 1
        for row_idx, rows in enumerate(A):
            row_sum = 0
            for col_idx in range(row_idx):
                row_sum += A[row_idx][col_idx] * x[col_idx]
            for col_idx in range(row_idx + 1, len(A)):
                row_sum += A[row_idx][col_idx] * prev_x[col_idx]
            x[row_idx] = (b[row_idx] - row_sum) / A[row_idx][row_idx]
        for idx, _ in enumerate(x):
            prev_x[idx] = x[idx]
        normalised_residual = normalise(get_residual(A, b, x))
        if normalised_residual <= max_convergence:
            residual = normalised_residual
            break
    total_time = time.time() - start_time
    if logStats:
        log_statistics('Gauss-Seidel', total_time, residual, iterations)
    return total_time


def lu_factorisation(A, b, logStats=True):
    U = deepcopy(A)
    L = generate_identity_matrix(len(A))
    matrix_size = len(A)
    start_time = time.time()

    for col_idx in range(matrix_size - 1):
        for row_idx in range(col_idx + 1, matrix_size):
            L[row_idx][col_idx] = U[row_idx][col_idx] / U[col_idx][col_idx]
            for k in range(col_idx, matrix_size):
                U[row_idx][k] = U[row_idx][k] - L[row_idx][col_idx] * U[col_idx][k]

    y = [0.0] * matrix_size
    for row_idx in range(matrix_size):
        row_sum = 0
        for col_idx in range(matrix_size):
            row_sum += L[row_idx][col_idx] * y[col_idx]
        y[row_idx] = (b[row_idx] - row_sum) / L[row_idx][row_idx]

    x = [0.0] * matrix_size
    for row_idx in range(matrix_size - 1, -1, -1):
        row_sum = 0
        for col_idx in range(row_idx + 1, matrix_size):
            row_sum += U[row_idx][col_idx] * x[col_idx]
        x[row_idx] = (y[row_idx] - row_sum) / U[row_idx][row_idx]
    end_time = time.time()
    total_time = end_time - start_time
    normalised_residual = normalise(get_residual(A, b, x))
    if logStats:
        log_statistics("LU factorisation", total_time, normalised_residual, None)
    return total_time


if __name__ == '__main__':
    print("Task A and B:")
    a_matrix = generate_a_matrix(size=constants.N, values=[5 + constants.E, -1, -1])
    b_vector = generate_b_vector(constants.N, constants.F)
    jacobi_method(a_matrix, b_vector, constants.ALT_CONVERGENCE)
    gauss_seidel_method(a_matrix, b_vector, constants.ALT_CONVERGENCE)

    print("Task C:")
    print("Because of increasing divergence during computations method will not yield solution")
    # a_matrix = generate_a_matrix(size=constants.N, values=[3, -1, -1])
    # print("JACOBI:")
    # jacobi_method(a_matrix, b_vector, constants.ALT_CONVERGENCE)
    # print("GAUSS-SEIDEL:")
    # gauss_seidel_method(a_matrix, b_vector, constants.ALT_CONVERGENCE)

    print("Task D:")
    lu_factorisation(a_matrix, b_vector)

    print("Task E:")
    sample_sizes = [100, 250, 500, 750, 1000, 1250, 1500, 1750, 2000]
    jacobi_times = []
    gauss_seidel_times = []
    lu_factorisation_times = []
    for size in sample_sizes:
        print(f'Current size: {size}')
        a_matrix = generate_a_matrix(size, values=[5 + constants.E, -1, -1])
        b_vector = generate_b_vector(size, constants.F)
        jacobi_times.append(jacobi_method(a_matrix, b_vector, constants.ALT_CONVERGENCE))
        gauss_seidel_times.append(gauss_seidel_method(a_matrix, b_vector, constants.ALT_CONVERGENCE))
        lu_factorisation_times.append(lu_factorisation(a_matrix, b_vector, logStats=False))

    plt.plot(sample_sizes, jacobi_times, label="Jacobi method")
    plt.plot(sample_sizes, gauss_seidel_times, label="Gauss-Seidel method")
    plt.plot(sample_sizes, lu_factorisation_times, label="LU factorisation method")
    plt.xlabel("Size of matrix")
    plt.ylabel("Time taken [sec]")
    plt.legend()
    plt.title("Solving systems of equations, by given methods")
    plt.savefig('systemsOfEquations.png')
    plt.show()
