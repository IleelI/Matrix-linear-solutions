from math import sqrt


def get_residual(A, b, x):
    output = []
    for row_idx, rows in enumerate(A):
        row_sum = 0
        for col_idx, col in enumerate(rows):
            row_sum += col * x[col_idx]
        output.append(row_sum - b[row_idx])
    return output


def normalise(vector):
    output = 0
    for value in vector:
        output += value * value
    return sqrt(output)


def log_statistics(method_name, time, residual_value, iterations_count):
    print(f'{method_name}:')
    print(f'Elapsed time: {time} [sec]')
    print(f'Residual value: {residual_value}')
    print(f'Iterations: {iterations_count}')
    print("###")

