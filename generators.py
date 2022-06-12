from math import sin


def generate_offset_value_pairs(values: []):
    output_values = [0] * (len(values) * 2 - 1)
    offset_values = [0] * (len(values) * 2 - 1)
    center_idx = (len(output_values)) // 2
    output_values[center_idx] = values[0]
    values.reverse()
    starting_offset = -(len(offset_values) // 2)
    for idx, _ in enumerate(output_values):
        if idx != center_idx:
            output_values[idx] = values[idx % 2]
        offset_values[idx] = starting_offset + idx
    values.reverse()
    return zip(offset_values, output_values)


def generate_a_matrix(size, values: []):
    a_matrix = [[0 for _ in range(size)] for _ in range(size)]
    for (row_idx, _) in enumerate(a_matrix):
        offset_values = generate_offset_value_pairs(values)
        for (offset, value) in offset_values:
            if 0 <= row_idx + offset < size:
                a_matrix[row_idx][row_idx + offset] = value
    return a_matrix


def generate_identity_matrix(size):
    output = [[0 for _ in range(size)] for _ in range(size)]
    for row_idx, _ in enumerate(output):
        output[row_idx][row_idx] = 1
    return output


def generate_b_vector(N, F):
    output = []
    for index in range(N):
        output.append(sin((index + 1) * (F + 1)))
    return output
