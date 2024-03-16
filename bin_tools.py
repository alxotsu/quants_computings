import numpy as np


def dec_to_bin_list(value, length):
    result = np.zeros(length, dtype=int)
    for number_ind in range(length):
        result[number_ind] = value % 2
        value //= 2
    result = result[::-1]
    return result


def int_to_bits(value, length):
    # not reversed list unlike another func
    result = list()
    while value > 0:
        result.append(value % 2)
        value //= 2
    while len(result) < length:
        result.append(0)
    return result


def bin_list_to_dec(value):
    result = 0
    pow = 1
    for pow_ind in value[::-1]:
        result += pow_ind * pow
        pow *= 2
    return result


def bin_register(value):
    return np.array([[1], [0]]) if value == 0 else np.array([[0], [1]])
