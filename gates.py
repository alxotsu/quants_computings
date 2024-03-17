import math
import numpy as np
import sympy as sp
from bin_tools import dec_to_bin_list, bin_list_to_dec, bin_register, int_to_bits

sigma_x_matrix = np.array([[0, 1], [1, 0]])
sigma_y_matrix = np.array([[0, 0 - 1j], [0 + 1j, 0]])
sigma_z_matrix = np.array([[1, 0], [0, -1]])
hadamard_matrix = np.array([[1, 1], [1, -1]])


def apply_matrix(input_array, matrix, qubit_ind):
    ###print('======= matrix =======')
    length = int(np.log2(input_array.size))
    result = np.zeros(input_array.size, input_array.dtype)
    ###print (' matrix ',matrix  )
    ###print ('qubit_ind  ',qubit_ind )
    for array_ind in range(input_array.size):
        bin_array_ind = dec_to_bin_list(array_ind, length)
        ###print('bin_array_ind in apply_matrix',bin_array_ind )
        ###print('bin_array_ind[qubit_ind]  ',bin_array_ind[qubit_ind] )
        register = bin_register(bin_array_ind[qubit_ind])  # если аргумнт 0 [[1] [0]] if 1 [[0] [1]]
        ###print (' register in apply_matrix ',register )

        output = np.ravel(np.dot(matrix, register))
        ###print ('output 1 ',output , 'type of autput',type(output ) )

        if bin_array_ind[qubit_ind] == 1:
            output = output[::-1]
        ###print('output 2  ',output ) 
        # print (' ==========end matrix if================ ' )
        result[array_ind] += "+((" + str(output[0]) + ")*(" + input_array[array_ind] + "))"
        ###print('result 1 result[array_ind] ',result[array_ind])
        bin_array_ind[qubit_ind] = not bin_array_ind[qubit_ind]
        ###print('bin_array_ind[qubit_ind]  ',bin_array_ind[qubit_ind] )
        result[array_ind] += "+((" + str(output[1]) + ")*(" + input_array[bin_list_to_dec(bin_array_ind)] + "))"
        ###print('result 2 result[array_ind] ', result[array_ind]) 
        result[array_ind] = sp.sympify(result[array_ind].replace('I', 'j'))
        ###print('result 3 result[array_ind] ', result[array_ind])
        ###print ('================== цикл  =  ', array_ind)
    ###print (' выход окончательный ' , result)
    ###print (' ==========end matrix ================ ' ) 
    return result


def apply_sigma_x(input_array, qubit_ind):
    return apply_matrix(input_array, sigma_x_matrix, qubit_ind)


def apply_sigma_y(input_array, qubit_ind):
    return apply_matrix(input_array, sigma_y_matrix, qubit_ind)


def apply_sigma_z(input_array, qubit_ind):
    return apply_matrix(input_array, sigma_z_matrix, qubit_ind)


def apply_hadamard(input_array, qubit_ind):
    result = apply_matrix(input_array, hadamard_matrix, qubit_ind)
    for row in range(0, len(result)):
        result[row] = "1/sqrt(2)*(" + str(result[row]).replace('I', 'j') + ")"
    return result


def apply_walsh_hadamard(input_array):
    length = int(np.log2(input_array.size))
    result = input_array
    print(result)  # моя вставка
    for array_ind in range(length):
        result = apply_matrix(result, hadamard_matrix, array_ind)
    for row in range(0, len(result)):
        result[row] = "(1/sqrt(2))**" + str(length) + " *(" + str(result[row]).replace('I', 'j') + ")"
    return result


###========================= CNOT

def apply_cnot(input_array, control_qubit_ind, target_qubit_ind):
    length = int(np.log2(input_array.size))  # input_array.size - размер входного массива
    # print(' число кубит = ', length)             # lenght  это стпень 2 (8 --3,  4--2) число кубит
    result = np.zeros(input_array.size, input_array.dtype)  # кажется нули все
    # print('первый result нулевой массив типа  dtype', result )
    # print('дан массив input_array = ',input_array)   #  моя вставка
    # print(' размер входного массива = ',input_array.size) #  моя вставка
    for array_ind in range(input_array.size):
        # print('========== цикл номер ======', array_ind )
        # print('array_ind = ',array_ind)               #  номер строки данных
        bin_array_ind = dec_to_bin_list(array_ind, length)  # перевод 10 - 2 номер цикла, число кубит
        # print('первый bin_array_ind',bin_array_ind)                     #  базисные состояния
        # print('тип bin_array_ind ', type(bin_array_ind ))
        # print('bin_array_ind.shape', bin_array_ind.shape )
        # print('bin_array_ind.dtype', bin_array_ind.dtype )
        if bin_array_ind[control_qubit_ind] == 1:
            bin_array_ind[target_qubit_ind] = not bin_array_ind[target_qubit_ind]
            # print('bin_array_ind[target_qubit_ind] = ',bin_array_ind[target_qubit_ind] )   # мишень
        bin_array_ind = bin_list_to_dec(bin_array_ind)
        # print( 'после учета условий bin_array_ind '   ,bin_array_ind )
        # print( 'после учета условий bin_array_ind type =  ',type(bin_array_ind )  )
        # print('второй bin_array_ind =',   bin_array_ind  )      # упорядочение?
        result[array_ind] = input_array[bin_array_ind]
        # print('result в цикле', result)    # результат array_ind в цикле
    return result


###-------------------------------CONTROLED HADAMAR
def apply_cnot1(input_array, control_qubit_ind, target_qubit_ind):
    matrix = np.array(hadamard_matrix, dtype="<U20")
    for i in range(hadamard_matrix.shape[0]):
        for j in range(hadamard_matrix.shape[1]):
            matrix[i, j] = f"({hadamard_matrix[i, j]})/sqrt(2)"
    result = apply_controlled_matrix(input_array, matrix, target_qubit_ind, [control_qubit_ind])
    return result


###========================= CCNOT

def apply_ccnot(input_array, first_control_qubit_ind, second_control_qubit_ind, target_qubit_ind):
    length = int(np.log2(input_array.size))
    result = np.zeros(input_array.size, input_array.dtype)
    for array_ind in range(input_array.size):
        bin_array_ind = dec_to_bin_list(array_ind, length)
        if bin_array_ind[first_control_qubit_ind] == 1 and bin_array_ind[second_control_qubit_ind] == 1:
            bin_array_ind[target_qubit_ind] = not bin_array_ind[target_qubit_ind]
        bin_array_ind = bin_list_to_dec(bin_array_ind)
        result[array_ind] = input_array[bin_array_ind]
    return result


###================= controled phase with sign + in exp
def apply_phase_with_sign(input_array, control_qubit_ind, target_qubit_ind, phase_k, sign):
    # sign = 1 if sign == "+" else -1

    # if sign == "+":
    #     sign = 1
    # else:
    #     sign = -1

    sign_values = {'+': 1, '-': -1}  # Define predefined values for each sign

    # Get the predefined value based on the selected sign
    sign_value = sign_values.get(sign)

    real_phase = sign_value * 2 * math.pi / (2 ** phase_k)
    length = int(np.log2(input_array.size))
    result = input_array
    for array_ind in range(input_array.size):
        if dec_to_bin_list(array_ind, length)[control_qubit_ind] == 1 and \
                dec_to_bin_list(array_ind, length)[target_qubit_ind] == 1:
            result[array_ind] = sp.sympify("(" + str(np.exp(real_phase * 1j)) + ")*(" + input_array[array_ind] + ")")
    return result


###====================== swap normal

def apply_swap(input_array, first_line_ind, second_line_ind):
    result = np.zeros(input_array.size, input_array.dtype)
    la = input_array.size
    nq = int(np.log2(la))  # numer of qubits
    # print( 'число кубит', nq )
    if first_line_ind > nq:
        print('first_line_ind is in correct')
    elif second_line_ind > nq:
        print('second_line_ind is in correct ')
    i1 = int(first_line_ind)
    # print(' i1 = ', i1)
    i2 = int(second_line_ind)
    # print(' i2 = ', i2)
    real = AllBazisStates(nq)
    # print('вход',real)
    time = AllBazisStates(nq)
    # print('time[i1] = ',time[i1],' time[i1][0] = ',time[i1][0],' time[i1][1] = ',time[i1][1]    )
    # print('time[12] = ',time[i2],' time[i2][0] = ',time[i2][0],' time[i2][1] = ',time[i2][1]   )
    for i in range(la):
        real[i][i1] = time[i][i2]
        real[i][i2] = time[i][i1]
        bin_array_ind = bin_list_to_dec(real[i])
        result[i] = input_array[bin_array_ind]
        # print('result в цикле', result)    # результат array_ind в цикле
    # print('выход',real)
    return result


def apply_swap1(input_array, first_line_ind, second_line_ind):  # make mistake of first qubit
    result = np.zeros(input_array.size, input_array.dtype)
    la = input_array.size  # size of input array
    nq = int(np.log2(la))  # numer of qubits
    # print( 'число кубит', nq )
    # print('из gates  ', ' first line',first_line_ind )
    # print('из gates  ', ' second_line',second_line_ind)
    if first_line_ind > nq:
        print('first_line_ind is in correct')
    elif second_line_ind > nq:
        print('second_line_ind is in correct ')
    i1 = int(first_line_ind)
    # print(' i1 = first_line_ind =', i1)
    i2 = int(second_line_ind)
    # print(' i2 = ', i2)
    real = AllBazisStates(nq)
    # print('вход',real)
    time = AllBazisStates(nq)
    for i in range(la):
        # print(' первый индекс сначала', time[i][first_line_ind])
        if time[i][first_line_ind] == int(0):
            time[i][first_line_ind] = int(1)
        else:
            time[i][first_line_ind] = int(0)
        # print(' первый индекс', time[i][first_line_ind])
        # print(' тип  ', type ( time[i][first_line_ind])  )
        real[i][i1] = time[i][first_line_ind]
        # print('i= ', i, 'real[i][first_line_ind] ',real[i][i1])
        # real[i][i2]=time[i][i1]
        bin_array_ind = bin_list_to_dec(real[i])
        # print( 'bin_array_ind',bin_array_ind  )
        result[i] = input_array[bin_array_ind]
        # print('result в цикле', result)    # результат array_ind в цикле
    # print('выход real',real)
    # print('выход result',result)
    # print('до if')
    # if i2 == 0:
    # output_array_to_table(result, ui.inputTableWidget)
    # else:
    # output_array_to_table(result, ui.outputTableWidget)
    # print('имя кнопки', ' кнопка swap1 ')
    return result


def apply_erroup(input_array, first_line_ind, second_line_ind):  # make mistake of first qubit
    result = np.zeros(input_array.size, input_array.dtype)
    la = input_array.size
    nq = int(np.log2(la))  # numer of qubits
    # print( 'число кубит', nq )
    if first_line_ind > nq:
        print('first_line_ind is in correct')
    elif second_line_ind > nq:
        print('second_line_ind is in correct ')
    i1 = int(first_line_ind)
    # print(' i1 = ', i1)
    # i2=int(second_line_ind)
    # print(' i2 = ', i2)
    real = AllBazisStates(nq)
    # print('вход',real)
    time = AllBazisStates(nq)
    for i in range(la):
        # print(' первый индекс сначала', time[i][first_line_ind])
        if time[i][first_line_ind] == int(0):
            time[i][first_line_ind] = int(1)
        else:
            time[i][first_line_ind] = int(0)
        # print(' первый индекс', time[i][first_line_ind])
        # print(' тип  ', type ( time[i][first_line_ind])  )
        real[i][i1] = time[i][first_line_ind]
        # real[i][i2]=time[i][i1]
        bin_array_ind = bin_list_to_dec(real[i])
        result[i] = input_array[bin_array_ind]
        # print('result в цикле', result)    # результат array_ind в цикле
    # print('выход',real)
    return result


def is_numeric(input_array):
    result = True
    for expr in input_array:
        if not sp.sympify(expr).is_number:
            result = False
            break
    return result


def length(input_array):
    result = ''
    for expr in input_array:
        result += ' + (' + expr + ') * (' + expr + ')'
    return sp.simplify('(' + result + ') ^ (1 / 2)')


def normalize(input_array):
    input_length = length(input_array)
    # todo: find way to check by eps here for very small coeffs
    if input_length == 0:
        return input_array

    result = input_array
    for array_ind in range(input_array.size):
        result[array_ind] = sp.simplify('(' + result[array_ind] + ') / (' + str(input_length) + ')')
    return result


def AllBazisStates(n):  # gives massive states
    la = 2 ** n
    result = np.zeros(shape=(2 ** n, n), dtype=int)
    for i in range(la):
        x = dec_to_bin_list(i, n)
        # print(dec_to_bin_list(i, n))
        # print(x)
        result[i] = x
        ###print('digital = ',bin_list_to_dec(x) )
    ###print(result)
    return result


def change(n):  # n число состояний
    mas = []
    bb = []
    for i in range(n):
        mas.append(str(int(i)))
        bb.append(bin(i))
        print('число = ', mas[i], '  бинарное =', bb[i])
        # print( 'что это', bb[i].reverse()  )
    result = mas
    return result


def apply_matrix1(input_array, matrix, target_qubit_ind):  # target_qubit_ind
    print(' =========  apply_matrix1 ========  ')
    print('input_array = ', input_array)
    print('matrix = ', matrix)
    print('target_qubit_ind = ', target_qubit_ind)
    length = int(np.log2(input_array.size))
    print('length ', length)
    result = np.zeros(input_array.size, input_array.dtype)
    print('result =  ', result)
    for array_ind in range(input_array.size):
        bin_array_ind = dec_to_bin_list(array_ind, length)
        print(1, 'bin_array_ind= ', bin_array_ind)
        register = bin_register(bin_array_ind[qubit_ind])
        print(2)
        output = np.ravel(np.dot(matrix, register))
        print(3)

        if bin_array_ind[qubit_ind] == 1:
            output = output[::-1]

        result[array_ind] += "+((" + str(output[0]) + ")*(" + input_array[array_ind] + "))"
        bin_array_ind[qubit_ind] = not bin_array_ind[qubit_ind]
        result[array_ind] += "+((" + str(output[1]) + ")*(" + input_array[bin_list_to_dec(bin_array_ind)] + "))"
        result[array_ind] = sp.sympify(result[array_ind].replace('I', 'j'))
    return result


def apply_controlled_matrix(input_array, matrix, target_qubit, control_qubits):
    full_matrix = np.diag(np.full(input_array.size, 1)).astype("<U64", copy=False)
    submatrix_size = 2 ** target_qubit
    bits_length = int(np.log2(input_array.size))
    for i in range(input_array.size):
        bits = int_to_bits(i, bits_length)
        if all([bits[qbit] for qbit in control_qubits]) and not bits[target_qubit]:
            full_matrix[i, i] = matrix[0, 0]
            full_matrix[i, i + submatrix_size] = matrix[0, 1]
            full_matrix[i + submatrix_size, i] = matrix[1, 0]
            full_matrix[i + submatrix_size, i + submatrix_size] = matrix[1, 1]

    result = np.zeros(input_array.size, input_array.dtype)
    for res_elem_id in range(input_array.size):
        result_elem_str = ""
        for input_elem_id in range(input_array.size):
            result_elem_str += f"(({full_matrix[res_elem_id, input_elem_id]}) * ({input_array[input_elem_id]})) + "
        result_elem_str = result_elem_str[:-3]
        result[res_elem_id] = sp.sympify(result_elem_str.replace('I', 'j'))
    return result
