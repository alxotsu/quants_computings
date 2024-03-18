from gates import *

SIGMA_X_NAME = 'sigX'
SIGMA_Y_NAME = 'sigY'
SIGMA_Z_NAME = 'sigZ'
HADAMARD_NAME = 'H'
WALSH_HADAMARD_NAME = 'WH'
CNOT_NAME = 'CN'
CNOT_NAME1 = 'CH'
CCNOT_NAME = 'CCN'
PHASE_NAME = 'PH'
PHASE_NAME1 = 'PHM'
SWAP_NAME = 'SW'
SWAP1_NAME = 'ER'
TRANSITION_NAME = '_'
T_GATE_NAME = 'T'
S_GATE_NAME = 'S'
ROOT_NAME = '√'


def add_apply_sigma_x(scheme_string, qubit_ind):
    return scheme_string + TRANSITION_NAME + SIGMA_X_NAME + '(' + str(qubit_ind + 1) + ')'


def add_apply_sigma_y(scheme_string, qubit_ind):
    return scheme_string + TRANSITION_NAME + SIGMA_Y_NAME + '(' + str(qubit_ind + 1) + ')'


def add_apply_sigma_z(scheme_string, qubit_ind):
    return scheme_string + TRANSITION_NAME + SIGMA_Z_NAME + '(' + str(qubit_ind + 1) + ')'


def add_apply_t_gate(scheme_string, qubit_ind, conjugation):
    if conjugation:
        return scheme_string + TRANSITION_NAME + T_GATE_NAME + '(' + \
               str(qubit_ind + 1) + ',' + '†' + ')'
    return scheme_string + TRANSITION_NAME + T_GATE_NAME + '(' + str(qubit_ind + 1) + ')'


def add_apply_s_gate(scheme_string, qubit_ind, conjugation):
    if conjugation:
        return scheme_string + TRANSITION_NAME + S_GATE_NAME + '(' + \
               str(qubit_ind + 1) + ',' + '†' + ')'
    return scheme_string + TRANSITION_NAME + S_GATE_NAME + '(' + str(qubit_ind + 1) + ')'


def add_apply_root(scheme_string, qubit_ind):
    return scheme_string + TRANSITION_NAME + ROOT_NAME + '(' + str(qubit_ind + 1) + ')'


def add_apply_hadamard(scheme_string, qubit_ind):
    return scheme_string + TRANSITION_NAME + HADAMARD_NAME + '(' + str(qubit_ind + 1) + ')'


def add_apply_walsh_hadamard(scheme_string):
    return scheme_string + TRANSITION_NAME + WALSH_HADAMARD_NAME + '()'


def add_apply_cnot(scheme_string, control_qubit_ind, target_qubit_ind):
    return scheme_string + TRANSITION_NAME + CNOT_NAME + '(' + \
           str(control_qubit_ind + 1) + ',' + str(target_qubit_ind + 1) + ')'

def add_apply_cnot1(scheme_string, control_qubit_ind, target_qubit_ind):
    return scheme_string + TRANSITION_NAME + CNOT_NAME1 + '(' + \
           str(control_qubit_ind + 1) + ',' + str(target_qubit_ind + 1) + ')'



def add_apply_ccnot(scheme_string, first_control_qubit_ind, second_control_qubit_ind, target_qubit_ind):
    return scheme_string + TRANSITION_NAME + CCNOT_NAME + '(' + str(first_control_qubit_ind + 1) + ',' + \
           str(second_control_qubit_ind + 1) + ',' + str(target_qubit_ind + 1) + ')'


def add_apply_phase(scheme_string, control_qubit_ind, target_qubit_ind, phase_in_fractions_of_pi, sign):
    return scheme_string + TRANSITION_NAME + PHASE_NAME + '(' + str(control_qubit_ind + 1) + ',' + \
           str(target_qubit_ind + 1) + ',' + str(phase_in_fractions_of_pi) + ',' + str(sign) + ')'

def add_apply_phase1(scheme_string, control_qubit_ind, target_qubit_ind, phase_in_fractions_of_pi):
    return scheme_string + TRANSITION_NAME + PHASE_NAME1 + '(' + str(control_qubit_ind + 1) + ',' + \
           str(target_qubit_ind + 1) + ',' + str(phase_in_fractions_of_pi) + ')'


def add_apply_swap(scheme_string, first_line_ind, second_line_ind):
    return scheme_string + TRANSITION_NAME + SWAP_NAME + '(' + str(first_line_ind + 1) + ',' + \
           str(second_line_ind + 1) + ')'

def add_apply_swap1(scheme_string, first_line_ind, second_line_ind):
    #print ('первый индекс из scheme', first_line_ind )
    #print ('второй индекс из scheme', second_line_ind )
    return scheme_string + TRANSITION_NAME + SWAP1_NAME + '(' + str(first_line_ind + 1) + ',' + \
           str(second_line_ind + 1) + ')'

def add_apply_erroup(scheme_string, first_line_ind, second_line_ind):
    return scheme_string + TRANSITION_NAME + ERROR_NAME + '(' + str(first_line_ind + 1) + ',' + \
           str(second_line_ind + 1) + ')'




def apply_scheme(scheme_string, input_array):
    gates_list = scheme_string.replace(' ', '').split('_')
    result_array = input_array
    for gate_string in gates_list:
        for array_ind in range(result_array.size):
            result_array[array_ind] = sp.simplify(result_array[array_ind])

        def is_gate(gate_name):
            return gate_string.startswith(gate_name + '(')

        def gate_args(gate_name):
            return gate_string.replace(gate_name + '(', '').replace(')', '').split(',')

        if gate_string == '':
            continue

        elif is_gate(SIGMA_X_NAME):
            args_list = gate_args(SIGMA_X_NAME)
            if len(args_list) != 1:
                return result_array, 'Wrong number of args: ' + gate_string
            result_array = apply_sigma_x(result_array, int(args_list[0]) - 1)

        elif is_gate(SIGMA_Y_NAME):
            args_list = gate_args(SIGMA_Y_NAME)
            if len(args_list) != 1:
                return result_array, 'Wrong number of args: ' + gate_string
            result_array = apply_sigma_y(result_array, int(args_list[0]) - 1)

        elif is_gate(SIGMA_Z_NAME):
            args_list = gate_args(SIGMA_Z_NAME)
            if len(args_list) != 1:
                return result_array, 'Wrong number of args: ' + gate_string
            result_array = apply_sigma_z(result_array, int(args_list[0]) - 1)

        elif is_gate(T_GATE_NAME):
            args_list = gate_args(T_GATE_NAME)
            if len(args_list) != 1 or len(args_list) != 2:
                return result_array, 'Wrong number of args: ' + gate_string
            result_array = apply_t_gate(result_array, int(args_list[0]) - 1, args_list[1])

        elif is_gate(S_GATE_NAME):
            args_list = gate_args(S_GATE_NAME)
            if len(args_list) != 1 or len(args_list) != 2:
                return result_array, 'Wrong number of args: ' + gate_string
            result_array = apply_s_gate(result_array, int(args_list[0]) - 1, args_list[1])

        elif is_gate(ROOT_NAME):
            args_list = gate_args(ROOT_NAME)
            if len(args_list) != 1:
                return result_array, 'Wrong number of args: ' + gate_string
            result_array = apply_root(result_array, int(args_list[0]) - 1)

        elif is_gate(HADAMARD_NAME):
            args_list = gate_args(HADAMARD_NAME)
            if len(args_list) != 1:
                return result_array, 'Wrong number of args: ' + gate_string
            result_array = apply_hadamard(result_array, int(args_list[0]) - 1)

        elif is_gate(WALSH_HADAMARD_NAME):
            args_list = gate_args(WALSH_HADAMARD_NAME)
            print(args_list)
            if len(args_list) != 1:
                return result_array, 'Wrong number of args: ' + gate_string
            result_array = apply_walsh_hadamard(result_array)

        elif is_gate(CNOT_NAME):
            args_list = gate_args(CNOT_NAME)
            if len(args_list) != 2:
                return result_array, 'Wrong number of args: ' + gate_string
            result_array = apply_cnot(result_array, int(args_list[0]) - 1, int(args_list[1]) - 1)

        elif is_gate(CCNOT_NAME):
            args_list = gate_args(CCNOT_NAME)
            if len(args_list) != 3:
                return result_array, 'Wrong number of args: ' + gate_string
            result_array = apply_ccnot(
                result_array, int(args_list[0]) - 1, int(args_list[1]) - 1, int(args_list[2]) - 1)

        elif is_gate(PHASE_NAME):
            args_list = gate_args(PHASE_NAME)
            if len(args_list) != 4:
                return result_array, 'Wrong number of args: ' + gate_string
            result_array = apply_phase_with_sign(result_array, int(args_list[0]) - 1, int(args_list[1]) - 1, float(args_list[2]), args_list[3])

        elif is_gate(SWAP_NAME):
            args_list = gate_args(SWAP_NAME)
            if len(args_list) != 2:
                return result_array, 'Wrong number of args: ' + gate_string
            result_array = apply_swap(result_array, int(args_list[0]) - 1, int(args_list[1]) - 1)

        elif is_gate(SWAP1_NAME):
            args_list = gate_args(SWAP1_NAME)
            if len(args_list) != 2:
                return result_array, 'Wrong number of args: ' + gate_string
            result_array = apply_swap1(result_array, int(args_list[0]) - 1, int(args_list[1]) - 1)

        elif is_gate(CNOT_NAME1):
            args_list = gate_args(CNOT_NAME1)
            if len(args_list) != 2:
                return result_array, 'Wrong number of args: ' + gate_string
            result_array = apply_cnot1(result_array, int(args_list[0]) - 1, int(args_list[1]) - 1)

        else:
            return result_array, 'Unknown operator: ' + gate_string

    return result_array, ''
