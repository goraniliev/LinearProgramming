# -*- coding: utf-8 -*-
from scipy.optimize import linprog

'''
Суровите резултати од пресметките ги печатам во конзола. Во датотеките output_1.csv и output_2.csv ги печатам резултатите
кои ни се од интерес - оптималниот (минимален) трошок и за кои вредности се добива истиот.

Притоа, во примерот 1 наоѓам цена за трансфер на стока.

Во примерот 2 податоците се однесуваат за одредување кој работник
на која задача треба да се фокусира, доколку знаеме колкава е цената за секој работник да работи на секоја од задачите и
знаеме дека бројот на работници и задачи е еднаков и притоа секој работник треба да работи на точно една задача.

Влезните датотеки имаа формат:

M N
a_1 ... a_M
b_1 ... b_N
c_1,1 ... c_1,N
...............
c_M,1 ... c_M,N
'''

__author__ = 'goran'

src_files = ['input_data1.csv', 'input_data2.csv']


def read_int_array(reader):
    return list(map(int, reader.readline().split(' ')))


def read_data(input_file):
    c = []
    cost_function = []
    with open(input_file) as reader:
        dimensions = read_int_array(reader)
        M = dimensions[0]
        N = dimensions[1]

        a = read_int_array(reader)
        b = read_int_array(reader)

        for i in range(M):
            coeff = read_int_array(reader)
            c.append(coeff)
            cost_function.extend(coeff)

    return M, N, a, b, c, cost_function

for i in range(len(src_files)):
    print('\nTest case: {}'.format(i + 1))
    input_file = src_files[i]
    M, N, a, b, c, cost_function = read_data(input_file)

    decision_variables = M * N

    matrix_coefficients = []
    upper_bounds = []
    bounds = None


    def define_lp_constraints_coefficients():
        # Ограничувања по производни капацитети:
        for i in range(M):
            coefficients = [0 for i in range(decision_variables)]
            for j in range(i * N, i * N + N):
                coefficients[j] = 1
            matrix_coefficients.append(coefficients)
            upper_bounds.append(a[i])

        # Ограничувања за продажните салони
        for j in range(N):
            coefficients = [0 for i in range(decision_variables)]
            for i in range(j, decision_variables, N):
                coefficients[i] = -1
            matrix_coefficients.append(coefficients)
            upper_bounds.append(-b[j])

        # Дефинирање на ненегативните ограничувања
        bounds = ((0, None) for i in range(decision_variables))

        print('Tableau:')
        for i in range(len(matrix_coefficients)):
            print('\t'.join([str(k) for k in matrix_coefficients[i]]), end="\t\t")
            print(upper_bounds[i])


    define_lp_constraints_coefficients()


    def solve():
        res = linprog(cost_function, A_ub=matrix_coefficients, b_ub=upper_bounds, bounds=bounds, options={"disp": True})
        print(res)
        out_file = 'output_{}.csv'.format(i + 1)
        with open(out_file, 'w') as writer:
            writer.write('X вредностите кои го даваат оптималното решение:\n')

            for j in range(0, decision_variables, N):
                writer.write(' '.join([str(k) for k in res.x[j:j + N]]) + '\n')

            writer.write('Вкупната цена е: {}'.format(res.fun))

    solve()
