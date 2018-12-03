import utility as util
import copy
import numpy as np
# from scipy.linalg import eigh_tridiagonal
from scipy import linalg as LA


# takes a matrix of integers as arguments
# each entry of the matrix in the data file is separated by single space

# Features:
# 1: Read txt file and convert it to a matrix
# 2: Perform transpose of a matrix
# 3: Multiply two matrices

matrix_1 = [
    [7, 8],
    [9, 10],
    [11, 12]
]

a_symmetric_matrix = [
    [1, 1, 2],
    [1, 2, 1],
    [2, 1, 1]
]

a_vector = [7, 9, 11]

vector_size_two = [1, 2]

CSR_1 = [
    [0, 1, 1],
    [0, 2, 2],
    [1, 0, 3]
]

# arg: takes a matrix in csr format
# transpose the matrix
def transpose_of_csr(matrix_in_csr):
    copy_of_matrix = copy.deepcopy(matrix_in_csr)
    for row in copy_of_matrix:
        new_col = row[0]
        row[0] = row[1]
        row[1] = new_col
    return copy_of_matrix


def transpose(matrix):
    transpose_matrix = []
    row_length = len(matrix)
    col_length = len(matrix[0])
    for i in range(col_length):
        temp_row = []
        for j in range(row_length):
            temp_row.append(matrix[j][i])
        transpose_matrix.append(temp_row)
    return transpose_matrix


def multiply_two_matrices(first_matrix, second_matrix):
    row_dimension = len(first_matrix)
    col_dimension = len(second_matrix)  # type: int
    times = len(second_matrix[0])
    # assert row_dimension == col_dimension

    result = []
    for i in range(row_dimension):
        temp_row = []
        for j in range(times):
            sum = 0
            for k in range(col_dimension):
                sum += first_matrix[i][k] * second_matrix[k][j]
            temp_row.append(sum)
        result.append(temp_row)
    return result


def multiply_matrix_times_vector(a_matrix, vector):
    row_dimension = len(a_matrix)
    col_dimension = len(vector)
    assert len(a_matrix[0]) == col_dimension
    result = []
    for i in range(row_dimension):
        total = 0
        for j in range(col_dimension):
            total += a_matrix[i][j] * vector[j]
        result.append([total])
    return result


def read_from_file(name_of_file):
    a_matrix_from_file = []
    with open(name_of_file) as f:
        for line in f:
            remove_end_line = line.rstrip()  # removes end of new line
            data = remove_end_line.split(' ')  # split line by space
            row = map(int, data)  # convert data to Int
            a_matrix_from_file.append(row)
            # print(line)
    return a_matrix_from_file


def write_to_file(name_of_file, a_list):
    f = open(name_of_file, "w")
    for elem in a_list:
        sanitized_data = str(elem)
        sanitized_data = sanitized_data.replace(',', '').replace('[', '').replace(']', '')
        f.write(sanitized_data + "\n")


def multi_two_csr(csr_1, csr_2):
    a_map = {}
    final_result = []
    for row in csr_1:
        column_val = row[1]
        for inner_row in csr_2:
            if inner_row[0] == column_val:
                temp_sum = row[2] * inner_row[2]
                row_index = row[0]
                col_index = inner_row[1]
                if (row_index, col_index) in a_map:
                    new_sum = temp_sum + a_map[(row_index , col_index)]
                    a_map[(row_index, col_index)] = new_sum
                else:
                    a_map[(row_index, col_index)] = temp_sum

    for key, value in a_map.iteritems():
        final_result.append([key[0], key[1], value])
    final_result.sort()
    write_to_file("result.txt", final_result)
    return a_map


def convert_matrix_to_csr(a_matrix):
    result_in_csr = []
    for i in range(len(a_matrix)):
        for j in range(len(a_matrix[0])):
            if a_matrix[i][j] != 0:
                temp = [i, j, a_matrix[i][j]]
                result_in_csr.append(temp)
    return result_in_csr


def check_if_symmetric(a_matrix, length):
    for i in range(length):
        for j in range(length):
            if a_matrix[i][j] != a_matrix[j][i]:
                return False
    return True


def wrapper_to_check_symmetric(assigned_matrix):
    if check_if_symmetric(assigned_matrix, len(assigned_matrix)):
        print("It is symmetric")
    else:
        print("It's not symmetric ")


matrix_from_file = read_from_file("data.txt")


# A: mxm array
# b: = initial vector (length m)
# n: dimension of the Krylov subspace
# returns a list of columns that are an orthonormal basis of the
# Krylov subspace
def arnoldi_iteration(A, b, n):

    A = np.array(A)
    m = A.shape[0]

    # h = np.zeros((n + 1, n), dtype=np.complex)
    # Q = np.zeros((m, n + 1), dtype=np.complex)

    h = np.zeros((n + 1, n))        # A on basis Q. It is upper Hessenberg
    Q = np.zeros((m, n + 1))

    q = b / np.linalg.norm(b)       # normalize the input vector
    Q[:, 0] = q                     # use it as the first Krylov vector

    for k in range(n):
        v = A.dot(q)                # generate a new candidate vector
        for j in range(k + 1):      # subtract the projection on the previous vector
            h[j, k] = np.dot(Q[:, j].conj(), v)  # <-- Q needs conjugation! q1 * v2
            v = v - h[j, k] * Q[:, j]            # q2' = v2 - (q1 * v2) * q1

        h[k + 1, k] = np.linalg.norm(v)          # q2
        eps = 1e-12                 # 10 ^ -12, if v is shorter than this threshold it is the zero vector
        if h[k + 1, k] > eps:       # Add the produced vector to the list, unless
            q = v / h[k + 1, k]     # the zero vector is produced
            Q[:, k + 1] = q
        else:                       # if that happens stop iterating
            return Q, h
    return Q, h

orthonormal, upper_hessenberg = arnoldi_iteration(a_symmetric_matrix, a_vector, 4)
print(util.print_twod_array(orthonormal))

# d = 3*np.ones(4)
# e = -1*np.ones(3)
# w, v = eigh_tridiagonal(d, e)
# A = np.diag(d) + np.diag(e, k=1) + np.diag(e, k=-1)
# print(w)
# e_vals, e_vecs = LA.eig(orthonormal)
# print(e_vals)
# np.allclose(A @ v - v @ np.diag(w), np.zeros((4, 4)))


# matlab code to compute eigenvalues of a tridiagonal matrix
# for n = 1:50
# A = diag(-2*ones(n,1))+ diag(ones(n-1,1),1) + diag(ones(n-1,1),-1);
# A(1,1) = -1;A(n,n) = -1;
# plot(eig(A),n,’o’)
# hold  on
# end

# Transpose a matrix
# print("Before: ")
# util.print_twod_array(matrix_1)
# new_matrix = transpose(matrix_1)
# print("After transpose: ")
# util.print_twod_array(new_matrix)


# Matrix Multiplication
# util.print_twod_array(matrix_from_file)
# print("*****")
# util.print_twod_array(matrix_1)
# print("Result: ")
# matrix_times_matrix = multiply_two_matrices(matrix_from_file, matrix_1)
# util.print_twod_array(matrix_times_matrix)
# wrapper_to_check_symmetric(matrix_times_matrix)
#

# Matrix times a vector
# matrix_times_vector = multiply_matrix_times_vector(matrix_from_file, a_vector)
# util.print_twod_array(matrix_times_vector)

# Transpose of a matrix times a vector
# transpose_matrix = transpose(matrix_from_file)
# new_matrix = multiply_matrix_times_vector(transpose_matrix, vector_size_two)
# util.print_twod_array(new_matrix)

# Transpose of CSR format
# new_csr = transpose_of_csr(CSR_1)
# util.print_twod_array(new_csr)

# Multiplication of two matrices in CSR format
# result in CSR format
csr_data_1 = read_from_file("csrData1.txt")
# util.print_twod_array(csr_data_1)

csr_data_2 = read_from_file("csrData2.txt")
# util.print_twod_array(csr_data_2)

# Multiply two csr and write the result to a file in the same format
result_of_two_csr = multi_two_csr(csr_data_1, csr_data_2)

# Converts a matrix into csr format
# csr_format = convert_matrix_to_csr(matrix_1)
# util.print_twod_array(csr_format)

# Produce A transpose in csr format
# a_transpose = transpose(matrix_1)
# csr_of_transpose = convert_matrix_to_csr(a_transpose)
# util.print_twod_array(csr_of_transpose)

# Check if a matrix is symmetric
