# SparseMatrixLinearAlgrabra-
Documentation:

def transpose_of_csr(matrix_in_csr):
arg: Takes a matrix in CSR format and transpose it by switching the rows and colums 
returns: a new matrix in CST format 

def transpose(matrix):
arg: takes a matrix and transpose it by converting the columns into rows 

def multiply_two_matrices(first_matrix, second_matrix):
arg:  takes two matrices and multiply them 
returns a new matrix 

def multiply_matrix_times_vector(a_matrix, vector):
arg: takes a matrix and a vector 
return a new matrix 

def read_from_file(name_of_file):
arg: takes the name of the file
return a list of lists where each list is a single line in the file 

def write_to_file(name_of_file, a_list):
arg: takes the name of the file and a list of lists 
writes each list in a single line of the file 

def multi_two_csr(csr_1, csr_2):
arg: takes two matrices in csr format 
and multiply them 
the row of the first csr and the column of the second csr are where it stores the result 
of the multiplication of the column of the csr_1 time the row of csr_2.
The result is stored in a dictionary. The key is the row and column indices and the value is 
The result. Each time we try to store a result in a dictionary, we check if that key exists already.
If so, we add the existing value plus a new one. 
And the end, we return the multiplication in csr format 

def convert_matrix_to_csr(a_matrix):
arg: takes a matrix and returns it in csr format 

def wrapper_to_check_symmetric(assigned_matrix):
arg: takes a matrix 
returns True if the matrix is symmetric 

def arnoldi_iteration(A, b, n):
arg: takes a matrix, a vector and the dimension of the Krylov subspace 
returns a list of columns that are an orthonormal basis of the Krylov subspace. 
For more details, check the code, it has in line comments 

