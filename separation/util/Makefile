CXX = gcc

all: matrix_multiply matrix_inverse

matrix_multiply:
	$(CXX) -DMATRIX_MUL_TEST -o test_matrix_multiply matrix.c io_util.c

matrix_inverse:
	$(CXX) -DMATRIX_INVERSE_TEST -o test_matrix_inverse matrix.c

clean:
	rm -f test_matrix_multiply test_matrix_inverse
