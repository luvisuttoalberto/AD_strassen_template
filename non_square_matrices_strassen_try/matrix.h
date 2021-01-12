#ifndef __MATRIX__
#include <stdlib.h>

void naive_matrix_multiplication(float **C, float const *const *const A,
                                float const *const *const B,
                                const size_t A_rows, const size_t A_columns, const size_t B_columns);

int same_matrix(float const *const *const A, float const *const *const B,
                const size_t rows, const size_t cols);

float **allocate_matrix(const size_t rows, const size_t cols);

float **allocate_matrix_with_calloc(const size_t rows, const size_t cols);

void deallocate_matrix(float **A, const size_t rows);


float **allocate_random_matrix(const size_t rows, const size_t cols);

float **pad_matrix(float const * const * const old_matrix, const size_t old_rows, const size_t old_columns, const size_t new_rows, const size_t new_columns, const size_t f_row, const size_t f_column);

void unpad_matrix(float ** unpadded_matrix, float const * const * const padded_matrix, const size_t unpadded_rows, const size_t unpadded_columns, const size_t f_row, const size_t f_column);

#endif //__MATRIX__
