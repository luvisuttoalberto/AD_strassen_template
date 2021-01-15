#ifndef __STRASSEN__

void strassen_matrix_multiplication(float **C, float const *const *const A,
                                    float const *const *const B, size_t A_rows, size_t A_columns, size_t B_columns);

void strassen_matrix_multiplication_improved_nonsquared(float **C, float const *const *const A,
                                    float const *const *const B, size_t A_rows, size_t A_columns, size_t B_columns);

void strassen_matrix_multiplication_improved(float **C, float const *const *const A,
                                    float const *const *const B, size_t A_rows, size_t A_columns, size_t B_columns);


#endif //__STRASSEN__
