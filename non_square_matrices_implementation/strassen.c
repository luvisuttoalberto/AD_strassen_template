#include "matrix.h"

/*
 * this function performs the element-wise
 * subtraction of A and B and put the resulting
 * sub-matrix in C. The parameters *_f_row and *_f_col
 * represents the first row and the first column,
 * respectively, of the sub-matrix we want to deal with.
 */
void sub_matrix_blocks(float **C, float const * const * const A, float const * const * const B, const size_t C_f_row, const size_t C_f_col, const size_t A_f_row, const size_t A_f_col, const size_t B_f_row, const size_t B_f_col, const size_t n)
{
    // see if you can optimize these indeces too as in the function naive_aux
    for(size_t y = 0; y < n; y ++){
        for(size_t x = 0; x < n; x++){
            C[y + C_f_row][x + C_f_col] = A[y + A_f_row][x + A_f_col] - B[y + B_f_row][x + B_f_col];
        }
    }
}

/*
 * this function performs the element-wise
 * sum of A and B and put the resulting
 * sub-matrix in C. The parameters *_f_row and *_f_col
 * represents the first row and the first column,
 * respectively, of the sub-matrix we want to deal with.
 */
void sum_matrix_blocks(float **C, float const * const * const A, float const * const * const B, const size_t C_f_row, const size_t C_f_col, const size_t A_f_row, const size_t A_f_col, const size_t B_f_row, const size_t B_f_col, const size_t n)
{
    // see if you can optimize these indeces too as in the function naive_aux
    for(size_t y = 0; y < n; y ++){
        for(size_t x = 0; x < n; x++){
            C[y + C_f_row][x + C_f_col] = A[y + A_f_row][x + A_f_col] + B[y + B_f_row][x + B_f_col];
        }
    }
}

/*
 * this function implements the naive algorithm
 * for the matrix multiplication between sub-matrices.
 * The result is placed in the sub-matrix C.
 * The parameters *_f_row and *_f_col
 * represents the first row and the first column,
 * respectively, of the sub-matrix we want to deal with.
 */
void naive_aux(float **C, float const * const * const A, float const * const * const B, const size_t C_f_row, const size_t C_f_col, const size_t A_f_row, const size_t A_f_col, const size_t B_f_row, const size_t B_f_col, const size_t n)
{
    // optimize this section changing the indeces in the for loops instead of summing them in the
    // matrices indeces selections
    size_t limit_row_A = n + A_f_row;
    size_t limit_col_B = n + B_f_col;
    for(size_t y = A_f_row; y < limit_row_A; y++){
        for(size_t x = B_f_col; x < limit_col_B; x++){
            float value = 0.0;
            for (size_t z = 0; z < n; z++){
                value += A[y][z + A_f_col]*B[z + B_f_row][x];
            }
            C[y - A_f_row + C_f_row][x - B_f_col + C_f_col] = value;
        }
    }
}

/*
 * This function implements the Strassen's algorithm
 * for matrix multiplication between sub-matrices.
 * The result is placed in the sub-matrix C.
 * The parameters *_f_row and *_f_col
 * represents the first row and the first column,
 * respectively, of the sub-matrx we want to deal with. 
 */
void strassen_aux(float **C, float const * const * const A, float const * const * const B, const size_t C_f_row, const size_t C_f_col, const size_t A_f_row, const size_t A_f_col, const size_t B_f_row, const size_t B_f_col, const size_t n)
{
    if(n <= (1<<5)){
        naive_aux(C, A, B, C_f_row, C_f_col, A_f_row, A_f_col, B_f_row, B_f_col, n);
        return;
    }

    size_t n2 = n/2; // this is the size of the blocks

    float ***S = (float ***) malloc(sizeof(float **) * 10);
    for(size_t i = 0; i < 10; i++){
        S[i] = allocate_matrix(n2, n2);
    }

    float ***P = (float ***) malloc(sizeof(float **) * 7);
    for(size_t i = 0; i < 7; i++){
        P[i] = allocate_matrix(n2, n2);
    }

    // S1 = B12 - B22
    sub_matrix_blocks(S[0], B, B, 0, 0, B_f_row, B_f_col + n2, B_f_row + n2, B_f_col + n2, n2);

    // P1 = A11 * S1
    strassen_aux(P[0], A, (const float * const* const) S[0], 0, 0, A_f_row, A_f_col, 0, 0, n2);

    // S2 = A11 + A12
    sum_matrix_blocks(S[1], A, A, 0, 0, A_f_row, A_f_col, A_f_row, A_f_col + n2, n2);

    // P2 = S2 * B22
    strassen_aux(P[1], (const float * const* const) S[1], B, 0, 0, 0, 0, B_f_row + n2, B_f_col + n2, n2);

    // S3 = A21 + A22
    sum_matrix_blocks(S[2], A, A, 0, 0, A_f_row + n2, A_f_col, A_f_row + n2, A_f_col + n2, n2);

    // P3 = S3 * B11
    strassen_aux(P[2], (const float * const* const) S[2], B, 0, 0, 0, 0, B_f_row, B_f_col, n2);

    // S4 = B21 - B11
    sub_matrix_blocks(S[3], B, B, 0, 0, B_f_row + n2, B_f_col, B_f_row, B_f_col, n2);

    // P4 = A22 * S4
    strassen_aux(P[3], A, (const float * const* const) S[3], 0, 0, A_f_row + n2, A_f_col + n2, 0, 0, n2);

    // S5 = A11 + A22
    sum_matrix_blocks(S[4], A, A, 0, 0, A_f_row, A_f_col, A_f_row + n2, A_f_col + n2, n2);

    // S6 = B11 + B22
    sum_matrix_blocks(S[5], B, B, 0, 0, B_f_row, B_f_col, B_f_row + n2, B_f_col + n2, n2);

    // P5 = S5 * S6
    strassen_aux(P[4], (const float * const* const) S[4], (const float * const* const) S[5], 0, 0, 0, 0, 0, 0, n2);

    // S7 = A12 - A22
    sub_matrix_blocks(S[6], A, A, 0, 0, A_f_row, A_f_col + n2, A_f_row + n2, A_f_col + n2, n2);

    // S8 = B21 + B22
    sum_matrix_blocks(S[7], B, B, 0, 0, B_f_row + n2, B_f_col, B_f_row + n2, B_f_col + n2, n2);

    // P6 = S7 * S8
    strassen_aux(P[5], (const float * const* const) S[6], (const float * const* const) S[7], 0, 0, 0, 0, 0, 0, n2);

    // S9 = A11 - A21
    sub_matrix_blocks(S[8], A, A, 0, 0, A_f_row, A_f_col, A_f_row + n2, A_f_col, n2);

    // S10 = B11 + B12
    sum_matrix_blocks(S[9], B, B, 0, 0, B_f_row, B_f_col, B_f_row, B_f_col + n2, n2);

    // P7 = S9 * S10
    strassen_aux(P[6], (const float * const* const) S[8], (const float * const* const) S[9], 0, 0, 0, 0, 0, 0, n2);

    // C11 = P5 + P4 - P2 + P6
    sum_matrix_blocks(C, (const float * const* const) P[4], (const float * const* const) P[3], C_f_row, C_f_col, 0, 0, 0, 0, n2);
    sub_matrix_blocks(C, (const float * const* const) C, (const float * const* const) P[1], C_f_row, C_f_col, C_f_row, C_f_col, 0, 0, n2);
    sum_matrix_blocks(C, (const float * const* const) C, (const float * const* const) P[5], C_f_row, C_f_col, C_f_row, C_f_col, 0, 0, n2);

    // C12 = P1 + P2
    sum_matrix_blocks(C, (const float * const* const) P[0], (const float * const* const) P[1], C_f_row, C_f_col + n2, 0, 0, 0, 0, n2);
    
    // C21 = P3 + P4
    sum_matrix_blocks(C, (const float * const* const) P[2], (const float * const* const) P[3], C_f_row + n2, C_f_col, 0, 0, 0, 0, n2);
    
    // C22 = P5 + P1 - P3 - P7
    sum_matrix_blocks(C, (const float * const* const) P[4], (const float * const* const) P[0], C_f_row + n2, C_f_col + n2, 0, 0, 0, 0, n2);
    sub_matrix_blocks(C, (const float * const* const) C, (const float * const* const) P[2], C_f_row + n2, C_f_col + n2, C_f_col + n2, C_f_col + n2, 0, 0, n2);
    sub_matrix_blocks(C, (const float * const* const) C, (const float * const* const) P[6], C_f_row + n2, C_f_col + n2, C_f_col + n2, C_f_col + n2, 0, 0, n2);
    
    for(size_t i = 0; i < 10; i++){
        deallocate_matrix(S[i], n2);
    }
    free(S);

    for(size_t i = 0; i < 7; i++){
        deallocate_matrix(P[i], n2);
    }
    free(P);
}


/*
 * Remove the unnecessary memory allocations
 */
void strassen_aux_mem_improved(float **C, float const * const * const A, float const * const * const B, const size_t C_f_row, const size_t C_f_col, const size_t A_f_row, const size_t A_f_col, const size_t B_f_row, const size_t B_f_col, const size_t n)
{
    if(n <= (1<<5)){
        naive_aux(C, A, B, C_f_row, C_f_col, A_f_row, A_f_col, B_f_row, B_f_col, n);
        return;
    }

    size_t n2 = n/2; // this is the size of the blocks

    float ***S = (float ***) malloc(sizeof(float **) * 2);
    for(size_t i = 0; i < 2; i++){
        S[i] = allocate_matrix(n2, n2);
    }

    float ***P = (float ***) malloc(sizeof(float **) * 3);
    for(size_t i = 0; i < 3; i++){
        P[i] = allocate_matrix(n2, n2);
    }

    // S1 = B12 - B22
    sub_matrix_blocks(S[0], B, B, 0, 0, B_f_row, B_f_col + n2, B_f_row + n2, B_f_col + n2, n2);

    // P1 = A11 * S1
    strassen_aux_mem_improved(P[0], A, (const float * const* const) S[0], 0, 0, A_f_row, A_f_col, 0, 0, n2);

    // S2 = A11 + A12
    sum_matrix_blocks(S[0], A, A, 0, 0, A_f_row, A_f_col, A_f_row, A_f_col + n2, n2);

    // P2 = S2 * B22
    strassen_aux_mem_improved(P[1], (const float * const* const) S[0], B, 0, 0, 0, 0, B_f_row + n2, B_f_col + n2, n2);

    // C12 = P1 + P2
    sum_matrix_blocks(C, (const float * const* const) P[0], (const float * const* const) P[1], C_f_row, C_f_col + n2, 0, 0, 0, 0, n2);

    // S3 = A21 + A22
    sum_matrix_blocks(S[0], A, A, 0, 0, A_f_row + n2, A_f_col, A_f_row + n2, A_f_col + n2, n2);

    // P3 = S3 * B11
    strassen_aux_mem_improved(P[2], (const float * const* const) S[0], B, 0, 0, 0, 0, B_f_row, B_f_col, n2);

    // C22 = P1 - P3
    sub_matrix_blocks(C, (const float * const* const) P[0], (const float * const* const) P[2], C_f_row + n2, C_f_col + n2, 0, 0, 0, 0, n2);

    // S4 = B21 - B11
    sub_matrix_blocks(S[0], B, B, 0, 0, B_f_row + n2, B_f_col, B_f_row, B_f_col, n2);

    // P4 = A22 * S4
    strassen_aux_mem_improved(P[0], A, (const float * const* const) S[0], 0, 0, A_f_row + n2, A_f_col + n2, 0, 0, n2);

    // C11 = P4 - P2
    sub_matrix_blocks(C, (const float * const* const) P[0], (const float * const* const) P[1], C_f_row, C_f_col, 0, 0, 0, 0, n2);

    // C21 = P3 + P4
    sum_matrix_blocks(C, (const float * const* const) P[2], (const float * const* const) P[0], C_f_row + n2, C_f_col, 0, 0, 0, 0, n2);

    // P[1] and P[2] are not needed anymore
    // so we can deallocate them
    deallocate_matrix(P[2], n2);
    deallocate_matrix(P[1], n2);

    // S5 = A11 + A22
    sum_matrix_blocks(S[0], A, A, 0, 0, A_f_row, A_f_col, A_f_row + n2, A_f_col + n2, n2);

    // S6 = B11 + B22
    sum_matrix_blocks(S[1], B, B, 0, 0, B_f_row, B_f_col, B_f_row + n2, B_f_col + n2, n2);

    // P5 = S5 * S6
    strassen_aux_mem_improved(P[0], (const float * const* const) S[0], (const float * const* const) S[1], 0, 0, 0, 0, 0, 0, n2);

    // C11 += P5
    sum_matrix_blocks(C, (const float * const* const) C, (const float * const* const) P[0], C_f_row, C_f_col, C_f_col, C_f_col, 0, 0, n2);
    
    // C22 += P5
    sum_matrix_blocks(C, (const float * const* const) C, (const float * const* const) P[0], C_f_row + n2, C_f_col + n2, C_f_col + n2, C_f_col + n2, 0, 0, n2);

    // S7 = A12 - A22
    sub_matrix_blocks(S[0], A, A, 0, 0, A_f_row, A_f_col + n2, A_f_row + n2, A_f_col + n2, n2);

    // S8 = B21 + B22
    sum_matrix_blocks(S[1], B, B, 0, 0, B_f_row + n2, B_f_col, B_f_row + n2, B_f_col + n2, n2);

    // P6 = S7 * S8
    strassen_aux_mem_improved(P[0], (const float * const* const) S[0], (const float * const* const) S[1], 0, 0, 0, 0, 0, 0, n2);

    // C11 += P6
    sum_matrix_blocks(C, (const float * const* const) C, (const float * const* const) P[0], C_f_row, C_f_col, C_f_col, C_f_col, 0, 0, n2);

    // S9 = A11 - A21
    sub_matrix_blocks(S[0], A, A, 0, 0, A_f_row, A_f_col, A_f_row + n2, A_f_col, n2);

    // S10 = B11 + B12
    sum_matrix_blocks(S[1], B, B, 0, 0, B_f_row, B_f_col, B_f_row, B_f_col + n2, n2);

    // P7 = S9 * S10
    strassen_aux_mem_improved(P[0], (const float * const* const) S[0], (const float * const* const) S[1], 0, 0, 0, 0, 0, 0, n2);

    // S[i] are not needed anymore, we can deallocate them
    for(size_t i = 0; i < 2; i++){
        deallocate_matrix(S[i], n2);
    }
    free(S);

    // C22 -= P7
    sub_matrix_blocks(C, (const float * const* const) C, (const float * const* const) P[0], C_f_row + n2, C_f_col + n2, C_f_col + n2, C_f_col + n2, 0, 0, n2);

    //finally we can deallocate the last P[i] in use
    deallocate_matrix(P[0], n2);
    free(P);
}

int custom_pow_of_two(const int exponent){
    int res = 1;
    for(int j = 1; j <= exponent; j++){
        res = res*2;
    }
    return res;
}

int find_next_power_of_two_of(const int n){
    int i = 0;
    int next_power = custom_pow_of_two(i);

    while(next_power < n){
        next_power = custom_pow_of_two(++i);
    }

    return next_power;
}

/*
 * This function is exclusively meant to provide an
 * easy to use API.
 */
void strassen_matrix_multiplication(float **C, float const *const *const A, float const *const *const B, 
                                    const size_t A_rows, const size_t A_columns, const size_t B_columns) 
{

    size_t A_padded_rows = find_next_power_of_two_of(A_rows);
    size_t A_padded_columns = find_next_power_of_two_of(A_columns);
    size_t B_padded_columns = find_next_power_of_two_of(B_columns);

    size_t n = A_padded_rows >= A_padded_columns ? A_padded_rows : A_padded_columns;
    n = n >= B_padded_columns ? n : B_padded_columns;

    float ** A_padded = pad_matrix(A, A_rows, A_columns, n, n, 0, 0);
    float ** B_padded = pad_matrix(B, A_columns, B_columns, n, n, 0, 0);
    float ** C_temp = allocate_matrix(n, n);

    strassen_aux_mem_improved(C_temp, (float const * const * const) A_padded, (float const * const * const) B_padded, 0, 0, 0, 0, 0, 0, n);

    unpad_matrix(C, (float const * const * const) C_temp, A_rows, B_columns, 0, 0);

    deallocate_matrix(A_padded, n);
    deallocate_matrix(B_padded, n);
    deallocate_matrix(C_temp, n);

}

