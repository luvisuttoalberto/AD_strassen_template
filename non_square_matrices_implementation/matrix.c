#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"

void naive_matrix_multiplication(float **C, float const *const *const A,
                                float const *const *const B,
                                const size_t A_rows, const size_t A_columns, const size_t B_columns) 
{
  for(size_t y = 0; y < A_rows; y++){
    for(size_t x = 0; x < B_columns; x++){
      float value = 0.0;
      for (size_t z = 0; z < A_columns; z++){
        value += A[y][z]*B[z][x];
      }
      C[y][x] = value;
    }
  }
}

int same_matrix(float const *const *const A, float const *const *const B,
                const size_t rows, const size_t cols) {
  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++) {
      if (A[i][j] != B[i][j]) {
        return 0;
      }
    }
  }

  return 1;
}

float **allocate_matrix(const size_t rows, const size_t cols) {
  float **M = (float **)malloc(sizeof(float *) * rows);

  for (size_t i = 0; i < rows; i++) {
    M[i] = (float *)malloc(sizeof(float) * cols);
  }

  return M;
}

float **allocate_matrix_with_calloc(const size_t rows, const size_t cols) {
  float **M = (float **)malloc(sizeof(float *) * rows);

  for (size_t i = 0; i < rows; i++) {
    M[i] = (float *)calloc(cols, sizeof(float));
  }

  return M;
}

void deallocate_matrix(float **A, const size_t rows) {
  for (size_t i = 0; i < rows; i++) {
    free(A[i]);
  }

  free(A);
}

float **allocate_random_matrix(const size_t rows, const size_t cols) {
  
  float **A = allocate_matrix(rows, cols);
  
  srand(10);
  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++) {
      A[i][j] = (rand() - RAND_MAX / 2) % 5;
    }
  }

  return A;
}

float **pad_matrix(float const * const * const old_matrix, const size_t old_rows, const size_t old_columns, const size_t new_rows, const size_t new_columns, const size_t f_row, const size_t f_column){
  //insert check on the number of rows and columns
  //to be sure a padding is correctly requested


  float ** new_matrix = allocate_matrix_with_calloc(new_rows, new_columns);
  for(size_t i = 0; i < old_rows; i++){
    for(size_t j = 0; j < old_columns; j++){
      new_matrix[i][j] = old_matrix[i + f_row][j + f_column];
    }
  }

  return new_matrix;
}

void unpad_matrix(float ** unpadded_matrix, float const * const * const padded_matrix, const size_t unpadded_rows, const size_t unpadded_columns, const size_t f_row, const size_t f_column){
  //maybe insert check on the number of rows and columns
  //to be sure an unpadding is correctly requested
  
  size_t limit_rows = f_row + unpadded_rows;
  size_t limit_columns = f_column + unpadded_columns;
  for(size_t i = f_row; i < limit_rows; i++){
    for(size_t j = f_column; j < limit_columns; j++){
      unpadded_matrix[i][j] = padded_matrix[i][j];
    }
  }
  return;
}

