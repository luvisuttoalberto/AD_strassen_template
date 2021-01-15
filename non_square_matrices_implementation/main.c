#include <stdio.h>

#include "test.h"
#include "matrix.h"
#include "strassen.h"

int main(int argc, char *argv[]) {

  size_t n = 1 << 12;

  // uncomment this for loop if you want to do more than one execution of the multiplication
  // remember to close the for loop after the deallocation of the matrices
  //for(int i = 0; i < 10; i++){ 

  float **A = allocate_random_matrix(n, n);
  float **B = allocate_random_matrix(n, n);
  float **C0 = allocate_matrix(n, n);
  float **C1 = allocate_matrix(n, n);
  float **C2 = allocate_matrix(n, n);

  
  printf("n\tStrassen's Alg.\tImproved Alg.\tNon squared Alg.\tSame result\n");
  for (size_t j = 1; j <= n; j *= 2) {

    printf("%ld\t", j);
    fflush(stdout);
    
    // classic Strassen's algorithm
    printf("%lf\t", test(strassen_matrix_multiplication, C1, A, B, j, j, j));
    fflush(stdout);
    //printf("%lf\t", test(naive_matrix_multiplication, C0, A, B, j, j, j));
    //fflush(stdout);

    // improved Strassen's algorithm
    printf("%lf\t", test(strassen_matrix_multiplication_improved, C0, A, B, j, j, j));
    fflush(stdout);

    // improved Strassen's algorithm for non-squared matrices
    printf("%lf\t", test(strassen_matrix_multiplication_improved_nonsquared, C2, A, B, j, j, j));
    fflush(stdout);

    // check if the matrices obtained by the two improved algorithms are equal
    printf("%d\n", same_matrix((float const *const *const)C0, (float const *const *const)C2, j, j));
    
  }

  deallocate_matrix(A, n);
  deallocate_matrix(B, n);
  deallocate_matrix(C0, n);
  deallocate_matrix(C1, n);
  deallocate_matrix(C2, n);
  
  // uncomment the following section if you want to test the functions with non-square matrices
/*
  

  int n1 = 9; // rows of A
  int n2 = 15; // columns of A / rows of B
  int n3 = 7; // columns of B

  float ** A = allocate_random_matrix(n1, n2);
  float ** B = allocate_random_matrix(n2, n3);
  float ** C0 = allocate_matrix(n1, n3);
  float ** C1 = allocate_matrix(n1, n3);

  //convert this into a print function to avoid duplication, to be implemented in matrix.c
  printf("Matrix A:\n");
  for(size_t i = 0; i < n1; i++){
    for(size_t j = 0; j < n2; j++){
      printf("%f ", A[i][j]);
    }
    printf("\n");
  }
  printf("Matrix B:\n");
  for(size_t i = 0; i < n2; i++){
    for(size_t j = 0; j < n3; j++){
      printf("%f ", B[i][j]);
    }
    printf("\n\n");
  }
  
  test(naive_matrix_multiplication, C0, A, B, n1, n2, n3);
  test(strassen_matrix_multiplication, C1, A, B, n1, n2, n3);
  
  printf("Result of naive\n");
  for(size_t i = 0; i < n1; i++){
    for(size_t j = 0; j < n3; j++){
      printf("%f ", C0[i][j]);
    }
    printf("\n\n");
  }
  printf("Result of Strassen\n");
  for(size_t i = 0; i < n1; i++){
    for(size_t j = 0; j < n3; j++){
      printf("%f ", C1[i][j]);
    }
    printf("\n\n");
  }

  printf("Same result: %d\n", same_matrix((float const *const *const)C0, (float const *const *const)C1, n1, n3));
  
  deallocate_matrix(A, n1);
  deallocate_matrix(B, n2);
  deallocate_matrix(C0, n1);
  deallocate_matrix(C1, n1);

  */
  return 0;
}
