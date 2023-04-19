/**
 * @file Parallel_LU_Factorization.cu
 * @author Shawn George
 * @author Adelin Owona
 * @author Michael Lenyszn
 * @author Miles Corn
 * @brief This file performs LU factorization on a matrix using partial pivoting
 * in parallel
 * @version 0.1
 * @date 2023-04-05
 *
 * @copyright Copyright (c) 2023
 *
 */

#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>

#include <iostream>
#include <random>

#include "clockcycle.h"

using namespace std;

__global__ void check_matrix_equivalence(double** A, double** B, bool* equal,
                                         int dimension) {
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  *equal = true;

  for (int r = index; r < dimension; r += stride) {
    for (int c = 0; c < dimension; c++) {
      if (fabs(A[r][c] - B[r][c]) > 0.0001) {
        printf("A[%d][%d] = %f, B[%d][%d] = %f\n", r, c, A[r][c], r, c,
               B[r][c]);
        *equal = false;
        return;
      }
    }
  }
}

// kernel to perform matrix multiplication
__global__ void matrix_mult(double** A, double** B, double** C, int dimension) {
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  for (int r = index; r < dimension; r += stride) {
    for (int c = 0; c < dimension; c++) {
      C[r][c] = 0;
      for (int k = 0; k < dimension; k++) {
        C[r][c] += A[r][k] * B[k][c];
      }
    }
  }
}

// kernel to swap rows in U
__global__ void swap_rows_U(int row, int max_index, int col, int dimension,
                            double** U) {
  double rowholder;
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;

  for (int k = col + index; k < dimension; k += stride) {
    rowholder = U[row][k];
    U[row][k] = U[max_index][k];
    U[max_index][k] = rowholder;
  }
}

// kernel to swap rows in L
__global__ void swap_rows_L(int row, int max_index, int col, double** L) {
  double rowholder;
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;

  for (int k = index; k < col - 1; k += stride) {
    rowholder = L[row][k];
    L[row][k] = L[max_index][k];
    L[max_index][k] = rowholder;
  }
}

// kernel to swap rows of P
__global__ void swap_rows_P(int row, int max_index, int dimension, double** P) {
  double rowholder;
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  for (int k = index; k < dimension; k += stride) {
    rowholder = P[row][k];
    P[row][k] = P[max_index][k];
    P[max_index][k] = rowholder;
  }
}

// kernel to perform row operations
__global__ void row_ops_kernel(int col, int dimension, double** L, double** U) {
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;

  for (int r = col + 1 + index; r < dimension; r += stride) {
    L[r][col] = U[r][col] / U[col][col];

    // distribute this loop across the y-dimension of the threadblock
    for (int k = col + threadIdx.y; k < dimension; k += blockDim.y) {
      U[r][k] = U[r][k] - L[r][col] * U[col][k];
    }
  }
}

// function to print a matrix
void print_matrix(double** matrix, int dimension) {
  for (int r = 0; r < dimension; r++) {
    for (int c = 0; c < dimension; c++) {
      printf("%f ", matrix[r][c]);
    }
    printf("\n");
  }
}



void LU_fact(double** matrix, double** L, double** U, double** P,
             int dimension) {
  // begin factorization with partial pivoting
  for (int c = 0; c < dimension - 1; c++) {
    double max = fabs(U[c][c]);
    int max_index = c;
    // find the max for the partial pivot
    for (int r = c; r < dimension - 1; r++) {
      if (fabs(U[r][c]) > max) {
        max = fabs(U[r][c]);
        max_index = r;
      }
    }

    swap_rows_U<<<2048, 2048>>>(c, max_index, c, dimension, U);

    swap_rows_L<<<2048, 2048>>>(c, max_index, c, L);

    swap_rows_P<<<2048, 2048>>>(c, max_index, dimension, P);

    row_ops_kernel<<<2048, 2048>>>(c, dimension, L, U);
    cudaDeviceSynchronize();
  }
}

extern void matrix_cuda_alloc(double ** matrix,dimension){
  cudaMallocManaged(matrix, dimension * dimension * sizeof(double));
}

extern void matrix_cuda_free(double** matrix, double **L, double **U, double **P){
  cudaFree(matrix);
  cudaFree(L);
  cudaFree(U);
  cudaFree(P);
}

extern void Lu_fact_wrapper(double** matrix, double** L, double** U, double** P,
  int dimension)
{
  double** L2D;
  double** U2D;
  double** P2D;
  double** matrix2D;

  cudaMallocManaged(&matrix2D, dimension * sizeof(double*));
  cudaMallocManaged(&L2D, dimension * sizeof(double*));
  cudaMallocManaged(&U2D, dimension * sizeof(double*));
  cudaMallocManaged(&P2D, dimension * sizeof(double*));

  cudaMallocManaged(L, dimension * dimension * sizeof(double));
  cudaMallocManaged(U, dimension * dimension * sizeof(double));
  cudaMallocManaged(P, dimension * dimension * sizeof(double));

  for (int i = 0; i < dimension; i++) {
    matrix2d[i] = matrix +(i*dimension);
    L2D[i] = L +(i*dimension)
    U2D[i] = U +(i*dimension);
    P2D[i] = P +(i*dimension);
  }

  for (int r = 0; r < dimension; r++) {
    for (int c = 0; c < dimension; c++) {
      if (r == c) {
        L2D[r][c] = 1;
        P2D[r][c] = 1;
      } else {
        L2D[r][c] = 0;
        P2D[r][c] = 0;
      }
      U2D[r][c] = matrix2D[r][c];
    }
  }

  LU_fact(matrix2D, L2D, U2D, P2D, dimension);

  cudaFree(matrix2D);
  cudaFree(L2D);
  cudaFree(U2D);
  cudaFree(P2D);

}

// cuda kernel to generate a strictly diagonally dominant matrix
__global__ void generateSDDMatrix(double** matrix, int n) {
  int col = blockIdx.x * blockDim.x + threadIdx.x;
  int row = blockIdx.y * blockDim.y + threadIdx.y;

  if (row < n && col < n) {
    if (row == col) {
      int sum = 0;
      for (int i = 0; i < n; i++) {
        if (i != row) {
          sum += matrix[row][i];
        }
      }
      matrix[row][col] = sum + 1;
    }
  }
}

__global__ void rand_init(curandState* state) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  curand_init(1337, idx, 0, &state[idx]);
}

__global__ void generateRandomValues(double** matrix, int n, curandState* state,
                                     double density) {
  int col = blockIdx.x * blockDim.x + threadIdx.x;
  int row = blockIdx.y * blockDim.y + threadIdx.y;

  if (row < n && col < n) {
    curandState localState = state[row * n + col];
    double rand = curand_uniform_double(&localState);
    state[row * n + col] = localState;

    // if (rand > density) {
    //   matrix[row][col] = 0;
    // } else {
    //   matrix[row][col] = rand * 9 + 1;
    // }

    // use rand to generate random values
    matrix[row][col] = rand * 9 + 1;
  }
}

int main(int argc, char* argv[]) {
  // initialize matrix A using cudaMallocManaged
  int dimension = 1000;
  int density = 0.1;
  double** A;
  double** L;
  double** U;
  double** P;
  double** LU;
  double** PA;
  curandState* state;
  bool* equal;

  unsigned long long start_time = clock_now();  // used to time functions
  unsigned long long end_time = clock_now();
  double cycles_per_second = 512000000;

  cudaMallocManaged(&A, dimension * sizeof(double*));
  cudaMallocManaged(&L, dimension * sizeof(double*));
  cudaMallocManaged(&U, dimension * sizeof(double*));
  cudaMallocManaged(&P, dimension * sizeof(double*));
  cudaMallocManaged(&LU, dimension * sizeof(double*));
  cudaMallocManaged(&PA, dimension * sizeof(double*));
  cudaMallocManaged(&equal, sizeof(bool));
  cudaMallocManaged(&state, dimension * dimension * sizeof(curandState));

  for (int r = 0; r < dimension; r++) {
    cudaMallocManaged(&A[r], dimension * sizeof(double));
    cudaMallocManaged(&L[r], dimension * sizeof(double));
    cudaMallocManaged(&U[r], dimension * sizeof(double));
    cudaMallocManaged(&P[r], dimension * sizeof(double));
    cudaMallocManaged(&LU[r], dimension * sizeof(double));
    cudaMallocManaged(&PA[r], dimension * sizeof(double));
  }

  // initialize curand state
  rand_init<<<2048, 2048>>>(state);

  // int block_size = 2048;

  // dim3 grid_size((dimension + block_size - 1) / block_size,
  //                (dimension + block_size - 1) / block_size);
  // dim3 blocksize(block_size, block_size);

  // initialize A to be a random matrix
  generateRandomValues<<<2048, 2048>>>(A, dimension, state, density);

  // generate a strictly diagonally dominant matrix
  generateSDDMatrix<<<2048, 2048>>>(A, dimension);
  cudaDeviceSynchronize();

  // // print A
  // printf("A = \n");
  // print_matrix(A, dimension);

  // LU factorization

  // make sure that P,L = I and U = matrix
  for (int r = 0; r < dimension; r++) {
    for (int c = 0; c < dimension; c++) {
      if (r == c) {
        L[r][c] = 1;
        P[r][c] = 1;
      } else {
        L[r][c] = 0;
        P[r][c] = 0;
      }
      U[r][c] = A[r][c];
    }
  }

  start_time = clock_now();
  LU_fact(A, L, U, P, dimension);
  end_time = clock_now();
  double time_elapsed = (double)((end_time - start_time) / cycles_per_second);

  // // print results
  // printf("L = \n");
  // print_matrix(L, dimension);

  // printf("\nU = \n");
  // print_matrix(U, dimension);

  // compute LU and PA
  matrix_mult<<<2048, 2048>>>(L, U, LU, dimension);
  // cudaDeviceSynchronize();
  matrix_mult<<<2048, 2048>>>(P, A, PA, dimension);
  cudaDeviceSynchronize();

  // check_matrix_equivalence<<<dimension, dimension>>>(U, A, equal, dimension);
  // cudaDeviceSynchronize();

  // // print results
  // if (*equal) {
  //   printf("\nU = A\n");
  // } else {
  //   printf("\nU != A\n");
  // }

  // check if LU = PA using check_matrix_equivalence
  check_matrix_equivalence<<<1024, 1024>>>(LU, PA, equal, dimension);
  cudaDeviceSynchronize();

  // print number of dimensions
  printf("Dimension: %d\n", dimension);

  // print results
  if (*equal) {
    printf("\nLU = PA\n");
  } else {
    printf("\nLU != PA\n");
  }

  printf("Time elapsed: %f seconds\n", time_elapsed);

  // printf("LU = \n");
  // print_matrix(LU, dimension);
  // printf("\nPA = \n");
  // print_matrix(PA, dimension);

  // free memory
  for (int r = 0; r < dimension; r++) {
    cudaFree(A[r]);
    cudaFree(L[r]);
    cudaFree(U[r]);
    cudaFree(P[r]);
    cudaFree(LU[r]);
    cudaFree(PA[r]);
  }

  cudaFree(A);
  cudaFree(L);
  cudaFree(U);
  cudaFree(P);
  cudaFree(LU);
  cudaFree(PA);
  cudaFree(equal);
  cudaFree(state);

  return 0;
}