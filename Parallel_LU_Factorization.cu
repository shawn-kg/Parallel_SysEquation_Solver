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

allocate cuda managed memory for matrix
extern void matrix_cuda_alloc(double** matrix, int dimension) {
  cudaMallocManaged(matrix, dimension * dimension * sizeof(double));
}

// free cuda managed memory
extern void matrix_cuda_free(double** matrix, double** L, double** U,
                             double** P) {
  cudaFree(matrix);
  cudaFree(L);
  cudaFree(U);
  cudaFree(P);
}

// wrapper function to call LU_fact from mpi code
extern void Lu_fact_wrapper(double** matrix, double** L, double** U, double**
P,
                            int dimension) {
  double** L2D;
  double** U2D;
  double** P2D;
  double** matrix2D;

  // allocate cuda managed memory for 2D arrays
  cudaMallocManaged(&matrix2D, dimension * sizeof(double*));
  cudaMallocManaged(&L2D, dimension * sizeof(double*));
  cudaMallocManaged(&U2D, dimension * sizeof(double*));
  cudaMallocManaged(&P2D, dimension * sizeof(double*));

  // allocate cuda managed memory for 1D arrays
  cudaMallocManaged(L, dimension * dimension * sizeof(double));
  cudaMallocManaged(U, dimension * dimension * sizeof(double));
  cudaMallocManaged(P, dimension * dimension * sizeof(double));

  // convert 1D arrays to 2D arrays
  for (int i = 0; i < dimension; i++) {
    matrix2D[i] = *matrix + (i * dimension);
    L2D[i] = *L + (i * dimension);
    U2D[i] = *U + (i * dimension);
    P2D[i] = *P + (i * dimension);
  }

  // initialize L and P to identity matrices and U to the matrix
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

// cuda kernel to generate curandState objects for index in matrix used to
// generate random values
__global__ void rand_init(curandState* state) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  curand_init(1337, idx, 0, &state[idx]);
}

__global__ void generateRandomValues(double** matrix, int n,
                                     curandState* state) {
  int col = blockIdx.x * blockDim.x + threadIdx.x;
  int row = blockIdx.y * blockDim.y + threadIdx.y;

  if (row < n && col < n) {
    curandState localState = state[row * n + col];
    double rand = curand_uniform_double(&localState);
    state[row * n + col] = localState;

    // use rand to generate random values
    matrix[row][col] = rand * 9 + 1;
  }
}

