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
#include "cuda_io.h"

#include <iostream>
#include <random>

using namespace std;

__global__ void check_matrix_equivalence(double** A, double** B, bool* equal,
                                         int dimension) {
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  *equal = true;

  for (int r = index; r < dimension; r += stride) {
    for (int c = 0; c < dimension; c++) {
      if (fabs(A[r][c] - B[r][c]) > 0.0001) {
        printf("A[%d][%d] = %f, B[%d][%d] = %f",r,c,A[r][c],r,c,B[r][c]);
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

    __syncthreads();  // make sure all threads have computed L[r][col]

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
      U[r][c] = matrix[r][c];
    }
  }

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

    swap_rows_U<<<1, 32>>>(c, max_index, c, dimension, U);

    swap_rows_L<<<1, 32>>>(c, max_index, c, L);

    swap_rows_P<<<1, 32>>>(c, max_index, dimension, P);

    row_ops_kernel<<<1, 32>>>(c, dimension, L, U);
    cudaDeviceSynchronize();
  }
}

// cuda kernel to generate random values for a 2d matrix
// __global__ void generateRandomValues(double** matrix, int n, int seed,
//                                      curandState_t* state) {
//   int row = blockIdx.x * blockDim.x + threadIdx.x;
//   int col = blockIdx.y * blockDim.y + threadIdx.y;

//   // Fill diagonal with random values between 1 and 10
//   random_device rd;
//   mt19937 gen(rd());
//   uniform_real_distribution<double> dis(1.0, 10.0);

//   // Fill off-diagonal with random values between -1 and 1
//   uniform_real_distribution<double> dis_off(-1.0, 1.0);
//   if (row < n && col < n) {
//     // use rand to generate random values
//     if (row == col) {
//       matrix[row][col] = dis(gen);
//     } else {
//       matrix[row][col] = dis_off(gen);
//     }
//   }
// }

// cuda kernel to generate a strictly diagonally dominant matrix
__global__ void generateSDDMatrix(double** matrix, int n) {
  int row = blockIdx.x * blockDim.x + threadIdx.x;
  int col = blockIdx.y * blockDim.y + threadIdx.y;

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

// __global__ void generateRandomValues(double** matrix, int n,
//                                      curandState* state) {
//   int row = blockIdx.x * blockDim.x + threadIdx.x;
//   int col = blockIdx.y * blockDim.y + threadIdx.y;
//   int index = row * n + col;

//   if (col < n && row < n) {
//     curandState localState = state[index];
//     matrix[row][col] = curand_uniform(&localState) * 100.0;
//     state[index] = localState;
//   }
// }

int main(int argc, char* argv[]) {
  // initialize matrix A using cudaMallocManaged
  int dimension = 3;
  double** A;
  double** L;
  double** U;
  double** P;
  double** LU;
  double** PA;
  double * randMatrix;
  curandState* state;
  bool* equal;

  cudaMallocManaged(&A, dimension * sizeof(double*));
  cudaMallocManaged(&L, dimension * sizeof(double*));
  cudaMallocManaged(&U, dimension * sizeof(double*));
  cudaMallocManaged(&P, dimension * sizeof(double*));
  cudaMallocManaged(&LU, dimension * sizeof(double*));
  cudaMallocManaged(&PA, dimension * sizeof(double*));
  cudaMallocManaged(&randMatrix, dimension * dimension * sizeof(double))
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

  // initialize A to be this matrix
  // A = [ 2  1  1 ]
  //     [ 4  3  3 ]
  //     [ 8  7  9 ]
  // A[0][0] = 2;
  // A[0][1] = 1;
  // A[0][2] = 1;
  // A[1][0] = 4;
  // A[1][1] = 3;
  // A[1][2] = 3;
  // A[2][0] = 8;
  // A[2][1] = 7;
  // A[2][2] = 9;

  curandGenerator_t prng;
  curandCreateGenerator(&prng, CURAND_RNG_PSEUDO_DEFAULT);
  curandSetPseudoRandomGeneratorSeed(prng, (unsigned long long)clock());
  curandGenerateUniformDouble(prng, randMatrix, dimension * dimension);

  cuda_unflatten<<<1,32>>>(A, randMatrix, dimension, dimension);

  // generate a strictly diagonally dominant matrix
  generateSDDMatrix<<<1, 32>>>(A, dimension);

  // print A
  printf("A = \n");
  print_matrix(A, dimension);

  // LU factorization
  // LU_fact(A, L, U, P, dimension);

  // // print results
  // printf("L = \n");
  // print_matrix(L, dimension);

  // printf("\nU = \n");
  // print_matrix(U, dimension);

  // // compute LU and PA
  // matrix_mult<<<1, 32>>>(L, U, LU, dimension);
  // // cudaDeviceSynchronize();
  // matrix_mult<<<1, 32>>>(P, A, PA, dimension);
  // cudaDeviceSynchronize();

  // // check if LU = PA using check_matrix_equivalence
  // check_matrix_equivalence<<<1, 32>>>(LU, PA, equal, dimension);
  // cudaDeviceSynchronize();

  // // print results
  // if (*equal) {
  //   printf("\nLU = PA\n");
  // } else {
  //   printf("\nLU != PA\n");
  // }

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
  cudaFree(randMatrix)
  cudaFree(equal);
  cudaFree(state);

  return 0;
}