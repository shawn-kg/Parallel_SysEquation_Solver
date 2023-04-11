/**
 * @file cuda_sdd.cu
 * @author Shawn George
 * @author Adelin Owona
 * @author Michael Lenyszn
 * @author Miles Corn
 * @brief This file develops a CUDA program to
    generate a strictly diagonally dominant matrix of size n
 * @version 0.1
 * @date 2023-04-05
 *
 * @copyright Copyright (c) 2023
 *
 */

#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>

#include <algorithm>
#include <iostream>
#include <random>
#include <vector>

using namespace std;

// cuda kernel to generate random values for a 2d matrix
__global__ void generateRandomValues(int **matrix, int n, int seed) {
  int row = blockIdx.x * blockDim.x + threadIdx.x;
  int col = blockIdx.y * blockDim.y + threadIdx.y;

  if (row < n && col < n) {
    curandState_t state;
    curand_init(seed, row, col, &state);
    matrix[row][col] = curand(&state) % 100;
  }
}

// cuda kernel to generate a strictly diagonally dominant matrix
__global__ void generateSDDMatrix(int **matrix, int n) {
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