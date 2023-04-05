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





#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>

using namespace std;

// CUDA kernel to generate random values for the matrix
__global__ void generateRandomValues(double* A, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i < n && j < n) {
        if (i == j) {
            // Diagonal element
            curandState state;
            curand_init(clock64(), i*n+j, 0, &state);
            A[i*n+j] = curand_uniform(&state) * 9.0 + 1.0; // Random value between 1 and 10
        } else {
            // Off-diagonal element
            curandState state;
            curand_init(clock64(), i*n+j, 0, &state);
            A[i*n+j] = curand_uniform(&state) * 2.0 - 1.0; // Random value between -1 and 1
        }
    }
}

// CUDA kernel to make the matrix strictly diagonally dominant
__global__ void makeStrictlyDiagonallyDominant(double* A, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        double row_sum = 0.0;
        for (int j = 0; j < n; j++) {
            if (i != j) {
                row_sum += abs(A[i*n+j]);
            }
        }
        if (row_sum >= abs(A[i*n+i])) {
            // If the row sum is greater than or equal to the diagonal element, shift the diagonal element
            A[i*n+i] += row_sum + 1;
        }
    }
}

// Function to generate a strictly diagonally dominant matrix of size n using CUDA
void generateSDD_CUDA(int n, double **A) {
    // Allocate memory on device
    double* d_A;
    cudaMalloc(&d_A, n*n*sizeof(double));
    // Generate random values for the matrix on device
    dim3 threadsPerBlock(32, 32);
    dim3 numBlocks((n+threadsPerBlock.x-1)/threadsPerBlock.x, (n+threadsPerBlock.y-1)/threadsPerBlock.y);
    generateRandomValues<<<numBlocks, threadsPerBlock>>>(d_A, n);
    // Make the matrix strictly diagonally dominant on device
    makeStrictlyDiagonallyDominant<<<(n+511)/512, 512>>>(d_A, n);
    // Copy matrix from device to host
    cudaMemcpy(A, d_A, n*n*sizeof(double), cudaMemcpyDeviceToHost);
    // Free memory on device
    cudaFree(d_A);
}

// int main() {
//     // Generate a random 5x5 strictly diagonally dominant matrix using CUDA and print it
//     vector<vector<double>> A = generateSDD_CUDA(5);
//     for (int i = 0; i < 5; i++) {
//         for (int j = 0; j < 5; j++) {
//             cout << A[i][j] << " ";
//         }
//         cout << endl;
//     }
//     return 0;
// }