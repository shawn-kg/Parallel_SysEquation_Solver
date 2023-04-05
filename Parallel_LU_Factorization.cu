/**
 * @file Parallel_LU_Factorization.cu
 * @author Shawn George
 * @author Adelin Owona
 * @author Michael Lenyszn
 * @author Miles Corn
 * @brief This file performs LU factorization on a matrix using partial pivoting in parallel
 * @version 0.1
 * @date 2023-04-05
 * 
 * @copyright Copyright (c) 2023
 * 
 */


__global__ void check_matrix_equivalence(double **A, double **B, int dimension, bool *equal) {
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

__global__ void matrix_mult(double **A, double **B, double **C, int dimension) {
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

__global__ void swap_rows_L(int row, int max_index, int col, double** L) {
  double rowholder;
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  for (int k = index; k < col; k += stride) {
    rowholder = L[row][k];
    L[row][k] = L[max_index][k];
    L[max_index][k] = rowholder;
  }
}

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

__global__ void row_ops_kernel(int col, int dimension, double** L, double** U) {
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;

  for (int r = col + 1 + index; r < dimension; r += stride) {
    L[r][col] = U[r][col] / U[col][col];

    __syncthreads(); // make sure all threads have computed L[r][col]

    // distribute this loop across the y-dimension of the threadblock
    for (int k = col + threadIdx.y; k < dimension; k += blockDim.y) {
      U[r][k] = U[r][k] - L[r][col] * U[col][k];
    }
  }
}

void LU(double** matrix, double** L, double** U, double** P, int dimension) {
  // make sure that P,L = I and U = matrix
  for (int r = 0; r < dimension; r++) {
    for (int c = 0; c < dimension; c++) {
      if (r == c) {
        L[r][c] = 1;
        P[r][c] = 1;
      } 
			else {
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

    swap_rows_L<<<2048, 2048>>>(c, max_index, c, U);
    // cudaDeviceSynchronize();
    swap_rows_U<<<2048, 2048>>>(c, max_index, c, dimension, L);
    // cudaDeviceSynchronize();
    swap_rows_P<<<2048, 2048>>>(c, max_index, dimension, P);
    // cudaDeviceSynchronize();

    row_ops_kernel<<<2048, 2048>>>(c, dimension, L, U);
    cudaDeviceSynchronize();
  }
}

int main(int argc, char* argv[]) {}