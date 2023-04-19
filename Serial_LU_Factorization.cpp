/**
 * @file Serial_LU_Factorization.cpp
 * @author Shawn George
 * @author Adelin Owona
 * @author Michael Lenyszn
 * @author Miles Corn
 * @brief This file performs the LU factorization of a matrix with partial
 * pivoting in serial
 * @version 0.1
 * @date 2023-04-05
 *
 * @copyright Copyright (c) 2023
 *
 */

#include <math.h>
#include <stdio.h>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <random>

#include "clockcycle.h"

using namespace std;

void LU(double** matrix, double** L, double** U, double** P, int dimension) {
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

    // get ready to pivot
    double rowholder;

    double lrowholder;

    double prowholder;

    // pivot - note, we have identified the row with the max index in col c
    // thus, we are swapping rows, so the col index is changing on the loop
    // since we are making triangular matrices, the row is fixed at c
    for (int k = c; k < dimension; k++) {
      rowholder = U[c][k];
      U[c][k] = U[max_index][k];
      U[max_index][k] = rowholder;
    }
    for (int k = 0; k < c - 1; k++) {
      lrowholder = L[c][k];
      L[c][k] = L[max_index][k];
      L[max_index][k] = lrowholder;
    }
    for (int k = 0; k < dimension; k++) {
      prowholder = P[c][k];
      P[c][k] = P[max_index][k];
      P[max_index][k] = prowholder;
    }

    // row operation
    for (int r = c + 1; r < dimension; r++) {
      L[r][c] = U[r][c] / U[c][c];
      for (int k = c; k < dimension; k++) {
        U[r][k] = U[r][k] - L[r][c] * U[c][k];
      }
    }
  }
}

/**
 * @brief Checks the answer of the LU factorization
 * 
 * @param matrix original matrix
 * @param answer lu decomposed matrices
 * @param dimension 
 * @return true if they are equal
 * @return false if they are not equal
 */
bool checkAnswer(double** matrix, double** answer, int dimension) {
  for (int r = 0; r < dimension; r++) {
    for (int c = 0; c < dimension; c++) {
      if (fabs(matrix[r][c] - answer[r][c]) > 0.0001) {
        return false;
      }
    }
  }

  return true;
}

/**
 * @brief Multiplies two matrices together
 * 
 * @param A First matric
 * @param B Second matrix
 * @param C Resulting matrix
 * @param dimension 
 */
void matrixMult(double** A, double** B, double** C, int dimension) {
  for (int i = 0; i < dimension; i++) {
    for (int j = 0; j < dimension; j++) {
      C[i][j] = 0;
      for (int k = 0; k < dimension; k++) {
        C[i][j] += A[i][k] * B[k][j];
      }
    }
  }
}

/**
 * @brief Function to generate a strictly diagonally dominant matrix
 * 
 * @param A  matrix to be generated
 * @param n 
 * @return 
 */
void generateSDD(double** A, int n) {
  // Initialize matrix with zeros

  // Fill diagonal with random values between 1 and 10
  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<double> dis(1.0, 10.0);
  for (int i = 0; i < n; i++) {
    A[i][i] = dis(gen);
  }
  // Fill off-diagonal with random values between -1 and 1
  uniform_real_distribution<double> dis_off(-1.0, 1.0);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i != j) {
        A[i][j] = dis_off(gen);
      }
    }
  }
  // Make the matrix strictly diagonally dominant
  for (int i = 0; i < n; i++) {
    double row_sum = 0.0;
    for (int j = 0; j < n; j++) {
      if (i != j) {
        row_sum += abs(A[i][j]);
      }
    }
    if (row_sum >= abs(A[i][i])) {
      // If the row sum is greater than or equal to the diagonal element, shift
      // the diagonal element
      A[i][i] += row_sum + 1;
    }
  }
}

int main() {
  // Testing
  int n = 800;
  double** A = new double*[n];
  for (int i = 0; i < n; i++) {
    A[i] = new double[n];
  }

  generateSDD(A, n);

  double** matrix;
  double** U;
  double** L;
  double** P;

  L = new double*[n];
  P = new double*[n];
  U = new double*[n];

  matrix = new double*[n];
  int num = 1;
  for (int r = 0; r < n; r++) {
    matrix[r] = new double[n];
    U[r] = new double[n];
    L[r] = new double[n];
    P[r] = new double[n];
    for (int c = 0; c < n; c++) {
      matrix[r][c] = A[r][c];
      U[r][c] = A[r][c];
      L[r][c] = (r == c) ? 1 : 0;
      P[r][c] = (r == c) ? 1 : 0;
      num++;
    }
  }

  unsigned long long start_time;
  unsigned long long end_time;
  double cycles_per_second = 512000000;
  double time_elapsed;

  start_time = clock_now();
  LU(matrix, L, U, P, n);
  end_time = clock_now();
  time_elapsed = (double)(end_time - start_time) / cycles_per_second;

  double** C = new double*[n];
  for (int i = 0; i < n; i++) {
    C[i] = new double[n];
  }

  double** ans = new double*[n];
  for (int i = 0; i < n; i++) {
    ans[i] = new double[n];
  }

  matrixMult(P, matrix, C, n);

  matrixMult(L, U, ans, n);

  bool sameMatrix = checkAnswer(C, ans, n);

  if (sameMatrix) {
    cout << "The matrices are the same" << endl;
  } else {
    cout << "The matrices are not the same" << endl;
  }

  printf("Time elapsed: %f seconds\n", time_elapsed);

  for (int r = 0; r < n; r++) {
    free(A[r]);
    free(matrix[r]);
    free(U[r]);
    free(L[r]);
    free(P[r]);
    free(C[r]);
    free(ans[r]);
  }

  free(A);
  free(matrix);
  free(U);
  free(L);
  free(P);
  free(C);
  free(ans);
}
