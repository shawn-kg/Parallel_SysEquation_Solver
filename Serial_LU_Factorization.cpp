#include <math.h>
#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <random>
#include <algorithm>

using namespace std;

void LU(double** matrix, double ** L, double ** U, double ** P, int dimension)
{
	//make sure that P,L = I and U = matrix
	for (int r = 0; r < dimension; r++)
	{
		for (int c = 0; c < dimension; c++)
		{
			if (r == c)
			{
				L[r][c] = 1;
				P[r][c] = 1;
			}
			else
			{
				L[r][c] = 0;
				P[r][c] = 0;
			}
			U[r][c]		= matrix[r][c];
		}
	}

	//begin factorization with partial pivoting
	for (int c = 0; c < dimension-1; c++)
	{
		double max		= fabs(U[c][c]);
		int max_index	= c;
		//find the max for the partial pivot
		for (int r = c; r < dimension - 1; r++)
		{
			if (fabs(U[r][c]) > max)
			{
				max			= fabs(U[r][c]);
				max_index	= r;
			}

		}

		//get ready to pivot
		double rowholder;

		double lrowholder;
		
		double prowholder;

		//pivot - note, we have identified the row with the max index in col c
		//thus, we are swapping rows, so the col index is changing on the loop
		//since we are making triangular matrices, the row is fixed at c
		for (int k = c; k < dimension; k++)
		{
			rowholder		= U[c][k];
			U[c][k]			= U[max_index][k];
			U[max_index][k]	= rowholder;
		}
		for (int k = 0; k < c-1; k++)
		{
			lrowholder  	= L[c][k];
			L[c][k] 		= L[max_index][k];
			L[max_index][k] = lrowholder;
		}
		for (int k = 0; k < dimension; k++)
		{
			prowholder 		= P[c][k];
			P[c][k]			= P[max_index][k];
			P[max_index][k] = prowholder;
		}
		
		//row operation
		for (int r = c + 1; r < dimension; r++)
		{
			L[r][c]			= U[r][c] / U[c][c];
			for (int k = c; k < dimension; k++)
			{
				U[r][k]		= U[r][k] - L[r][c]*U[c][k];
			}
		}
	}
}

bool checkAnswer(double ** matrix, double ** answer, int dimension)
{
	for (int r = 0; r < dimension; r++)
	{
		for (int c = 0; c < dimension; c++)
		{
			if (fabs(matrix[r][c] - answer[r][c]) > 0.0001)
			{
				return false;
			}
		}
	}

	return true;
}


void matrixMult(double ** A, double ** B, double ** C, int dimension)
{

	for (int i = 0; i < dimension; i++)
	{
		for (int j = 0; j < dimension; j++)
		{
			C[i][j] = 0;
			for (int k = 0; k < dimension; k++)
			{
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}
	
}

// Function to generate a strictly diagonally dominant matrix of size n
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
            // If the row sum is greater than or equal to the diagonal element, shift the diagonal element
            A[i][i] += row_sum + 1;
        }
    }
}

int main()
{
	int n = 3;
    double **A = new double*[n];
    for (int i = 0; i < n; i++) {
        A[i] = new double[n];
    }
    generateSDD(A, n);

	double ** matrix;
	double ** U;
	double ** L;
	double ** P;
	
	L = new double*[3];
	P = new double*[3];
	U = new double*[3];
	
	
	matrix = new double*[3];
	int num = 1;
	for (int r = 0; r < 3; r++)
	{
		matrix[r] = new double[3];
		U[r] = new double[3];
		L[r] = new double[3];
		P[r] = new double[3];
		for (int c =0;c<3;c++)
		{
			matrix[r][c] = A[r][c];
			U[r][c] = A[r][c];
			L[r][c] = (r == c) ? 1 : 0;
			P[r][c] = (r == c) ? 1 : 0;
			num++;
		}
	}

	LU(matrix, L, U, P, 3);
	
	double ** C = new double*[3];
	for (int i = 0; i < 3; i++)
	{
		C[i] = new double[3];
	}

	double ** ans = new double*[3];
	for (int i = 0; i < 3; i++)
	{
		ans[i] = new double[3];
	}
	
	matrixMult(P, matrix, C, 3);

	matrixMult(L, U, ans, 3);

	bool sameMatrix = checkAnswer(matrix, ans, 3);

	cout << "Matrix: " << endl;
	for (int i=0; i<3; i++)
	{
		for (int j=0; j<3; j++)
		{
			cout << matrix[i][j] << " ";
		}
		cout << endl;
	}

	cout << "Answer Matrix: " << endl;
	for (int i=0; i<3; i++)
	{
		for (int j=0; j<3; j++)
		{
			cout << ans[i][j] << " ";
		}
		cout << endl;
	}
	
	if (sameMatrix)
	{
		cout << "The matrices are the same" << endl;
	}
	else
	{
		cout << "The matrices are not the same" << endl;
	}



	for (int r =0;r<3;r++)
	{
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
