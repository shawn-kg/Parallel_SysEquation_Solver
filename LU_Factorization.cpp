#include <math.h>
#include <cstdlib>
#include <stdio.h>

void LU(double** matrix, double ** L, double ** U, double ** P, int dimension)
{
	//make sure that L = I and U = matrix
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
			lrowholder		= L[c][k];
			prowholder		= P[c][k];

			U[c][k]			= U[max_index][k];
			L[c][k]			= L[max_index][k];
			P[c][k]			= P[max_index][k];

			U[max_index][k]	= U[c][k];
			L[max_index][k]	= L[c][k];
			P[max_index][k]	= P[c][k];
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