

__global__
void swap_rows_U(int row, int max_index, int col, double** U)
{
    double rowholder;
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int k = col + index; k < dimension; k+=stride)
    {
        rowholder 		= U[row][k];
		U[row][k]		= U[max_index][k];
		U[max_index][k] = rowholder;
    }


}
__global__
void swap_rows_L(int row, int max_index, int col, double** L)
{
    double rowholder;
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int k = index; k < col; k+=stride)
    {
        rowholder 		= L[row][k];
		L[row][k]		= L[max_index][k];
		L[max_index][k] = rowholder;
    }


}
__global__
void swap_rows_P(int row, int max_index, int dimension, double** P)
{
    double rowholder;
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    for (int k = index; k < dimension; k+=stride)
    {
        rowholder 		= P[row][k];
		P[row][k]		= P[max_index][k];
		P[max_index][k] = rowholder;
    }
}

__global__
void row_ops_kernel(double ** L, double ** U, int c, int dimension)
{
    int row = blockIdx.x * blockDim.x + threadIdx.x;

    if (row > c && row < dimension){
        L[row][c]			= U[row][c] / U[c][c];
        for (int k = c; k < dimension; k++)
        {
            U[row][k]		= U[row][k] - L[row][c]*U[c][k];
        }
    }	
}




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

        swap_rows_L<<<2048,2048>>>(c, max_index, c, U);
        //cudaDeviceSynchronize();
        swap_rows_U<<<2048,2048>>>(c, max_index, c, L);
        //cudaDeviceSynchronize();
        swap_rows_P<<<2048,2048>>>(c, max_index, dimension, P);
        cudaDeviceSynchronize();

		//get ready to pivot
/*		double rowholder;

		double lrowholder;
		
		double prowholder;
*/
		//pivot - note, we have identified the row with the max index in col c
		//thus, we are swapping rows, so the col index is changing on the loop
		//since we are making triangular matrices, the row is fixed at c
/*		for (int k = c; k < dimension; k++)
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
*/
		
		//row operation
		// for (int r = c + 1; r < dimension; r++)
		// {
		// 	L[r][c]			= U[r][c] / U[c][c];
		// 	for (int k = c; k < dimension; k++)
		// 	{
		// 		U[r][k]		= U[r][k] - L[r][c]*U[c][k];
		// 	}
		// }
        row_ops_kernel<<2048, 2048>>(L, U, c, dimension);
        cudaDeviceSynchronize();
	}
}



int main(int argc, char *argv[])
{


}