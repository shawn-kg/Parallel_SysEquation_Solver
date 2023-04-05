
__global__
void flatten(double** matrix, double* flattened_matrix,int num_rows,int num_cols)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for(int i = index; i < num_rows; i+=stride)
    {
        for(int j = index; j < num_cols; j+=stride)
        {
            flattened_matrix[i*j+j] = matrix[i][j];
        }
    }
}

__global__
void unflatten(double** matrix, double* flattened_matrix,int num_rows,int num_cols)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for(int i = index; i < num_rows; i+=stride)
    {
        for(int j = index; j < num_cols; j+=stride)
        {
            matrix[i][j] = flattened_matrix[i*j+j]
        }
    }
}

__global__
void generateSDD()