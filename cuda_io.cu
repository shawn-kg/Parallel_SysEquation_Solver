/**
 * @file cuda_io.cu
 * @author Shawn George
 * @author Adelin Owona
 * @author Michael Lenyszn
 * @author Miles Corn
 * @brief This file has helper functions for parallel read and write to files
 * @version 0.1
 * @date 2023-04-05
 * 
 * @copyright Copyright (c) 2023
 * 
 */

void flatten(double** matrix, double* flattened_matrix,int num_rows,int num_cols)
{
    double ** c_matrix;
    cudaMallocManaged(c_matrix, sizeof(double*)*num_rows);
    for(int i = 0; i < num_rows;++i)
    {
        cudaMallocManaged(c_matrix[i], sizeof(double)*num_cols);
        cudaMemcpy(c_matrix[i], matrix[i],sizeof(double)*num_cols, cudaMemcpyHostToDevice);
    }
        
    double c_flattened_matrix;
    cudaMallocManaged(c_flattened_matrix, sizeof(double*)*num_rows*num_cols);
    cuda_flatten<<<ceil(num_rows*num_cols/32),32>>>(c_matrix, c_flattened_matrix,int num_rows,int num_cols);
    cudaMemcpy(flattened_matrix, c_flattened_matrix, sizeof(double*)*num_rows*num_cols, cudaMemcpyDeviceToHost);
    cudaFree(c_flattened_matrix);
    for(int i = 0; i < num_rows;++i)
    {
        cudaFree(c_matrix[i]);
    }
    cudaFree(c_matrix);
} 

void unflatten(double** matrix, double* flattened_matrix,int num_rows,int num_cols)
{
    double ** c_matrix;
    cudaMallocManaged(c_matrix, sizeof(double*)*num_rows);
    for(int i = 0; i < num_rows;++i)
    {
        cudaMallocManaged(c_matrix[i], sizeof(double)*num_cols);
    }        
    double c_flattened_matrix;
    cudaMallocManaged(c_flattened_matrix, sizeof(double*)*num_rows*num_cols);
    cudaMemcpy(c_flattened_matrix, flattened_matrix,sizeof(double)*num_rows*num_cols, cudaMemcpyHostToDevice);
    cuda_unflatten<<<ceil(num_rows*num_cols/32),32>>>(c_matrix, c_flattened_matrix,int num_rows,int num_cols);    
    cudaFree(c_flattened_matrix);
    for(int i = 0; i < num_rows;++i)
    {
        cudaMemcpy(c_matrix[i], matrix[i],sizeof(double)*num_cols, cudaMemcpyDeviceToHost);
        cudaFree(c_matrix[i]);
    }
    cudaFree(c_matrix);
} 




__global__
void cuda_flatten(double** matrix, double* flattened_matrix,int num_rows,int num_cols)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for(int i = index; i < num_rows; i+=stride)
    {
        for(int j = index; j < num_cols; j+=stride)
        {
            flattened_matrix[i*num_cols+j] = matrix[i][j];
        }
    }
}

__global__
void cuda_unflatten(double** matrix, double* flattened_matrix,int num_rows,int num_cols)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for(int i = index; i < num_rows; i+=stride)
    {
        for(int j = index; j < num_cols; j+=stride)
        {
            matrix[i][j] = flattened_matrix[i*num_cols+j]
        }
    }
}

