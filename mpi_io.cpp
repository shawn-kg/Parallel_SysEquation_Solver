#include <mpi.h>
#include <vector>
#include <random>
#include <iostream>
#include <algorithm>
using namespace std;


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

void write(vector<vector<double>> matrix, int num_rows, int num_cols,MPI_File fh, int num_ranks,int num_doubles,int size)
{
    if(rank == 0)
    {   
        double* flattend_matrix;
        flatten(matrix,flattend_matrix,num_rows,num_cols);
        for(int i = 1; i < num_ranks;++i)
        {
            MPI_Send(flattend_matrix+num_doubles*i,num_doubles,MPI_DOUBLE,i,0,MPI_COMM_WORLD);
        }
        MPI_Status status; 
        MPI_File_write_all(fh,flattend_matrix,num_doubles, MPI_DOUBLE, &status)
    }
    else
    {
        double* buf = new double[size/sizeof(double)num_ranks];
        MPI_Recv(size,num_doubles,MPI_DOUBLE,0,0,
            MPI_COMM_WORLD,&status)
        MPI_Status status; 
        MPI_File_write_all(fh,buf,num_doubles,MPI_DOUBLE,&status)
    }
    
    
}


int main(int argc, char** argv)
{
    MPI_Init(NULL, NULL);
    MPI_Info info;
    MPI_File fh;
    int num_ranks;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_File_open(MPI_COMM_WORLD, argv[0], MPI_MODE_RDONLY , info, fh);
    MPI_Status status;
    MPI_Offset size;
    int num_rows = soti(argv[1]);
    int num_cols =  stoi(argv[2]);
    MPI_File_get_size(fh,&size);
    
    //read in matrix from file
    if(rank != 0)
    {
        double* buf = new double[size/sizeof(double)num_ranks];
        int num_doubles = size/sizeof(double);
        MPI_File_read_all(fh,buf,num_doubles,MPI_DOUBLE,MPI_Status &status);//fix offsets
        MPI_Send(buf,num_doubles,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
    }
    else
    {
        double** matrix = new double*[num_rows];
        double* flattened_matrix = new double[num_rows*num_cols];
        int num_doubles = size/sizeof(double);
        MPI_File_read_all(fh,flattened_matrix,num_doubles,MPI_DOUBLE,MPI_Status &status);//fix offsets
        for(int i = 1; i < num_ranks;++i)
        {
            MPI_Recv(flattened_matrix+i*num_doubles,num_doubles,MPI_DOUBLE,i,0,
            MPI_COMM_WORLD,&status)
        }
    
        unflatten(matrix,flattened_matrix);
        double** L = new double*[num_rows];
        double** U = new double*[num_rows];
        double** P = new double*[num_rows];
        LU(matrix);
    }

    
    
    
    
    

    
    //need to make parallel
    vector<vector<double>> matrix;
    for(int i = 0; i < num_rows,++i)
    {
        vector<double> row;
        for(int j = 0; j < num_cols,++j)
        {
            matrix.push(buf[i+j]);
        }
        matrix.push(row);
    }

    // generates a strictly diagonally dominant matrix of size n, need to make parallel
    int n = 3;
    double **A = new double*[n];
    for (int i = 0; i < n; i++) {
        A[i] = new double[n];
    }
    generateSDD(A, n);
}