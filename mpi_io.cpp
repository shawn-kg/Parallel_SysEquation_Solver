#include <mpi.h>
#include <vector>
#include <random>
#include <iostream>
#include <algorithm>
#include <string>
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

void serial_flatten(double** matrix, double* flattened_matrix,int num_rows,int num_cols)
{
    for(int i = 0; i < num_rows; ++i)
    {
        for(int j = 0; j < num_cols; ++j)
        {
            flattened_matrix[(i*num_cols)+j] = matrix[i][j];
        }
    }
}

void serial_unflatten(double** matrix, double* flattened_matrix,int num_rows,int num_cols)
{
    for(int i = 0; i < num_rows; ++i)
    {
        for(int j = 0; j < num_cols; ++j)
        {
            matrix[i][j] = flattened_matrix[(i*num_cols)+j];
        }
    }
}

void write(double ** matrix, int num_rows, int num_cols,MPI_File fh, int num_ranks,int doubles_per_rank,int rank)
{
    if(rank == 0)
    {   
        double* flattend_matrix = new double[num_rows*num_cols];
        // flatten(matrix,flattend_matrix,num_rows,num_cols);
        // for(int i = 0; i < num_cols; ++i)
        // {
        //     for(int j = 0; j < num_cols; ++j)
        //         cout << matrix[i][j] << " ";
        //     cout << "\n";
        // }
        
        serial_flatten(matrix,flattend_matrix,num_rows,num_cols);
        for(int i = 0; i < num_cols*num_cols; ++i)
            cout << flattend_matrix[i] << endl;

        for(int i = 1; i < num_ranks;++i)
        {
            MPI_Send(flattend_matrix+doubles_per_rank*i,doubles_per_rank,MPI_DOUBLE,i,0,MPI_COMM_WORLD);
        }
        MPI_Status status; 
        MPI_File_write_at_all(fh,rank*doubles_per_rank*sizeof(double),flattend_matrix,doubles_per_rank, MPI_DOUBLE, &status);
    }
    else
    {
        double* buf = new double[doubles_per_rank];
        MPI_Status status; 
        MPI_Recv(buf,doubles_per_rank,MPI_DOUBLE,0,0,
            MPI_COMM_WORLD,&status);
        
        cout << "wrintg " << doubles_per_rank << endl;
        cout << "from";
        for(int i = 0; i < doubles_per_rank;++i)
            cout << buf[i] << " ";
        cout << endl;
        MPI_File_write_at_all(fh,rank*doubles_per_rank*sizeof(double),buf,doubles_per_rank,MPI_DOUBLE,&status);
    }
    
    
}



void make_test_matrix(double** matrix, int row, int col)
{
    for(int i = 0; i < row; ++i)
    {
        for(int j = 0; j < col;++j)
        {
            matrix[i][j] = 99;
        }
    }
}



// int main(int argc,char** argv)
// {
//     MPI_Init(NULL, NULL);
//     MPI_File fh;
//     int num_ranks;
//     int rank;
//     // cout << argc << endl;
//     MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_CREATE |MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
//     int num_rows = stoi(argv[2],nullptr,10);
//     int num_cols =  stoi(argv[3],nullptr,10);
//     MPI_Status status;
//     int size;
//     // MPI_File_get_size(fh,&size);
//     int num_doubles = num_cols*num_rows;
//     int doubles_per_rank = num_doubles/num_ranks;
//     if(rank == 0)
//     {
        
        
//         cout << "final\n";
        
//         double** matrix = new double*[num_rows];
//         for(int i = 0; i < num_rows; ++i)
//             matrix[i] = new double[num_cols]; 
//         make_test_matrix(matrix,num_rows,num_cols);

//         cout << "bruh\n";
//         write(matrix,num_rows,num_cols,fh,num_ranks,doubles_per_rank,rank);
//     }
    
//     if(rank != 0)
//         write(nullptr,num_rows,num_cols,fh,num_ranks,doubles_per_rank,rank);
//     MPI_Finalize();

//     return 0;

// }



int main(int argc, char** argv)
{
    MPI_Init(NULL, NULL);
    MPI_File fh;
    int num_ranks;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY , MPI_INFO_NULL, &fh);
    MPI_Status status;
    MPI_Offset size;
    int num_rows = stoi(argv[2],nullptr,10);
    int num_cols =  stoi(argv[3],nullptr,10);
    MPI_File_get_size(fh,&size);
    // cout << "szie was " << (unsigned int)size << endl;
    //read in matrix from file
    if(rank != 0)
    {
        int num_doubles = size/sizeof(double);
        int double_per_rank = num_doubles/num_ranks;
        // cout << "m\n";
        double* buf = new double[double_per_rank];
        // cout << "n\n";
        MPI_File_read_all(fh,buf,double_per_rank,MPI_DOUBLE,&status);//fix offset
        // cout << "v\n";
        MPI_Send(buf,double_per_rank,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
        // cout << "vv\n";
        MPI_Barrier(MPI_COMM_WORLD);
    }
    else
    {
        double** matrix = new double*[num_rows];
        for(int i = 0; i < num_rows; ++i)
            matrix[i] = new double[num_cols];
        // matrix[0][0] = 0;
        // cout << "hmm\n";
        int num_doubles = size/sizeof(double);
        double* flattened_matrix = new double[num_rows*num_cols];
        
        int double_per_rank = num_doubles/num_ranks;
        MPI_File_read_all(fh,flattened_matrix,double_per_rank,MPI_DOUBLE, &status);//fix offsets
        for(int i = 1; i < num_ranks;++i)
        {
            MPI_Recv(flattened_matrix+i*double_per_rank,double_per_rank,MPI_DOUBLE,i,0,
            MPI_COMM_WORLD,&status);
        }
        for(int i =0 ; i < num_rows*num_cols;++i)
            cout << flattened_matrix[i] << " ";
        cout << endl;
        serial_unflatten(matrix,flattened_matrix,num_rows,num_cols);
        for(int i = 0; i < num_cols; ++i)
        {
            for(int j = 0; j < num_cols; ++j)
                cout << matrix[i][j] << " ";
            cout << "\n";
        }
        MPI_Barrier(MPI_COMM_WORLD);
        
        // double** L = new double*[num_rows];
        // double** U = new double*[num_rows];
        // double** P = new double*[num_rows];
        // LU(matrix,L,U,P);
        // write(vector<vector<double>> matrix, int num_rows, int num_cols,MPI_File fh, int num_ranks,int num_doubles,int size);
        // write(vector<vector<double>> matrix, int num_rows, int num_cols,MPI_File fh, int num_ranks,int num_doubles,int size);
        // write(vector<vector<double>> matrix, int num_rows, int num_cols,MPI_File fh, int num_ranks,int num_doubles,int size);



    }
    MPI_Finalize();

    return 0;
}

    
    
    
    
    

    
//     //need to make parallel
//     vector<vector<double>> matrix;
//     for(int i = 0; i < num_rows,++i)
//     {
//         vector<double> row;
//         for(int j = 0; j < num_cols,++j)
//         {
//             matrix.push(buf[i+j]);
//         }
//         matrix.push(row);
//     }
// }