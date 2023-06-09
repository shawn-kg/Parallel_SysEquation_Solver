/**
 * @file mpi_io.cpp
 * @author Michael Lenyszn
 * @brief This file performs parallel io operations for the input and output
 * matrices using MPI
 * @version 0.1
 * @date 2023-04-05
 *
 * @copyright Copyright (c) 2023
 *
 */

#include <mpi.h>

#include <algorithm>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using namespace std;

/*
This file is handles all the parallel io operations for the input and
output matrices. It is also the starting point from the program.

Parallel I/O using MPI is performed as follows
To start, all MPI ranks read in an equal portion of 
the data except for the last rank which reads in its 
portion of the data plus whatever portion of the data 
did not divide evenly among the ranks using MPI_read_at_all. 
Next all the data is given to MPI rank 0 using the MPI_Send and 
MPI_Recv functions. Once the LU factorization is completed, 
Rank 0 re-disputes all the data back to the other MPI ranks using 
MPI_Send and MPI_Recv. Once all ranks have their portion of the data 
they perform an MPI_write_at_all call to write the data back out to a file. 
*/




void Lu_fact_wrapper(double** matrix, double** L, double** U, double** P,
                     int dimension);
void matrix_cuda_alloc(double** matrix, int dimension);
void matrix_cuda_free(double** matrix, double** L, double** U, double** P);
void write(double* flattened_matrix, int num_rows, int num_cols, MPI_File fh,
           int num_ranks, unsigned long long doubles_per_rank, unsigned long long rank) {
  int num_doubles = num_cols * num_rows;
  if (rank == 0) {
    for (int i = 1; i < num_ranks; ++i) {
      if (i == num_ranks - 1)
        MPI_Send(flattened_matrix + (i * doubles_per_rank),
                 doubles_per_rank + (num_doubles % doubles_per_rank),
                 MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
      else
        MPI_Send(flattened_matrix + (i * doubles_per_rank), doubles_per_rank,
                 MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
    }
    MPI_Status status;
    MPI_File_write_at_all(fh, rank * doubles_per_rank * sizeof(double),
                          flattened_matrix, doubles_per_rank, MPI_DOUBLE,
                          &status);

  } else if (rank == num_ranks - 1) {
    double* buf =
        new double[doubles_per_rank + (num_doubles % doubles_per_rank)];
    MPI_Status status;
    MPI_Recv(buf, doubles_per_rank + (num_doubles % doubles_per_rank),
             MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
    MPI_File_write_at_all(fh, rank * doubles_per_rank * sizeof(double), buf,
                          doubles_per_rank + (num_doubles % doubles_per_rank),
                          MPI_DOUBLE, &status);
    delete[] buf;
  } else {
    double* buf = new double[doubles_per_rank];
    MPI_Status status;
    MPI_Recv(buf, doubles_per_rank, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
    MPI_File_write_at_all(fh, rank * doubles_per_rank * sizeof(double), buf,
                          doubles_per_rank, MPI_DOUBLE, &status);
    delete[] buf;
  }
}

int main(int argc, char** argv) {
  MPI_Init(NULL, NULL);
  MPI_File fh;
  int num_ranks;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  int num_rows = stoi(argv[2], nullptr, 10);
  int num_cols = stoi(argv[3], nullptr, 10);
  MPI_Status status;
  int num_doubles = num_cols * num_rows;
  int doubles_per_rank = num_doubles / num_ranks;

  if (rank == num_ranks - 1) {
    int num_doubles = num_rows * num_cols;
    unsigned long long double_per_rank = num_doubles / num_ranks;
    double* buf = new double[double_per_rank + (num_doubles % double_per_rank)];
    MPI_File_read_at_all(fh, (unsigned long long)rank * double_per_rank * sizeof(double), buf,
                         double_per_rank + (num_doubles % double_per_rank),
                         MPI_DOUBLE, &status);
    MPI_Send(buf, double_per_rank + (num_doubles % double_per_rank), MPI_DOUBLE,
             0, 0, MPI_COMM_WORLD);
    delete[] buf;

    MPI_File_close(&fh);
    MPI_File_open(MPI_COMM_WORLD, "L", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL,&fh);
    write(nullptr, num_rows, num_cols, fh, num_ranks, doubles_per_rank,(unsigned long long)rank);
    
    MPI_File_close(&fh);
    MPI_File_open(MPI_COMM_WORLD, "U", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL,&fh);
    write(nullptr, num_rows, num_cols, fh, num_ranks, doubles_per_rank,(unsigned long long)rank);

    MPI_File_close(&fh);
    MPI_File_open(MPI_COMM_WORLD, "P", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL,&fh);
    write(nullptr, num_rows, num_cols, fh, num_ranks, doubles_per_rank,(unsigned long long)rank);

  } else if (rank != 0) {
    int num_doubles = num_rows * num_cols;
    unsigned long long double_per_rank = num_doubles / num_ranks;
    double* buf = new double[double_per_rank];
    MPI_File_read_at_all(fh, (unsigned long long)rank * double_per_rank * sizeof(double), buf,
                         double_per_rank, MPI_DOUBLE, &status);
    MPI_Send(buf, double_per_rank, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    delete[] buf;
    
    MPI_File_close(&fh);
    MPI_File_open(MPI_COMM_WORLD, "L", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL,&fh);
    write(nullptr, num_rows, num_cols, fh, num_ranks, doubles_per_rank,(unsigned long long)rank);
    
    MPI_File_close(&fh);
    MPI_File_open(MPI_COMM_WORLD, "U", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL,&fh);
    write(nullptr, num_rows, num_cols, fh, num_ranks, doubles_per_rank,(unsigned long long)rank);

    MPI_File_close(&fh);
    MPI_File_open(MPI_COMM_WORLD, "P", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL,&fh);
    write(nullptr, num_rows, num_cols, fh, num_ranks, doubles_per_rank,(unsigned long long)rank);
    
  } else {
    double* flattened_matrix;
    matrix_cuda_alloc(&flattened_matrix, num_cols);
    unsigned long long double_per_rank = num_doubles / num_ranks;
    MPI_File_read_at_all(fh, (unsigned long long)rank * double_per_rank * sizeof(double),
                         flattened_matrix, double_per_rank, MPI_DOUBLE,
                         &status);

    for (int i = 1; i < num_ranks; ++i) {
      if (i == num_ranks - 1)
        MPI_Recv(flattened_matrix + i * double_per_rank,
                 double_per_rank + (num_doubles % double_per_rank), MPI_DOUBLE,
                 i, 0, MPI_COMM_WORLD, &status);
      else
        MPI_Recv(flattened_matrix + i * double_per_rank, double_per_rank,
                 MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
    }



    double* L;
    double* U;
    double* P;
    Lu_fact_wrapper(&flattened_matrix, &L, &U, &P, num_cols);
    
    MPI_File_close(&fh);
    MPI_File_open(MPI_COMM_WORLD, "L", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL,&fh);
    write(L, num_rows, num_cols, fh, num_ranks, doubles_per_rank,(unsigned long long)rank);
    
    MPI_File_close(&fh);
    MPI_File_open(MPI_COMM_WORLD, "U", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL,&fh);
    write(U, num_rows, num_cols, fh, num_ranks, doubles_per_rank,(unsigned long long)rank);

    MPI_File_close(&fh);
    MPI_File_open(MPI_COMM_WORLD, "P", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL,&fh);
    write(P, num_rows, num_cols, fh, num_ranks, doubles_per_rank,(unsigned long long)rank);
    
    
    matrix_cuda_free(&flattened_matrix, &L, &U, &P);
  }
  MPI_Finalize();

  return 0;
}
