/**
 * @file mpi_io.cpp
 * @author Michael Lenyszn
 * @brief This file performs parallel io operations for the input and output matrices using MPI
 * @version 0.1
 * @date 2023-04-05
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#include <mpi.h>
#include <vector>
#include <random>
#include <iostream>
#include <algorithm>
#include <string>
using namespace std;

/*
This file is handles all the parallel io operations for the input and
output matrices. It is also the starting point from the program.
*/



void LU_fact(double** matrix, double** L, double** U, double** P,
             int dimension);
void matrix_cuda_alloc(double ** matrix,dimension);
void write(double * flattened_matrix, int num_rows, int num_cols, MPI_File fh, int num_ranks, int doubles_per_rank, int rank)
{
    int num_doubles = num_cols*num_rows;
    if (rank == 0)
    {
        for (int i = 1; i < num_ranks; ++i)
        {
            if(i == num_ranks-1)
              MPI_Send(flattened_matrix+(i*doubles_per_rank), doubles_per_rank+(num_doubles%doubles_per_rank), MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            else
              MPI_Send(flattened_matrix+(i*doubles_per_rank), doubles_per_rank, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }
        MPI_Status status;
       	cout << "me with\n";
        MPI_File_write_at_all(fh, rank * doubles_per_rank * sizeof(double), flattened_matrix, doubles_per_rank, MPI_DOUBLE, &status);
        
    }
    else if(rank == num_ranks-1)
    {
        double *buf = new double[doubles_per_rank+(num_doubles%doubles_per_rank)];
        MPI_Status status;
        MPI_Recv(buf, doubles_per_rank+(num_doubles%doubles_per_rank), MPI_DOUBLE, 0, 0,
                 MPI_COMM_WORLD, &status);
	cout << "da\n";
  MPI_File_write_at_all(fh, rank * doubles_per_rank * sizeof(double), buf, doubles_per_rank+(num_doubles%doubles_per_rank), MPI_DOUBLE, &status);
      cout << "unstuck\n";
      delete [] buf;
    }
    else
    {
        double *buf = new double[doubles_per_rank];
        MPI_Status status;
        MPI_Recv(buf, doubles_per_rank, MPI_DOUBLE, 0, 0,
                 MPI_COMM_WORLD, &status);
        cout << "bois\n";
        MPI_File_write_at_all(fh, rank * doubles_per_rank * sizeof(double), buf, doubles_per_rank, MPI_DOUBLE, &status);
        cout << "unstuck\n";
      delete [] buf;
    }
}



int main(int argc,char** argv)
{
    MPI_Init(NULL, NULL);
    MPI_File fh;
    int num_ranks;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    int num_rows = stoi(argv[2],nullptr,10);
    int num_cols =  stoi(argv[3],nullptr,10);
    MPI_Status status;
    int num_doubles = num_cols*num_rows;
    int doubles_per_rank = num_doubles/num_ranks;
    cout << "ummmmm\n";
    
    if(rank == num_ranks-1)
    {
        int num_doubles = num_rows * num_cols;
        int double_per_rank = num_doubles / num_ranks;
        double * buf = new double[double_per_rank+ (num_doubles%double_per_rank)]; 
        cout << "so\n";
        MPI_File_read_at_all(fh, rank * double_per_rank*sizeof(double), buf, double_per_rank + (num_doubles%double_per_rank) , MPI_DOUBLE, &status);
        MPI_Send(buf, double_per_rank+(num_doubles%double_per_rank), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        delete buf;
        MPI_File_close(&fh);
        MPI_File_open(MPI_COMM_WORLD,"outfile",MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
        cout << "me\n";
        cout << "when"<< rank <<":"<<fh <<"\n";
        write(nullptr,num_rows,num_cols,fh,num_ranks,doubles_per_rank,rank);
        
        cout << "here" << rank << "\n";
    }
    else if(rank != 0)
    {
        int num_doubles = num_rows * num_cols;
        int double_per_rank = num_doubles / num_ranks;
        double * buf = new double[double_per_rank];
        cout << "like\n";
        MPI_File_read_at_all(fh, rank * double_per_rank*sizeof(double), buf, double_per_rank, MPI_DOUBLE, &status);
        MPI_Send(buf, double_per_rank, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        delete buf;
        MPI_File_close(&fh);
        MPI_File_open(MPI_COMM_WORLD,"outfile",MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
        
        cout << "when"<< rank <<":"<<fh <<"\n";
        write(nullptr,num_rows,num_cols,fh,num_ranks,doubles_per_rank,rank);
        cout << "here" << rank << "\n";
    }
    else
    {
        //double* flattened_matrix = new double[num_cols*num_rows];
        double* flattened_matrix;
        matrix_cuda_alloc(&flattened_matrix,num_cols);
        int double_per_rank = num_doubles / num_ranks;
        cout << "wtf\n";
        MPI_File_read_at_all(fh, rank * double_per_rank*sizeof(double), flattened_matrix, double_per_rank, MPI_DOUBLE, &status);
        MPI_File_close(&fh);
        MPI_File_open(MPI_COMM_WORLD,"outfile",MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
        cout << "when"<< rank <<":"<<fh <<"\n";
        
  for (int i = 1; i < num_ranks; ++i)
        {
            if(i == num_ranks-1)
              MPI_Recv(flattened_matrix + i * double_per_rank, double_per_rank+ (num_doubles%double_per_rank), MPI_DOUBLE, i, 0,
                     MPI_COMM_WORLD, &status);
            else
              MPI_Recv(flattened_matrix + i * double_per_rank, double_per_rank, MPI_DOUBLE, i, 0,
                     MPI_COMM_WORLD, &status);
        }

        
        // double** matrix = new double*[num_rows];
        // for(int i = 0; i < num_rows; ++i)
        //   matrix[i] =  flattened_matrix+(i*num_cols);

        double* L;
        double* U;
        double* P;
        // LU_fact(matrix,L,U,P,num_cols);
        Lu_fact_wrapper(&flattened_matrix,&L,&U,&P,num_cols);
        cout << "bruh\n";
        write(flattened_matrix,num_rows,num_cols,fh,num_ranks,doubles_per_rank,rank);
        // write(L,num_rows,num_cols,fh,num_ranks,doubles_per_rank,rank);
        // write(U,num_rows,num_cols,fh,num_ranks,doubles_per_rank,rank);
        // write(P,num_rows,num_cols,fh,num_ranks,doubles_per_rank,rank);
        // delete [] flattened_matrix;
        // delete [] matrix;
        matrix_cuda_free(&flattened_matrix,&L,&U,&P);
    }
    cout << "bye "<<rank<<"\n";
    MPI_Finalize();

    return 0;

}
