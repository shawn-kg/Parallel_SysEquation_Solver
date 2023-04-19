all: mpi_io.cpp Parallel_LU_Factorization.cu
	mpic++ mpi_io.cpp -c -o mpi_io.o -std=c++11
	nvcc Parallel_LU_Factorization.cu -c -o Parallel_LU_Factorization.o
	mpic++ mpi_io.o Parallel_LU_Factorization.o -o parallelSolver \
		-L/usr/local/cuda-11.2/lib64/ -lcudadevrt -lcudart -lstdc++
