all: mpi_io.cpp cuda_io.cu cuda_sdd.cu Parallel_LU_Factorization.cu
	mpi++ mpi_io.cpp -c -o mpi_io.o
	nvcc cuda_io.cu -c -o cuda_io.o
	nvcc cuda_sdd.cu -c -o cuda_sdd.o
	nvcc Parallel_LU_Factorization.cu -c -o Parallel_LU_Factorization.o
	mpi++ mpi_io.o cuda_io.o cuda_sdd.o Parallel_LU_Factorization.o -o a.out \
		-L/usr/local/cuda-11.2/lib64/ -lcudadevrt -lcudart -lstdc++