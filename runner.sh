#!/bin/bash -x


# This is the script that runs the experiments. It takes two arguments, the number of nodes and the number of GPUs per node. It then runs the experiments 30 times for each of the 10 different matrix sizes. The results are stored in a directory named after the number of nodes and GPUs used. The results are stored in a file named after the matrix size.

module load spectrum-mpi cuda/11.2

NUM_NODES=$1
NUM_GPUS=$2

# make sure the directory exists
mkdir -p ${NUM_NODES}Node${NUM_GPUS}GPUsExperimentsNew

# run the experiments
for i in {1..30}; do /gpfs/u/home/PCPC/PCPCmbdw/scratch/Parallel_SysEquation_Solver/parallel50 >> ${NUM_NODES}Node${NUM_GPUS}GPUsExperimentsNew/50.txt; done  
for i in {1..30}; do /gpfs/u/home/PCPC/PCPCmbdw/scratch/Parallel_SysEquation_Solver/parallel100 >> ${NUM_NODES}Node${NUM_GPUS}GPUsExperimentsNew/100.txt; done  
for i in {1..30}; do /gpfs/u/home/PCPC/PCPCmbdw/scratch/Parallel_SysEquation_Solver/parallel200 >> ${NUM_NODES}Node${NUM_GPUS}GPUsExperimentsNew/200.txt; done  
for i in {1..30}; do /gpfs/u/home/PCPC/PCPCmbdw/scratch/Parallel_SysEquation_Solver/parallel400 >> ${NUM_NODES}Node${NUM_GPUS}GPUsExperimentsNew/400.txt; done  
for i in {1..30}; do /gpfs/u/home/PCPC/PCPCmbdw/scratch/Parallel_SysEquation_Solver/parallel800 >> ${NUM_NODES}Node${NUM_GPUS}GPUsExperimentsNew/800.txt; done
for i in {1..30}; do /gpfs/u/home/PCPC/PCPCmbdw/scratch/Parallel_SysEquation_Solver/parallel1000 >> ${NUM_NODES}Node${NUM_GPUS}GPUsExperimentsNew/1000.txt; done
for i in {1..30}; do /gpfs/u/home/PCPC/PCPCmbdw/scratch/Parallel_SysEquation_Solver/parallel2000 >> ${NUM_NODES}Node${NUM_GPUS}GPUsExperimentsNew/2000.txt; done  
for i in {1..30}; do /gpfs/u/home/PCPC/PCPCmbdw/scratch/Parallel_SysEquation_Solver/parallel4000 >> ${NUM_NODES}Node${NUM_GPUS}GPUsExperimentsNew/4000.txt; done  
for i in {1..30}; do /gpfs/u/home/PCPC/PCPCmbdw/scratch/Parallel_SysEquation_Solver/parallel8000 >> ${NUM_NODES}Node${NUM_GPUS}GPUsExperimentsNew/8000.txt; done  
for i in {1..30}; do /gpfs/u/home/PCPC/PCPCmbdw/scratch/Parallel_SysEquation_Solver/parallel10000 >> ${NUM_NODES}Node${NUM_GPUS}GPUsExperimentsNew/10000.txt; done
for i in {1..30}; do /gpfs/u/home/PCPC/PCPCmbdw/scratch/Parallel_SysEquation_Solver/parallel20000 >> ${NUM_NODES}Node${NUM_GPUS}GPUsExperimentsNew/20000.txt; done