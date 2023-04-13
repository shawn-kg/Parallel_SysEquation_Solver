#!/bin/bash -x

for i in {1..30}; do /gpfs/u/home/PCPC/PCPCrgsh/scratch/ParallelComputingProject/parallel50.o >> 8Node6GPUExp/50.txt; done  
for i in {1..30}; do /gpfs/u/home/PCPC/PCPCrgsh/scratch/ParallelComputingProject/parallel100.o >> 8Node6GPUExp/100.txt; done  
for i in {1..30}; do /gpfs/u/home/PCPC/PCPCrgsh/scratch/ParallelComputingProject/parallel200.o >> 8Node6GPUExp/200.txt; done  
for i in {1..30}; do /gpfs/u/home/PCPC/PCPCrgsh/scratch/ParallelComputingProject/parallel400.o >> 8Node6GPUExp/400.txt; done  
for i in {1..30}; do /gpfs/u/home/PCPC/PCPCrgsh/scratch/ParallelComputingProject/parallel800.o >> 8Node6GPUExp/800.txt; done  

