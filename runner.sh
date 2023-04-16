#!/bin/bash -x

for i in {1..30}; do /gpfs/u/home/PCPC/PCPCmbdw/scratch/Parallel_SysEquation_Solver/parallel50 >> 1NodeExperimentsNew/50.txt; done  
for i in {1..30}; do /gpfs/u/home/PCPC/PCPCmbdw/scratch/Parallel_SysEquation_Solver/parallel100 >> 1NodeExperimentsNew/100.txt; done  
for i in {1..30}; do /gpfs/u/home/PCPC/PCPCmbdw/scratch/Parallel_SysEquation_Solver/parallel200 >> 1NodeExperimentsNew/200.txt; done  
for i in {1..30}; do /gpfs/u/home/PCPC/PCPCmbdw/scratch/Parallel_SysEquation_Solver/parallel400 >> 1NodeExperimentsNew/400.txt; done  
for i in {1..30}; do /gpfs/u/home/PCPC/PCPCmbdw/scratch/Parallel_SysEquation_Solver/parallel800 >> 1NodeExperimentsNew/800.txt; done
for i in {1..30}; do /gpfs/u/home/PCPC/PCPCmbdw/scratch/Parallel_SysEquation_Solver/parallel1000 >> 1NodeExperimentsNew/1000.txt; done