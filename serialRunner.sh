#!/bin/bash -x

for i in {1..30}; do /gpfs/u/home/PCPC/PCPCmbdw/scratch/Parallel_SysEquation_Solver/serial_lu50 >> serialExperiments/50.txt; done  
for i in {1..30}; do /gpfs/u/home/PCPC/PCPCmbdw/Parallel_SysEquation_Solver/serial_lu100 >> serialExperiments/100.txt; done  
for i in {1..30}; do /gpfs/u/home/PCPC/PCPCmbdw/scratch/Parallel_SysEquation_Solver/serial_lu200 >> serialExperiments/200.txt; done  
for i in {1..30}; do /gpfs/u/home/PCPC/PCPCmbdw/scratch/Parallel_SysEquation_Solver/serial_lu400 >> serialExperiments/400.txt; done  
for i in {1..30}; do /gpfs/u/home/PCPC/PCPCmbdw/scratch/Parallel_SysEquation_Solver/serial_lu800 >> serialExperiments/800.txt; done