#!/bin/bash

make mpi_bench
rm data/para_bench.dat

for i in {1..4}
do
	mpirun --mca opal_warn_on_missing_libcuda 0 -np $i bin/mpi_bench 0 0 12 12 >> data/para_bench.dat
done
