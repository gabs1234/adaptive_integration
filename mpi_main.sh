#!/bin/bash
PROGRAM=bin/mpi_main

NPROC=$1
A=$2
B=$3
PREC=$4
ORDER=$5

if [ $# -lt 5 ]; then
    echo "usage: ./main.sh <A> <B> <PREC> <ORDER>"
    echo "A,B: integration interval"
    echo "PREC: number of decimels of precision"
    echo "ORDER: order of gauss legendre quadrature"
    exit 1
fi

make mpi_main

mpirun -np $NPROC $PROGRAM $A $B $PREC $ORDER