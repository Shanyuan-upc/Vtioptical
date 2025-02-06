#!/bin/bash


cc=mpiicc
exe=acmod2d

$cc  -qopenmp -g  -o ../bin/$exe  fd2d_mod_staggered_omp_mpi_cpml.c -lm


