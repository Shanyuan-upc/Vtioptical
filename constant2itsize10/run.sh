#!/bin/bash

srcdir=../src/
bindir=../bin/
curdir=`pwd`



exe=../bin/acmod2d

cd $srcdir
sh make.sh
cd $curdir

rm snapshot.dat

$exe

#mpirun -n 10  ${exe}
