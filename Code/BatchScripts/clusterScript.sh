#!/bin/zsh

### jobname
#BSUB -J Test

### output file
#BSUB -o testOut.%J

### time constraint
#BSUB -W 0:05

### memory constraint
#BSUB -M 4096

### send email when finished
#BSUB -N

### exclusive node
#BSUB -x

### openmp
#BSUB -a openmp

### number of cores
#BSUB -n 4

cd $HOME/HPMC/MatrixChainProduct/Code/

time ./main.x

