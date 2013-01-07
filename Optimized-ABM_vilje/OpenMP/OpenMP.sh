#!/bin/bash
#PBS -A nn2849k
#PBS -N app
#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=32:ompthreads=16
 
# Tips to improving performance
#
# 1. Adjusting MPI_BUFS_PER_PROC and MPI_BUFS_PER_HOST.
#
# Use the "-stats" option to mpiexec_mpt to get additional information in the
# output file. Included in that information is the number of retries for
# allocating MPI buffers. After you have executed your program with the "-stats"
# option, you can see this by typing something similar to:
#
# $ cat my_mpi_job.o55809 | grep retries | grep -v " 0 retries"
#
# You can then increase the values of MPI_BUFS_PER_PROC (default 32) and
# MPI_BUFS_PER_HOST (default 96) until the number of retries is sufficiently
# low, e.g. by uncommenting these lines:
#
# export MPI_BUFS_PER_PROC=256
# export MPI_BUFS_PER_HOST=1024
#
# See "man mpi" for more information.
#
#
# 2. Adjusting MPI_BUFFER_MAX
#
# For some codes it gives a significant increase in performance to specify a
# value for MPI_BUFFER_MAX. According to "man mpi" this value "Specifies a
# minimum message size, in bytes, for which the message will be considered a
# candidate for single-copy transfer." The value of MPI_BUFFER_MAX varies from
# program to program, but typical values are between 2048 and 32768. You can
# therefore test if this improves the performance of your program by executing
# it like this:
#
# export MPI_BUFFER_MAX=2048
# time -p mpiexec_mpt -n 1024 ./myprog
#
# See "man mpi" for more information.
 
module load intelcomp/11.1.073
module load mpt
 
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=4
export KMP_AFFINITY=disabled
#omplace -vv -c 0,15 ./app 0.1 200000
omplace -vv -c 0,1,11,12  ./app 0.1 200000
#omplace -vv -c 0,1,2,3,10,11,12,13   ./app  0.1 200000
#omplace -vv -c 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15  ./app 0.1 200000
