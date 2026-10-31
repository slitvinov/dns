#!/bin/bash
# build.sh SRC OUT: boxsymmp on FASRC (nvhpc 24.11, HPC-X, NCCL)
. /n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec/src/mpenv.sh
NC=$NVR/comm_libs/12.6/nccl
nvcc -O3 -gencode arch=compute_80,code=sm_80 -gencode arch=compute_90,code=sm_90 -gencode arch=compute_90,code=compute_90 -o $2 $1 -I$NVML/include -I$NC/include -I$MPI_HOME/include \
  -L$NVML/lib -L$NC/lib -L$MPI_HOME/lib -lmpi -lnccl -lcufft -Xlinker -rpath=$NVML/lib:$NC/lib:$MPI_HOME/lib
