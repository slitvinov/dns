#!/bin/bash
# build.sh: boxsymmp (b6) for x86 H100 (sm_90) with the system OpenMPI and NCCL (Ubuntu packages)
set -e
cd $(dirname $0)
NVCC=$(command -v nvcc || echo /usr/local/cuda/bin/nvcc)
MI=$(mpicc --showme:incdirs | tr ' ' '\n' | sed 's/^/-I/' | tr '\n' ' ')
ML=$(mpicc --showme:libdirs | tr ' ' '\n' | sed 's/^/-L/' | tr '\n' ' ')
$NVCC -O3 -arch=sm_90 -o boxsymmp boxsymmp.cu $MI $ML -lmpi -lnccl -lcufft
ls -la boxsymmp
