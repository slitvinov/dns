#!/bin/bash
# l4.sh: multi-node timing for option A (NCCL_NCHANNELS_PER_NET_PEER=4 NCCL_MIN_NCHANNELS=16 in bench.sh), b4, -S:
# 3072^3 and 2048^3 on 16 ranks, full step and exchanges alone; then 4096^3 on 16 is skipped (does not fit). Production
# stays held (gates); the slot cancels itself. Results: yspec/bench.txt (tags A_*).
Y=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec; export BIN=$Y/src/boxsymmp_b4
hostname; scontrol show hostnames $SLURM_JOB_NODELIST | tr "\n" " "; echo
$Y/bench.sh 3072 16 20 A_full
$Y/bench.sh 3072 16 1 A_a2a BOXSYM_A2AONLY=5
$Y/bench.sh 3072 16 10 A_comp BOXSYM_NOCOMM=1
$Y/bench.sh 2048 16 20 A_full
$Y/bench.sh 2048 16 1 A_a2a BOXSYM_A2AONLY=5
scancel $SLURM_JOB_ID
