#!/bin/bash
# qbench.sh NP N...: short multi-node timing (b3): -S, exchanges alone (-S), compute alone, and the old path (b2 switches)
R=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec; NP=$1; shift; K=${K:-10}
hostname; scontrol show hostnames $SLURM_JOB_NODELIST | tr "\n" " "; echo
for N in "$@"; do
  BIN=$R/src/boxsymmp_b3 $R/bench.sh $N $NP $K q_newS_b3 -- -S
  BIN=$R/src/boxsymmp_b3 $R/bench.sh $N $NP 1 q_a2anewS BOXSYM_A2AONLY=3 -- -S
  BIN=$R/src/boxsymmp_b3 $R/bench.sh $N $NP $K q_comp_b3 BOXSYM_NOCOMM=1
  BIN=$R/src/boxsymmp_b2 $R/bench.sh $N $NP $K q_old BOXSYM_SYN6=1 BOXSYM_NOFUSE=1 BOXSYM_SELFCOPY=1 BOXSYM_NL=0
done
scancel $SLURM_JOB_ID
