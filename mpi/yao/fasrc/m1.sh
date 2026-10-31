#!/bin/bash
# m1.sh: option B measurement: 2560^3 on one node (4 H200), b5, memory per rank and s/step with the diet
# (BOXSYM_NL=0, BOXSYM_RB=2, double exchange) and variants; then cancel the slot
Y=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec; export BIN=$Y/src/boxsymmp_b5
hostname
$Y/bench.sh 2560 4 10 B_nl0rb2 BOXSYM_NL=0 BOXSYM_RB=2
$Y/bench.sh 2560 4 10 B_nl0rb3 BOXSYM_NL=0 BOXSYM_RB=3
$Y/bench.sh 2560 4 10 B_rb2 BOXSYM_RB=2
$Y/bench.sh 2048 4 10 B_ref2048_nl0rb2 BOXSYM_NL=0 BOXSYM_RB=2
$Y/bench.sh 2048 4 10 B_ref2048 BOXSYM_RB=1
BIN=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/mlc/boxsymmp_b5_c129 $Y/bench.sh 2560 4 10 B_c129_nl0rb2 BOXSYM_NL=0 BOXSYM_RB=2
BIN=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/mlc/boxsymmp_b5_c129 $Y/bench.sh 2560 4 10 B_c129_nl0 BOXSYM_NL=0
grep -h "rank 0/" $Y/bench/B_*.2560.4/err | sed "s/^/  /"
scancel $SLURM_JOB_ID
