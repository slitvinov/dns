#!/bin/bash
# bslot.sh NP N...: in a benchmark slot: bmat.sh (boxsymmp_b2) then bmat3.sh (boxsymmp_b3), then cancel the slot
R=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec; NP=$1; shift
hostname; scontrol show hostnames $SLURM_JOB_NODELIST | tr "\n" " "; echo
BIN=$R/src/boxsymmp_b2 $R/bmat.sh $NP "$@"
$R/bmat3.sh $NP "$@"
scancel $SLURM_JOB_ID
