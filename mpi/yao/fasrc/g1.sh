#!/bin/bash
# g1.sh: L3 check, Re 20000 at 1536 then 2048 from t 0 to t* 2.6 (t 5.582), C 2, -F 0.00125*512/n, -S
Y=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec; export BIN=$Y/src/boxsymmp_b3
$Y/run.sh r20k1536 20000 1536 2.0 5.582 4 -S -F 0.000416667
$Y/run.sh r20k2048 20000 2048 2.0 5.582 4 -S -F 0.0003125
scancel $SLURM_JOB_ID
