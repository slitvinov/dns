#!/bin/bash
# g2base.sh: Re 40000 at 1536 from t 0 to t* 1.3 (t 2.7911) with the k. dump (run/l40kbase), then the 2048 rung
Y=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec
BIN=$Y/src/boxsymmp_b4 $Y/seg.sh $Y/run/l40kbase 1536 4 40000 2.0 2.7911 - -S -F 0.000416667 && touch $Y/run/l40kbase/done
$Y/lad.sh 2048 4
$Y/lad.sh 1536 4
scancel $SLURM_JOB_ID
