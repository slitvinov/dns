#!/bin/bash
# l4.sh: G2 rung 2560 (to t* 1.95), then the 3072 production window (to t* 2.3); cancel the slot at the end
Y=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec
[ -e $Y/run/l40k2560/lad.log ] && grep -q exit $Y/run/l40k2560/lad.log || EXPSPS=1.2 $Y/lad.sh 2560 16
$Y/prod3072.sh
scancel $SLURM_JOB_ID
