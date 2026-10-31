#!/bin/bash
# l4b.sh: second fine-phase chunk (resumes prod3072.sh from the newest k. dump), then cancel the slot
Y=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec
$Y/prod3072.sh
scancel $SLURM_JOB_ID
