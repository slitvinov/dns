#!/bin/bash
# q1.sh [NGPU]: worker for a multi-GPU slot: NGPU loops pop lines "TAG RE N C T NP [extra]" from queue1 (lock in $HOME)
# and run.sh them on one device each (entries of the slot CUDA_VISIBLE_DEVICES: GPU indices or MIG UUIDs)
R=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec
IFS=, read -ra DEV <<< "${CUDA_VISIBLE_DEVICES:-0,1,2,3}"
pop() { flock $HOME/yspec_q1.lock bash -c "head -n 1 $R/queue1; sed -i 1d $R/queue1"; }
w() { while :; do l=$(pop); [ -z "$l" ] && break; echo "$(date +%T) dev $1: $l"; CUDA_VISIBLE_DEVICES=$1 $R/run.sh $l; done; }
for g in $(seq 0 $((${1:-4} - 1))); do w ${DEV[$g]} & sleep 2; done; wait
echo "q1.sh: queue empty $(date)"
