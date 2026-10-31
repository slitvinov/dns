#!/bin/bash
# bone.sh TAG RE BASEDIR T1 KT: one-node 2560^3 run (4 H200, b6, -S, BOXSYM_NL=0: ~145 GB per rank) from the 1536 k. dump in
# BASEDIR (t* 1.3), zero-padded, to time T1; C 2, -F 0.00025, spectra every 0.02, diagnostics every 0.002 in t* over
# 1.75-1.95; k. dumps at KT and T1 (newest kept). Resumable (seg2.sh). Cancels the slot at the end.
Y=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec; TAG=$1; RE=$2; BD=$3; T1=$4; KT=$5; D=$Y/run/$TAG
mkdir -p $D; [ -e $D/start ] || ln -s $(ls $BD/k.[0-9]*[0-9] | sort | tail -1) $D/start
( while sleep 600; do ls $D/k.[0-9]*[0-9] 2>/dev/null | while read f; do [ -L $f ] || echo $f; done | sort | head -n -1 | while read f; do rm -f $f $f.t; done; done ) &
w=$!
export BOXSYM_NL=0; BIN=$Y/src/boxsymmp_b6 $Y/seg2.sh $D 2560 4 $RE 2.0 $T1 $KT -S -F 0.00025 -e 0.02 -D 0.004294,3.75720,4.18660
s=$?; kill $w; [ $s = 0 ] && touch $D/finished; scancel $SLURM_JOB_ID
