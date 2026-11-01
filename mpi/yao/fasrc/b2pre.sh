#!/bin/bash
# b2pre.sh: the 4096 start field (Sergey 2026-10-31: start y40k4096 at t* 1.5): 2560^3 from l40kbase (1536, t* 1.3) zero-padded,
# to t* 1.5 exactly (t = 1.5 / 0.4657717102700246 = 3.220466142) with a k. dump there; flags identical to B2 (bone.sh):
# b6, -S, C 2, -F 0.00025, BOXSYM_NL=0, nu = Gamma0 / 40000. On success: verifies the dump size, touches finished, frees the slot.
Y=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec; D=$Y/run/l40k2560s
mkdir -p $D; [ -e $D/start ] || ln -s $Y/run/l40kbase/k.00003549 $D/start
export BOXSYM_NL=0; BIN=$Y/src/boxsymmp_b6 $Y/seg3.sh $D 2560 4 40000 2.0 3.220466142 - -S -F 0.00025 -e 0.02
s=$?; k=$(ls $D/k.[0-9]*[0-9] 2>/dev/null | while read f; do [ -L $f ] || echo $f; done | sort | tail -1)
if [ $s = 0 ] && [ -n "$k" ] && [ -s $k.t ] && [ "$(stat -c %s $k)" = "$((3 * 1281 * 1281 * 1281 * 16))" ]; then
  echo "$(date +%F\ %T) b2pre: start dump $k ($(cat $k.t)) verified" >> $D/seg.log; touch $D/finished
else
  echo "$(date +%F\ %T) b2pre: FAILED (exit $s, dump [$k])" >> $D/seg.log
fi
scancel $SLURM_JOB_ID
