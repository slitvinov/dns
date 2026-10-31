#!/bin/bash
# spscheck.sh OUT EXPSPS PID: first-steps speed check of a multi-node run: after the run has stepped for ~4 min, its
# s/step (from the step column of OUT and wall time) must be < 3 x EXPSPS; otherwise kill PID, record the nodes in
# yspec/badnodes.txt (and the team ban list candidates) and exit 1. EXPSPS 0 disables the check.
O=$1; E=$2; P=$3; Y=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec
awk -v e=$E "BEGIN{exit !(e > 0)}" || exit 0
st() { awk "NR > 1 {s = \$1} END {print s + 0}" $O 2> /dev/null; }
until [ "$(st)" -gt 0 ] 2> /dev/null; do kill -0 $P 2> /dev/null || exit 0; sleep 10; done
s0=$(st); t0=$(date +%s); sleep 240; kill -0 $P 2> /dev/null || exit 0
s1=$(st); t1=$(date +%s); ds=$((s1 - s0))
sps=$(awk -v d=$ds -v t=$((t1 - t0)) "BEGIN{printf \"%.3f\", d > 0 ? t / d : 1e9}")
echo "$(date +%F\ %T) spscheck $O: $sps s/step over $ds steps (expected $E)" >> $Y/spscheck.log
if awk -v s=$sps -v e=$E "BEGIN{exit !(s > 3 * e)}"; then
  echo "$(date +%F\ %T) job $SLURM_JOB_ID nodes $(scontrol show hostnames $SLURM_JOB_NODELIST | tr "\n" " ") $sps s/step vs expected $E: killed" >> ${SPSLOG:-$Y/badnodes.txt}
  kill $P; exit 1
fi
