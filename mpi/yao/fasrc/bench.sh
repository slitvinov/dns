#!/bin/bash
# bench.sh N NP STEPS TAG [VAR=val ...] [-- boxsymmp options]: timing of one variant of the yao-spec boxsymmp on an N^3 grid
# (quarter box, -I yao, Re 40000 nu, fixed dt 1e-4, Hou-Li), NP ranks; appends "TAG N NP s/step mem" to bench.txt
. /etc/profile > /dev/null 2>&1
R=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec; S=$R/src
. $S/mpenv.sh
N=$1; NP=$2; K=$3; TAG=$4; shift 4
E=""; while [ $# -gt 0 ] && [ "$1" != "--" ]; do export "$1"; E="$E -x ${1%%=*}"; shift; done; [ "$1" = "--" ] && shift
B=${BIN:-$S/boxsymmp}
D=$R/bench/$TAG.$N.$NP; mkdir -p $D && cd $D || exit 1
T=$(awk -v k=$K "BEGIN{printf \"%.6g\", k * 1e-4}")
/usr/bin/time -f "%e s wall" mpirun -np $NP --oversubscribe --bind-to none $E -x NCCL_DEBUG=WARN $B -x $N -y $N -z $N -X 2 -Y 2 -Z 2 \
  -n 1.920096311e-4 -s 1e-4 -t $T -f -I yao "$@" > out 2> err
s=$?
sp=$(grep -h "s/step" err | tail -1 | sed "s/.*: //"); mem=$(grep -h "rank 0/" err | sed "s/.*gpu used //")
echo "$TAG n $N np $NP exit $s: $sp | $(grep -h "exchanges only" err | sed "s/.*): //") | rank0 mem $mem | $(hostname) $(date +%T)" | tee -a $R/bench.txt
