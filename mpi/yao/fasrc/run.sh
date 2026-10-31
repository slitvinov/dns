#!/bin/bash
# run.sh TAG RE N C T [NP [extra boxsymmp options]]: Yao & Hussain 2020a tubes, boxsymmp quarter box (2 pi cube), IC -I yao,
# adaptive dt (-c C, dt <= 0.01), Hou-Li filter, spectra every 0.25, run dir yspec/run/TAG. NP ranks (default 1).
. /etc/profile > /dev/null 2>&1
R=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec; S=$R/src
. $S/mpenv.sh
TAG=$1; RE=$2; N=$3; C=$4; T=$5; NP=${6:-1}; shift 6 2>/dev/null; X="$@"
NU=$(awk -v r=$RE "BEGIN{printf \"%.17g\", 7.680385244/r}")
D=$R/run/$TAG; mkdir -p $D && cd $D || exit 1
echo "run.sh: $TAG Re $RE nu $NU n $N C $C T $T np $NP extra [$X] host $(hostname) gpus $CUDA_VISIBLE_DEVICES $(date)" | tee -a run.info
B=${BIN:-$S/boxsymmp}
mpirun -np $NP --oversubscribe --bind-to none $B -x $N -y $N -z $N -X 2 -Y 2 -Z 2 -n $NU -s 0.01 -c $C -f -I yao -t $T -e 0.25 $X > out 2> err
s=$?; echo "run.sh: exit $s $(date)" | tee -a run.info; tail -2 err
