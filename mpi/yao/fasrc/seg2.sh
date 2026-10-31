#!/bin/bash
# seg.sh DIR N NP RE C T1 KT [boxsymmp options...]: one resumable segment of a Yao & Hussain run in DIR (2 pi cube,
# quarter box, Hou-Li filter) on an N^3 grid with NP ranks: starts from the newest k.STEP in DIR (its grid from its size;
# zero-padded if smaller than N), else from the file DIR/start (a k. file, e.g. a link to another run), else -I yao.
# Runs to time T1 with adaptive dt (-c C), k. dumps at the comma list KT (times; "-" none) and at T1, spectra every 0.1.
# Out/err per segment: out.STEP, err.STEP. Extra options (e.g. -S) are passed to boxsymmp.
. /etc/profile > /dev/null 2>&1
Y=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec; S=$Y/src
. $S/mpenv.sh
D=$1; N=$2; NP=$3; RE=$4; C=$5; T1=$6; KT=$7; shift 7; X="$@"
B=${BIN:-$S/boxsymmp}
NU=$(awk -v r=$RE "BEGIN{printf \"%.17g\", 7.680385244/r}")
mkdir -p $D && cd $D || exit 1
grid() { awk -v s=$(stat -L -c %s $1) "BEGIN{m = s / 48; c = exp(log(m) / 3); printf \"%d\", 2 * (int(c + 0.5) - 1)}"; }
k=$(ls k.[0-9]*[0-9] 2>/dev/null | sort | tail -1)
[ -z "$k" ] && [ -e start ] && k=start
if [ -n "$k" ]; then
  n0=$(grid $k); step=$(basename $(readlink -f $k) | sed "s/^k\.//")
  [ "$k" = start ] && { ln -sf $(readlink -f start) k.$step; ln -sf $(readlink -f start).t k.$step.t; k=k.$step; }
  IN="-r $k"; [ "$n0" != "$N" ] && IN="$IN -G $n0,$n0,$n0"
else
  IN="-I yao"; step=0
fi
T="-T $T1"; [ "$KT" != "-" ] && T="-T $KT,$T1"
echo "$(date +%F\ %T) seg.sh: $D n $N np $NP Re $RE nu $NU C $C to t $T1 from [$IN] dumps [$T] [$X] host $(hostname)" | tee -a seg.log
mpirun -np $NP --oversubscribe --bind-to none -x NCCL_NCHANNELS_PER_NET_PEER=4 -x NCCL_MIN_NCHANNELS=16 $B -x $N -y $N -z $N -X 2 -Y 2 -Z 2 -n $NU -s 0.01 -c $C -f -k -e 0.1 $IN -t $T1 $T $X \
  > out.$step 2> err.$step &
p=$!; $Y/src/spscheck.sh out.$step ${EXPSPS:-0} $p &
wait $p; s=$?; echo "$(date +%F\ %T) seg.sh: exit $s; $(grep -h s/step err.$step)" | tee -a seg.log; exit $s
