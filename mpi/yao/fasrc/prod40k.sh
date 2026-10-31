#!/bin/bash
# prod40k.sh DIR NP "N1:T1 N2:T2 ... Nk:Tk" [C [extra]]: Yao & Hussain 2020a Re 40000 to t = Tk in one slot of NP GPUs, in phases:
# phase i runs on an Ni^3 grid (quarter box) from the end of phase i - 1 (zero-padded k. restart; phase 1 from -I yao) to
# time Ti, in DIR/pi. Resumable: rerun the same command in a new slot and every phase continues from its newest k. dump.
# k. dumps every 0.25 in time inside a phase (older ones deleted as newer ones are written, the phase end ones kept), adaptive
# dt (Courant C, default 1), -S (float exchange) unless extra says otherwise. t* = 0.465772 t; t* 3 = t 6.441.
Y=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec
D=$1; NP=$2; PH=$3; C=${4:-1}; X=${5:--S}
mkdir -p $D && cd $D || exit 1
prev=""; i=0
for ph in $PH; do
  i=$((i + 1)); N=${ph%%:*}; T=${ph##*:}; P=$D/p$i
  mkdir -p $P
  if [ -n "$prev" ] && [ ! -e $P/start ] && [ -z "$(ls $P/k.[0-9]*[0-9] 2>/dev/null)" ]; then
    k=$(ls $prev/k.[0-9]*[0-9] | sort | tail -1); ln -s $k $P/start
  fi
  done_=$(cat $P/done 2>/dev/null)
  if [ "$done_" != "$T" ]; then
    t0=$(cat $(ls $P/k.[0-9]*[0-9].t 2>/dev/null | sort | tail -1) 2>/dev/null | cut -d" " -f1); t0=${t0:-0}
    KT=$(awk -v a=$t0 -v b=$T "BEGIN{s=\"\"; for (t = 0.25 * (int(a / 0.25) + 1); t < b - 1e-9; t += 0.25) s = s (s == \"\" ? \"\" : \",\") t; print s == \"\" ? \"-\" : s}")
    # keep only the newest two k. dumps of this phase while it runs (each 4096^3 dump is 413 GB)
    ( while sleep 600; do ls $P/k.[0-9]*[0-9] 2>/dev/null | sort | head -n -2 | while read f; do [ -L $f ] || rm -f $f; done; done ) &
    w=$!
    $Y/seg.sh $P $N $NP 40000 $C $T "$KT" $X || { kill $w; exit 1; }
    kill $w
    echo $T > $P/done
  fi
  prev=$P
done
echo "prod40k.sh: done $(date)"
