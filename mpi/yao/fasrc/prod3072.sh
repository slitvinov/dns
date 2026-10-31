#!/bin/bash
# prod3072.sh: Yao Re 40000 production window (= G2 3072 rung): 3072^3 on 16 H200 from the 1536 k. dump at t* 1.3
# (run/l40kbase, zero-padded) to t* 2.3 (t 4.9381), C 2, -F 0.000208333, -S, spectra every 0.02 in t, diagnostics also every
# 0.002 in t* over t* 1.75-1.95 (as the rungs), snapshots (-W: |w|, w_z, lambda2/1000 half, quarter crop x < 1.5, 20.8 GB
# each) at t* 1.5 1.6 1.7 1.75 1.8 1.85 1.9 2.0 2.15 2.3; k. dumps every 0.25 in t, only the newest kept (174 GB).
# Resumable: rerun in a new 4-node slot. Starts only with yspec/gate.ok.
Y=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec; D=$Y/run/y40k3072
PY=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/venv/bin/python
g2k() { awk "NR > 1 {t = \$2} END {exit !(t * 0.4657717 > 2.0)}" $Y/run/r20k2048/out 2> /dev/null; }
until [ -e $Y/run/l40kbase/done ] && [ -e $Y/gate.ok ]; do
  [ -e $Y/gate.fail ] && exit 1
  [ ! -e $Y/gate.ok ] && [ -e $Y/run/r9kpad/done ] && [ -e $Y/run/r9k768noS/done ] && g2k && $PY $Y/src/gate.py > $Y/gate.log 2>&1
  sleep 60
done
mkdir -p $D
[ -e $D/start ] || ln -s $(ls $Y/run/l40kbase/k.[0-9]*[0-9] | sort | tail -1) $D/start
KT=$(awk "BEGIN{s=\"\"; for (t = 3.0; t < 4.93; t += 0.25) s = s (s == \"\" ? \"\" : \",\") t; print s}")
W=3.2205,3.4352,3.6499,3.7572,3.8646,3.9719,4.0793,4.2939,4.6160,4.9381
( while sleep 300; do
    te=$(awk "NR > 1 {t = \$2} END {print t * 0.4657717}" $(ls $D/out.* 2>/dev/null | sort | tail -1) 2> /dev/null)
    [ ! -e $Y/w4096.go ] && [ ! -e $Y/w4096.no ] && awk -v t=${te:-0} "BEGIN{exit !(t >= 1.95)}" && \
      grep -q exit $Y/run/l40k1536/lad.log $Y/run/l40k2048/lad.log $Y/run/l40k2560/lad.log 2> /dev/null && \
      $PY $Y/src/gate4096.py > $Y/gate4096.log 2>&1
    ls $D/k.[0-9]*[0-9] 2>/dev/null | while read f; do [ -L $f ] && continue; \
      awk -v t=$(cut -d" " -f1 $f.t) "BEGIN{exit !(t > 3.749 && t < 3.751)}" || echo $f; done | sort | head -n -1 | \
      while read f; do rm -f $f $f.t; done
  done ) &
w=$!
BIN=$Y/src/boxsymmp_b4 $Y/seg.sh $D 3072 16 40000 2.0 4.9381 $KT -S -F 0.000208333 -e 0.02 -D 0.004294,3.75720,4.18660 -W $W
s=$?; kill $w; exit $s
