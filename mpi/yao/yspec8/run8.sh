#!/bin/bash
# run8.sh: production-path check on 8 GPUs at 2560^3 (Re 40000, -I yao, -S, C 2, -F 0.00025, NCCL 8/32):
#  A: ~20 steps with a k. dump after ~10 and a snapshot (-W);  B: a fresh mpirun restarts from the dump to the same end;
#  rows compared; memory per rank; NCCL init lines. Summary in res.txt.
cd $(dirname $0); R=$PWD/res.txt; : > $R
NP=${NP:-8}; N=${N:-2560}
mp() { mpirun -np $NP --bind-to none --allow-run-as-root -x NCCL_NCHANNELS_PER_NET_PEER=8 -x NCCL_MIN_NCHANNELS=32 -x NCCL_DEBUG=INFO -x NCCL_DEBUG_SUBSYS=INIT "$@"; }
G="-x $N -y $N -z $N -X 2 -Y 2 -Z 2 -n 1.920096311e-4 -s 0.01 -c 2 -F 0.00025 -f -S"
echo "$(date) $(hostname) np $NP n $N; nvidia-smi: $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader | sort | uniq -c | tr '\n' ' ')" >> $R
df -h . | tail -1 >> $R
rm -rf A B; mkdir -p A B
(cd A && mp ../boxsymmp $G -k -e 1 -I yao -T 0.005 -W 0.007 -t 0.0098 > out 2> err); echo "A exit $?" >> $R
grep -h "rank [0-9]*/\|s/step\|wrote\|snapshot" A/err | cut -c1-200 >> $R
grep -h "NCCL INFO.*\(version\|Using network\|NET/\|comm .* nRanks\)" A/err A/out 2>/dev/null | sort | uniq -c | head -8 >> $R
k=$(cd A && ls k.[0-9]*[0-9] | grep -v k.00000000 | sort | tail -1)
echo "dump $k $(stat -c %s A/$k) B (expected $((3 * (N / 2 + 1) ** 3 * 16)))" >> $R
(cd B && ln -s ../A/$k $k && ln -s ../A/$k.t $k.t && mp ../boxsymmp $G -k -e 1 -r $k -t 0.0098 > out 2> err); echo "B exit $?" >> $R
grep -h "read\|s/step" B/err >> $R
s=${k#k.}; s=$((10#$s))
awk 'NR == FNR {if (FNR > 1) a[$1] = $0; next} FNR > 1 && $1 > s && ($1 in a) {n++; if (a[$1] != $0) {d++; print "A: " a[$1]; print "B: " $0}} END {printf "restart: %d common rows after step %d, %d differ\n", n, s, d + 0}' s=$s A/out B/out >> $R
ls -la A/w.* 2>/dev/null | awk '{print $5, $9}' >> $R
rm -f A/k.[0-9]* B/k.[0-9]*
cat $R
