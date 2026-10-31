#!/bin/bash
# bmat.sh NP N... : timing matrix in a slot with NP GPUs, per grid N:
#   newSnc = newS with NCCL_NCHANNELS_PER_NET_PEER=4 NCCL_P2P_NET_CHUNKSIZE=512k (multi-node only)
#   newShc2/4 = newS with chunked stage boundaries (BOXSYM_HC)
#   old   = the old code path (6-exchange synthesis, own diagnostic synthesis, local copy in the exchange, no nl_a/nl_b)
#   new   = 5 exchanges + nl_a/nl_b + fused diagnostics + direct own section;  newS = new -S;  newSnl0 = newS without nl
#   a2aold / a2anewS = exchanges alone (45 double / 40 float per step);  comp = new without exchanges
R=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec
NP=$1; shift; K=${K:-20}
OLD="BOXSYM_SYN6=1 BOXSYM_NOFUSE=1 BOXSYM_SELFCOPY=1 BOXSYM_NL=0"
for N in "$@"; do
  $R/bench.sh $N $NP $K old $OLD
  $R/bench.sh $N $NP $K new
  $R/bench.sh $N $NP $K newS -- -S
  $R/bench.sh $N $NP $K newSnl0 BOXSYM_NL=0 -- -S
  $R/bench.sh $N $NP $K newShc2 BOXSYM_HC=2 -- -S
  $R/bench.sh $N $NP $K newShc4 BOXSYM_HC=4 -- -S
  [ $NP -gt 4 ] && $R/bench.sh $N $NP $K newSnc NCCL_NCHANNELS_PER_NET_PEER=4 NCCL_P2P_NET_CHUNKSIZE=524288 -- -S
  $R/bench.sh $N $NP 1 a2aold $OLD BOXSYM_A2AONLY=5
  $R/bench.sh $N $NP 1 a2anewS BOXSYM_A2AONLY=5 -- -S
  $R/bench.sh $N $NP $K comp BOXSYM_NOCOMM=1
done
