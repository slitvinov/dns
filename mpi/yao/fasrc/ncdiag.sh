#!/bin/bash
# ncdiag.sh NP: NCCL multi-node diagnosis for boxsymmp (exchanges alone, 2048^3, float, b3), variants of the transport env;
# lines in yspec/bench.txt (tags nd_*) and the NCCL INFO transport lines in yspec/ncdiag.<jobid>.txt; cancels the slot.
R=/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec; NP=$1; O=$R/ncdiag.$SLURM_JOB_ID.txt
export BIN=$R/src/boxsymmp_b3
echo "$(date) nodes $(scontrol show hostnames $SLURM_JOB_NODELIST | tr '\n' ' ') np $NP" > $O
for n in $(scontrol show hostnames $SLURM_JOB_NODELIST | head -2); do ssh $n "hostname; ibstat 2>/dev/null | grep -E 'CA |State|Rate|Link layer'" >> $O 2>&1; done
v() { tag=$1; shift; $R/bench.sh 2048 $NP 1 nd_$tag BOXSYM_A2AONLY=3 NCCL_DEBUG=INFO NCCL_DEBUG_SUBSYS=INIT,NET "$@" -- -S
      echo "== $tag [$*]: $(grep "nd_$tag n 2048 np $NP" $R/bench.txt | tail -1 | cut -d'|' -f2)" >> $O
      grep -h "NET/\|via\|GDR\|Using network\|NCCL version" $R/bench/nd_$tag.2048.$NP/err $R/bench/nd_$tag.2048.$NP/out 2>/dev/null | sort | uniq -c | sort -rn | head -12 >> $O; }
v default
v hca NCCL_IB_HCA=mlx5_2,mlx5_3,mlx5_4,mlx5_5
v noplugin NCCL_NET_PLUGIN=none NCCL_IB_HCA=mlx5_2,mlx5_3,mlx5_4,mlx5_5
v gdr NCCL_IB_HCA=mlx5_2,mlx5_3,mlx5_4,mlx5_5 NCCL_NET_GDR_LEVEL=SYS NCCL_IB_QPS_PER_CONNECTION=4
v mpi BOXSYM_MPI=1
scancel $SLURM_JOB_ID
