python tgv.py -o tgv.raw -l 8
cat <<! | awk '{printf "%04d %.16e\n", $1, 1/$1}' | \
    xargs --process-slot-var I -n 2 -P `nproc` sh -xc \
          '
sleep $I
i=$((64+10*I))
j=$((64+10*(I+1)-1))
exec taskset --cpu-list $i-$j ./dns -v -t 10 -n $1 -s 0.0025 -i tgv.raw > $0
'
100
200
400
800
1600
3000
!
