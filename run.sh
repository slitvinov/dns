cat <<! | awk '{printf "%04d %.16e\n", $1, 1/$1}' | \
    xargs --process-slot-var I -n 2 -P `nproc` sh -xc \
          'exec taskset --cpu-list $I ./dns -t 10 -n $1 -s 0.1 -i tgv.raw > $0'
100
200
400
800
1600
3000
!
