<h2>Build</h2>

Serial
<pre>
$ ./tgv.py -l 4 -o tgv.raw
$ c99 main.c -lfftw3_mpi -lfftw3 -O2 -g -lm
$ ./a.out -i tgv.raw -t 1.0 -n 1.0
dns:        0  0.0000e+00  1.2500e-01
dns:       10  1.0000e-01  6.8593e-02
dns:       20  2.0000e-01  3.7629e-02
dns:       30  3.0000e-01  2.0640e-02
dns:       40  4.0000e-01  1.1322e-02
dns:       50  5.0000e-01  6.2114e-03
dns:       60  6.0000e-01  3.4080e-03
dns:       70  7.0000e-01  1.8701e-03
dns:       80  8.0000e-01  1.0262e-03
dns:       90  9.0000e-01  5.6315e-04
dns:      100  1.0000e+00  3.0905e-04
</pre>

with MPI
<pre>
$ mpicc mpi/main.c -lfftw3_mpi -lfftw3 -O2 -g -lm
$ mpiexec ./a.out
eng: 1.2499062522321387e-01
eng: 1.2498125036872199e-01
eng: 1.2497187465530142e-01
eng: 1.2496249730189801e-01
eng: 1.2495311752707980e-01
</pre>
