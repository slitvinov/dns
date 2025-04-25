<h2>Build</h2>

Serial
<pre>
$ ./tgv.py -l 6 -o tgv.raw
$ c99 main.c -O2 -g -lfftw3 -lm
$ ./a.out -i tgv.raw -t 10 -n 0.01
         0  0.0000000000000000e+00  1.2500000000000000e-01  3.7500000000000000e-01
        10  9.9999999999999992e-02  1.2425198809829003e-01  3.7314144557524406e-01
        20  2.0000000000000004e-01  1.2350692709732673e-01  3.7204458637899375e-01
        30  3.0000000000000010e-01  1.2276332514957790e-01  3.7168281219301155e-01
        40  4.0000000000000019e-01  1.2201974653636426e-01  3.7203557060776904e-01
        50  5.0000000000000022e-01  1.2127480635967738e-01  3.7308729039195336e-01
        60  6.0000000000000031e-01  1.2052716417583942e-01  3.7482625233875647e-01
        70  7.0000000000000040e-01  1.1977551493310098e-01  3.7724338978315803e-01
        ...
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

<h3>Validataion</h2>

<p align="center"><img src="img/tgv.svg" width=600></p>
Figure: Energy dissipation rate vs. time for the Taylor–Green
vortex. Reference data (points) from Brachet et al. is shown alongside
simulation results (lines). From top to bottom at time = 0: Re = 3000,
1600, 800, 400, 200, 100.

<h2>References</h2>

- Brachet, M. E., Meiron, D. I., Orszag, S. A., Nickel, B. G., Morf,
  R. H., & Frisch, U. (1983). Small-scale structure of the
  Taylor–Green vortex. Journal of Fluid Mechanics, 130, 411-452.

- Orszag, S. A., & Patterson Jr, G. S. (1972). Numerical simulation of
  three-dimensional homogeneous isotropic turbulence. Physical review
  letters, 28(2), 76.

- Mortensen, M. (2016). Massively parallel implementation in Python of
  a pseudo-spectral DNS code for turbulent flows. arXiv preprint
  arXiv:1607.00850.
