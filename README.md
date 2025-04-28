<h2>Build</h2>

Compile, prepear initial conditions and run. Outputs step, time,
energy, and enstropy.
<pre>
$ c99 main.c -fopenmp -O3 -march=native -lfftw3 -lfftw3_omp -lm -o dns
$ ./tgv.py -l 6 -o tgv.raw
$ tgv.py: n=64
$ ./dns -i tgv.raw -t 10 -n 0.01 -s 0.01
dns: omp_get_max_threads: 8
dns: n = 64
         0  0.0000000000000000e+00  6.2500000000000000e-02  1.8750000000000000e-01
        10  9.9999999999999992e-02  6.2109882886321448e-02  1.8653150768334470e-01
        20  2.0000000000000004e-01  6.1689599926948888e-02  1.8586250548841311e-01
        30  3.0000000000000010e-01  6.1239137609790719e-02  1.8547074698875845e-01
        40  4.0000000000000019e-01  6.0758404926193962e-02  1.8533236230172120e-01
        50  5.0000000000000022e-01  6.0247339559280366e-02  1.8542325605932955e-01
        ...
</pre>

```
Usage: dns [-v] [-d] -i <input.raw> -n <viscosity> -t <end time> -s <time step>

Options:
  -i <input.raw>    Input file
  -n <viscosity>    Viscosity
  -t <end time>     End time
  -s <time step>    Time step
  -v                Verbose output
  -d                Dump snapshots
  -h                Show this help message

Example:
  dns -i tgv.raw -n 0.01 -t 1.0 -s 0.001 -v
```

<h3>Validataion</h2>

<p align="center"><img src="img/tgv.svg" width=600></p>
Figure: Energy dissipation rate vs. time for the Taylor–Green
vortex. Reference data (points) from Brachet et al. is shown alongside
simulation results (lines). From top to bottom at time = 0: Re = 100,
200, 400, 800, 1600, 3000.

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
