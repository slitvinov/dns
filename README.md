<h2>Build</h2>

Two pseudo-spectral solvers for incompressible flow in a periodic box.

- `fourier.c`: general complex Fourier series, rotational form u × ω,
  2/3 dealiasing, RK4. Any initial condition.
- `tg.c`: the Taylor–Green vortex with the symmetric representation of
  Brachet et al. (1983): sine/cosine series, all-even and all-odd modes
  kept apart, fundamental box [0, π/2]^3, leapfrog for the nonlinear
  term and Crank–Nicolson for the viscous term.

Needs FFTW ≥ 3.3.9 with OpenMP (`fftw_threads_set_callback`).
<pre>
$ make
$ ./tgv.py -l 6 -o tgv.raw
tgv.py: n=64
$ ./fourier -i tgv.raw -t 10 -n 0.01 -s 0.01
         0  0.0000000000000000e+00  1.2500000000000000e-01  3.7500000000000000e-01
        10  9.9999999999999992e-02  1.2425198808904905e-01  3.7314144554701101e-01
        ...
$ ./tg -M 32 -t 10 -n 0.01 -s 0.01
</pre>

Both print step, time, energy ½⟨|u|²⟩ and enstrophy ½⟨|ω|²⟩. `tg -M
128` is the (256)^3 run of the paper. `tg` also prints the skewness and
flatness factors S3..S8 of ∂vx/∂x, S̄4, S̄6 of ∂²vx/∂x², and the
palinstrophy, and with `-e <interval>` writes the energy spectrum
E(k) in bins of width 1/2 to `e.<step>`.

```
Usage: fourier [-v] [-d] -i <input.raw> -n <viscosity> -t <end time> -s <time step>
Usage: tg -M <modes> -n <viscosity> -t <end time> -s <time step> [-e <spectrum interval>]
```

<h3>Validation</h2>

<p align="center"><img src="img/tgv.svg" width=600></p>
Figure: Energy dissipation rate vs. time for the Taylor–Green
vortex, `tg` at 64^3, 128^3, 256^3 and 512^3 (lines, `tgv.gp`, runs in
`data/tg/<grid>/<Re>`), and figure 7 of Brachet et al. (circles,
re-digitized in `img/ref.txt`: t, ε, Re, 1 or 0 where curves cross;
crosses mark where curves cross). The circles take the colour of the
grid the paper used. At the paper's grid the peak agrees to 0.3% (0.9% at Re 800):

| Re | grid | tg | figure 7 | 512^3 |
|---|---|---|---|---|
| 400 | 128^3 | 0.01098 | 0.01100 | 0.01098 |
| 800 | 128^3 (inferred) | 0.01188 | 0.01199 | 0.01172 |
| 1600 | 256^3 | 0.01292 | 0.01296 | 0.01286 |
| 3000 | 256^3 | 0.01530 | 0.01528 | 0.01505 |

The agreement at Re 3000 needs the truncation of the paper: the even
and the odd modes are kept by array index, n <= M/3, so at M = 128
even wavenumbers go to 84 and odd ones to 85. Cutting both at 84
gives 0.01472.

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
