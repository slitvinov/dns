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
vortex. Reference data (points) from Brachet et al., `fourier` at
256^3 (black lines) and `tg -M 128` (red dashed lines). From top to
bottom at time = 0: Re = 100, 200, 400, 800, 1600, 3000.

`brachet.py` reproduces more of the paper from `inv042/`, `inv084/`
(inviscid runs, k_max = 42 and 84) and `tg0256/`:

- figures 3–4: inviscid spectra, `img/spectrum.svg`;
- table 1 and figure 5: the width of the analyticity strip δ(t),
  `img/delta.svg`; δ(1.5) = 0.186 (paper 0.192), δ(2.5) = 0.031
  (0.034) at k_max = 84;
- figure 12: skewness S3(0)(t), `img/skewness.svg`. The figure of the
  paper agrees with the isotropic relation (5.8), not with the average
  of (∂vx/∂x)^3 in (5.6).

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
