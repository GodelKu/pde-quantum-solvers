### stationary_harmonic

| Method | Nominal order | Fitted order | Last-pair order | Error, smallest N | Error, largest N | Time, largest N (ms) |
|---|---|---|---|---|---|---|
| fd2 | 2 | 2.01 | 2.00 | 1.24e-01 (N=40) | 1.22e-04 (N=1280) | 3.3 |
| fd4 | 4 | 3.98 | 4.00 | 1.41e-02 (N=40) | 1.60e-08 (N=1280) | 463.2 |
| numerov | 4 | 4.01 | 4.00 | 5.98e-03 (N=40) | 5.99e-09 (N=1280) | 5004.6 |
| fourier | spectral | - | - | 1.00e-01 (N=16) | 2.66e-14 (N=96) | 28.7 |
| chebyshev | spectral | - | - | 6.47e-01 (N=16) | 2.49e-14 (N=96) | 47.1 |
| shooting | 4 | 4.00 | 4.00 | 9.90e-06 (N=200) | 2.45e-09 (N=1600) | 1102.7 |

### stationary_double_well

| Method | Nominal order | Fitted order | Last-pair order | Error, smallest N | Error, largest N | Time, largest N (ms) |
|---|---|---|---|---|---|---|
| fd2 | 2 | 2.01 | 2.00 | 8.07e-02 (N=40) | 8.06e-05 (N=1280) | 2.1 |
| fd4 | 4 | 3.98 | 4.00 | 7.32e-03 (N=40) | 8.30e-09 (N=1280) | 1691.5 |
| numerov | 4 | 4.01 | 4.00 | 3.10e-03 (N=40) | 3.11e-09 (N=1280) | 549.3 |
| fourier | spectral | - | - | 1.41e-01 (N=16) | 5.42e-13 (N=96) | 1.0 |
| chebyshev | spectral | - | - | 4.16e-01 (N=16) | 4.81e-13 (N=96) | 5.7 |
| shooting | 4 | 4.00 | 4.00 | 5.14e-06 (N=200) | 1.28e-09 (N=1600) | 546.1 |

### stationary_morse

| Method | Nominal order | Fitted order | Last-pair order | Error, smallest N | Error, largest N | Time, largest N (ms) |
|---|---|---|---|---|---|---|
| fd2 | 2 | 2.03 | 2.00 | 8.57e-01 (N=40) | 7.83e-04 (N=1280) | 2.3 |
| fd4 | 4 | 3.87 | 4.00 | 2.53e-01 (N=40) | 4.59e-07 (N=1280) | 390.6 |
| numerov | 4 | 3.94 | 4.00 | 1.15e-01 (N=40) | 1.72e-07 (N=1280) | 480.1 |
| fourier | spectral | - | - | 2.93e+00 (N=16) | 3.95e-10 (N=96) | 0.7 |
| chebyshev | spectral | - | - | 1.02e+00 (N=16) | 5.51e-14 (N=96) | 4.1 |
| shooting | 4 | 4.00 | 4.00 | 2.87e-04 (N=200) | 7.06e-08 (N=1600) | 528.7 |

### stationary_wavefunctions

| Method | max error of psi_0 | max error of psi_3 | max energy error (6 levels) |
|---|---|---|---|
| fd2 | 1.87e-04 | 1.82e-03 | 1.21e-02 |
| fd4 | 7.52e-07 | 1.32e-05 | 9.61e-05 |
| numerov | 2.84e-07 | 5.00e-06 | 3.63e-05 |
| fourier | 2.27e-14 | 5.49e-12 | 5.11e-14 |
| chebyshev | 9.44e-15 | 5.45e-12 | 5.60e-14 |
| shooting | 2.84e-07 | 5.00e-06 | 3.63e-05 |

### double_well_splitting

| Method | N | E_0 | E_1 - E_0 (tunnel splitting) |
|---|---|---|---|
| fd2 | 400 | 1.33836049 | 9.717665e-03 |
| fd4 | 400 | 1.33855177 | 9.716344e-03 |
| numerov | 400 | 1.33855187 | 9.716348e-03 |
| shooting | 400 | 1.33855187 | 9.716348e-03 |
| fourier | 96 | 1.33855193 | 9.716350e-03 |
| chebyshev | 96 | 1.33855193 | 9.716350e-03 |

### timedep_orders

| Method | Nominal dt-order | 1D fitted | 1D last pair | 2D fitted | 2D last pair | h-order fitted (1D) | h-order last pair |
|---|---|---|---|---|---|---|---|
| crank-nicolson | 2 | 1.96 | 2.00 | 1.89 | 2.00 | 2.00 | 2.00 |
| split-step | 2 | 2.00 | 2.00 | 2.00 | 2.00 | spectral | - |
| rk4 | 4 | 4.00 | 4.00 | 4.00 | 4.00 | 2.00 | 2.00 |

### conservation_1d

| Method | dt | max abs(norm - 1) | max rel. energy drift | L2 error vs analytic | Time (s) |
|---|---|---|---|---|---|
| crank-nicolson | 2.0e-03 | 1.33e-14 | 2.74e-14 | 1.43e-02 | 0.178 |
| split-step | 2.0e-03 | 4.06e-13 | 6.43e-07 | 3.02e-06 | 0.195 |
| rk4 | 2.0e-03 | 1.91e-12 | 1.08e-11 | 1.42e-02 | 0.284 |

### conservation_2d

| Method | dt | max abs(norm - 1) | max rel. energy drift | L2 error vs analytic | Time (s) |
|---|---|---|---|---|---|
| crank-nicolson | 5.0e-03 | 5.33e-15 | 1.01e-14 | 1.05e-01 | 1.825 |
| split-step | 5.0e-03 | 2.84e-13 | 4.23e-06 | 2.75e-05 | 0.206 |
| rk4 | 5.0e-03 | 1.90e-09 | 8.09e-09 | 1.03e-01 | 0.609 |

### rk4_stability

| dt | dt x rho(H) | norm after t = 2 |
|---|---|---|
| 2.92e-03 | 2.00 | 1.000e+00 |
| 5.12e-03 | 3.50 | inf |

### double_well_dynamics

| Quantity | Value |
|---|---|
| splitting E1-E0 from stationary Fourier solver | 9.716350e-03 |
| splitting from time-dependent first minimum (pi / t_min) | 9.716381e-03 |
| relative difference | 3.25e-06 |
| max abs(norm - 1) over the run | 3.07e-13 |

### scattering_2d

| Quantity | Value |
|---|---|
| grid | 128 x 128, L = 24.0 |
| dt, steps | 0.01, 400 |
| transmitted probability (split-step) | 0.60777 |
| transmitted probability (Crank-Nicolson) | 0.58906 |
| relative L2 difference of the two wave functions | 1.18e-01 |
| same difference with the grid spacing doubled (half as many points) | 4.41e-01 |
| max abs(norm - 1): split-step / CN | 2.71e-14 / 2.22e-16 |
| wall time: split-step / CN (s) | 0.16 / 1.72 |

### tunneling_fd

| Points in barrier | h | T (FD) | T (exact) | abs error | order | abs(R + T - 1) |
|---|---|---|---|---|---|---|
| 5 | 3.10e-01 | 0.277414 | 0.280540 | 3.13e-03 | - | 4.4e-16 |
| 11 | 1.41e-01 | 0.279917 | 0.280540 | 6.23e-04 | 2.05 | 1.1e-15 |
| 21 | 7.38e-02 | 0.280371 | 0.280540 | 1.70e-04 | 2.01 | 1.3e-15 |
| 41 | 3.78e-02 | 0.280496 | 0.280540 | 4.45e-05 | 2.00 | 2.2e-15 |
| 81 | 1.91e-02 | 0.280529 | 0.280540 | 1.14e-05 | 2.00 | 8.9e-16 |
| 161 | 9.63e-03 | 0.280537 | 0.280540 | 2.88e-06 | 2.00 | 1.3e-15 |

### tunneling_packet

| k0 | mean energy k0^2/2 | T(mean energy) | T averaged over packet (exact) | T from split-step | abs difference |
|---|---|---|---|---|---|
| 0.60 | 0.180 | 0.04426 | 0.05107 | 0.05094 | 1.3e-04 |
| 0.80 | 0.320 | 0.09004 | 0.10021 | 0.10011 | 1.1e-04 |
| 1.00 | 0.500 | 0.16500 | 0.17887 | 0.17879 | 8.2e-05 |
| 1.20 | 0.720 | 0.28054 | 0.29560 | 0.29548 | 1.2e-04 |
| 1.60 | 1.280 | 0.63041 | 0.62534 | 0.62496 | 3.7e-04 |
| 2.00 | 2.000 | 0.92366 | 0.90329 | 0.90307 | 2.2e-04 |

### poisson_iterations

| Solver | n = 31 | n = 63 | n = 127 | n = 255 |
|---|---|---|---|---|
| direct | 0.003 s (direct) | 0.010 s (direct) | 0.046 s (direct) | 0.310 s (direct) |
| cg | 106 it / 0.003 s | 215 it / 0.010 s | 434 it / 0.059 s | 874 it / 0.530 s |
| cg-jacobi | 106 it / 0.003 s | 215 it / 0.011 s | 434 it / 0.062 s | 874 it / 0.528 s |
| cg-ic0 | 35 it / 0.004 s | 67 it / 0.017 s | 132 it / 0.105 s | 263 it / 0.780 s |
| mg | 13 it / 0.009 s | 13 it / 0.012 s | 14 it / 0.024 s | 14 it / 0.075 s |
| mg-cg | 8 it / 0.006 s | 8 it / 0.015 s | 8 it / 0.016 s | 8 it / 0.055 s |

### poisson_discretisation

| n | h | L2 error of the 5-point solution | observed order |
|---|---|---|---|
| 31 | 3.12e-02 | 4.75e-05 | - |
| 63 | 1.56e-02 | 1.19e-05 | 2.00 |
| 127 | 7.81e-03 | 2.97e-06 | 2.00 |
| 255 | 3.91e-03 | 7.42e-07 | 2.00 |

### heat_orders

| Scheme | dt = 0.0100 | dt = 0.0050 | dt = 0.0025 | dt = 0.0013 | last-pair order |
|---|---|---|---|---|---|
| backward-euler | 5.49e-03 | 2.49e-03 | 1.18e-03 | 5.71e-04 | 1.05 |
| crank-nicolson | 3.57e-04 | 9.02e-05 | 2.26e-05 | 5.66e-06 | 2.00 |

### heat_solvers

| Linear solver | inner iterations (20 steps) | wall time (s) | max difference to direct |
|---|---|---|---|
| direct | 1 solve per step | 0.201 | 0.00e+00 |
| cg | 973 | 0.053 | 1.89e-13 |
| cg-ic0 | 338 | 0.167 | 1.29e-13 |
| mg-cg | 160 | 0.181 | 5.57e-14 |

### environment

Python 3.13.17, NumPy 2.5.3, SciPy 1.18.1, 2 CPU cores, Linux-6.18.44-fc-v165-x86_64-with-glibc2.39
