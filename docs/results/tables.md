### stationary_harmonic

| Method | Nominal order | Fitted order | Last-pair order | Error, smallest N | Error, largest N | Time at largest N (ms) |
|---|---|---|---|---|---|---|
| fd2 | 2 | 2.01 | 2.00 | 1.24e-01 (N=40) | 1.22e-04 (N=1280) | 2.6 (2.4-2.7) |
| fd4 | 4 | 3.98 | 4.00 | 1.41e-02 (N=40) | 1.60e-08 (N=1280) | 409.7 (400.7-438.1) |
| numerov | 4 | 4.01 | 4.00 | 5.98e-03 (N=40) | 5.99e-09 (N=1280) | 485.9 (455.7-501.8) |
| fourier | spectral | - | - | 1.00e-01 (N=16) | 2.66e-14 (N=96) | 0.7 (0.6-0.8) |
| chebyshev | spectral | - | - | 6.47e-01 (N=16) | 2.49e-14 (N=96) | 3.3 (3.2-3.5) |
| shooting | 4 | 4.00 | 4.00 | 9.90e-06 (N=200) | 2.45e-09 (N=1600) | 553.3 (433.4-614.7) |

### stationary_double_well

| Method | Nominal order | Fitted order | Last-pair order | Error, smallest N | Error, largest N | Time at largest N (ms) |
|---|---|---|---|---|---|---|
| fd2 | 2 | 2.01 | 2.00 | 8.07e-02 (N=40) | 8.06e-05 (N=1280) | 2.4 (2.3-2.5) |
| fd4 | 4 | 3.98 | 4.00 | 7.32e-03 (N=40) | 8.30e-09 (N=1280) | 381.9 (379.3-433.9) |
| numerov | 4 | 4.01 | 4.00 | 3.10e-03 (N=40) | 3.11e-09 (N=1280) | 489.9 (470.3-565.3) |
| fourier | spectral | - | - | 1.41e-01 (N=16) | 5.42e-13 (N=96) | 0.8 (0.7-2.9) |
| chebyshev | spectral | - | - | 4.16e-01 (N=16) | 4.81e-13 (N=96) | 3.5 (3.3-3.8) |
| shooting | 4 | 4.00 | 4.00 | 5.14e-06 (N=200) | 1.28e-09 (N=1600) | 499.8 (426.5-581.6) |

### stationary_morse

| Method | Nominal order | Fitted order | Last-pair order | Error, smallest N | Error, largest N | Time at largest N (ms) |
|---|---|---|---|---|---|---|
| fd2 | 2 | 2.03 | 2.00 | 8.57e-01 (N=40) | 7.83e-04 (N=1280) | 2.8 (2.5-3.2) |
| fd4 | 4 | 3.87 | 4.00 | 2.53e-01 (N=40) | 4.59e-07 (N=1280) | 399.5 (376.6-470.8) |
| numerov | 4 | 3.94 | 4.00 | 1.15e-01 (N=40) | 1.72e-07 (N=1280) | 504.3 (458.1-526.5) |
| fourier | spectral | - | - | 2.93e+00 (N=16) | 3.95e-10 (N=96) | 0.6 (0.6-0.7) |
| chebyshev | spectral | - | - | 1.02e+00 (N=16) | 5.51e-14 (N=96) | 3.1 (3.0-3.2) |
| shooting | 4 | 4.00 | 4.00 | 2.87e-04 (N=200) | 7.06e-08 (N=1600) | 457.9 (428.7-547.0) |

### hard_walls

| Problem | Method | Error, smallest N | Error, largest N | Last-pair order |
|---|---|---|---|---|
| infinite-well | fd2 | 6.47e-01 (N=39) | 1.01e-02 (N=319) | 2.00 |
| infinite-well | fd4 | 8.47e-03 (N=39) | 2.09e-06 (N=319) | 4.00 |
| infinite-well | numerov | 3.22e-03 (N=39) | 7.82e-07 (N=319) | 4.00 |
| infinite-well | shooting | 3.22e-03 (N=39) | 7.82e-07 (N=319) | 4.00 |
| infinite-well | chebyshev | 3.98e-01 (N=8) | 2.70e-13 (N=24) | spectral |
| linear-ramp | fd2 | 3.94e-02 (N=149) | 6.05e-04 (N=1199) | 2.00 |
| linear-ramp | fd4 | 1.53e-03 (N=149) | 3.86e-07 (N=1199) | 4.00 |
| linear-ramp | numerov | 5.83e-04 (N=149) | 1.41e-07 (N=1199) | 4.00 |
| linear-ramp | shooting | 5.83e-04 (N=149) | 1.41e-07 (N=1199) | 4.00 |
| linear-ramp | chebyshev | 3.53e-04 (N=32) | 1.92e-13 (N=96) | spectral |

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

### morse_box

| Box | Chebyshev error, N = 96 | Chebyshev error, N = 192 |
|---|---|---|
| [-2, 18] | 5.51e-14 | 1.42e-14 |
| [-1.5, 18] | 1.55e-07 | 1.55e-07 |

### fd4_cost

| N | fd2, tridiagonal driver, 4 eigenpairs (ms) | fd4, banded driver, 4 eigenvalues (ms) | fd4, banded driver, 4 eigenpairs (ms) |
|---|---|---|---|
| 320 | 0.6 (0.5-0.6) | 0.8 (0.8-0.8) | 2.8 (2.8-2.9) |
| 640 | 1.1 (1.1-1.1) | 2.2 (2.2-2.2) | 35.1 (30.7-35.2) |
| 1280 | 2.3 (2.2-2.3) | 8.9 (8.7-9.0) | 396.0 (363.3-409.3) |

### timedep_orders

| Method | Nominal dt-order | 1D fitted | 1D last pair | 2D fitted | 2D last pair | h-order fitted (1D) | h-order last pair |
|---|---|---|---|---|---|---|---|
| crank-nicolson | 2 | 1.96 | 2.00 | 1.89 | 2.00 | 2.00 | 2.00 |
| split-step | 2 | 2.00 | 2.00 | 2.00 | 2.00 | spectral | - |
| rk4 | 4 | 4.00 | 4.00 | 4.00 | 4.00 | 2.00 | 2.00 |

### conservation_1d

| Method | dt | max abs(norm(t) - norm(0)) | max rel. energy drift | L2 error vs analytic | Time (s) |
|---|---|---|---|---|---|
| crank-nicolson | 2.0e-03 | 1.33e-14 | 2.74e-14 | 1.43e-02 | 0.151 (0.124-0.180) |
| split-step | 2.0e-03 | 4.06e-13 | 6.43e-07 | 3.02e-06 | 0.154 (0.131-0.171) |
| rk4 | 2.0e-03 | 1.91e-12 | 1.08e-11 | 1.42e-02 | 0.291 (0.262-0.371) |

### conservation_2d

| Method | dt | max abs(norm(t) - norm(0)) | max rel. energy drift | L2 error vs analytic | Time (s) |
|---|---|---|---|---|---|
| crank-nicolson | 5.0e-03 | 5.33e-15 | 1.01e-14 | 1.05e-01 | 2.181 (1.926-2.300) |
| split-step | 5.0e-03 | 2.84e-13 | 4.23e-06 | 2.75e-05 | 0.260 (0.226-0.316) |
| rk4 | 5.0e-03 | 1.90e-09 | 8.09e-09 | 1.03e-01 | 0.604 (0.591-0.629) |

### rk4_stability

| dt | dt x rho(H) | growth per step of the top mode, abs(R(-i dt rho)) | norm after t = 20 |
|---|---|---|---|
| 2.921e-03 | 2.000 | 0.745 | 1.000e+00 |
| 4.089e-03 | 2.800 | 0.931 | 1.000e+00 |
| 4.176e-03 | 2.860 | 1.081 | 3.939e+121 |
| 5.111e-03 | 3.500 | 3.816 | inf |

### double_well_dynamics

Reference E1 - E0 (chebyshev, N = 200): 9.716349909e-03

| dt | E1 - E0 from the fit of P(x < 0, t) | relative difference | order | E1 - E0 of the Strang step (eigenphases) | relative difference | max abs(norm(t) - norm(0)) |
|---|---|---|---|---|---|---|
| 0.04 | 9.724052849e-03 | 7.93e-04 | - | 9.724053026e-03 | 7.93e-04 | 6.01e-13 |
| 0.02 | 9.718268974e-03 | 1.98e-04 | 2.01 | 9.718269126e-03 | 1.98e-04 | 6.86e-13 |
| 0.01 | 9.716829157e-03 | 4.93e-05 | 2.00 | 9.716829304e-03 | 4.93e-05 | 1.40e-13 |

### scattering_2d

| Quantity | Value |
|---|---|
| grid | 128 x 128, L = 24.0 |
| dt, steps | 0.01, 400 |
| transmitted probability (split-step) | 0.60777 |
| transmitted probability (Crank-Nicolson) | 0.58906 |
| relative L2 difference of the two wave functions | 1.18e-01 |
| same difference with the grid spacing doubled (half as many points) | 4.41e-01 |
| max abs(norm(t) - norm(0)): split-step / CN | 2.71e-14 / 2.22e-16 |
| wall time, split-step (s) | 0.23 (0.20-0.26) |
| wall time, Crank-Nicolson (s) | 1.72 (1.52-2.29) |

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

| k0 | mean energy k0^2/2 | T(mean energy) | T averaged over packet (exact) | T from split-step | split-step minus averaged exact |
|---|---|---|---|---|---|
| 0.60 | 0.180 | 0.04426 | 0.05107 | 0.05094 | -1.3e-04 |
| 0.80 | 0.320 | 0.09004 | 0.10021 | 0.10011 | -1.1e-04 |
| 1.00 | 0.500 | 0.16500 | 0.17887 | 0.17879 | -8.2e-05 |
| 1.20 | 0.720 | 0.28054 | 0.29560 | 0.29548 | -1.2e-04 |
| 1.60 | 1.280 | 0.63041 | 0.62534 | 0.62496 | -3.7e-04 |
| 2.00 | 2.000 | 0.92366 | 0.90329 | 0.90307 | -2.2e-04 |

### tunneling_packet_later

| k0 | box length | h | difference at 1 t1 | difference at 2 t1 | difference at 3 t1 |
|---|---|---|---|---|---|
| 0.60 | 800 | 0.05 | -1.6e-04 | +1.6e-04 | +1.7e-04 |
| 1.20 | 800 | 0.05 | +2.2e-05 | -4.6e-05 | -1.7e-05 |

### poisson_iterations

| Solver | n = 31 | n = 63 | n = 127 | n = 255 |
|---|---|---|---|---|
| cg | 106 | 215 | 434 | 874 |
| cg-jacobi | 106 | 215 | 434 | 874 |
| cg-ic0 | 35 | 67 | 132 | 263 |
| mg | 13 | 13 | 14 | 14 |
| mg-cg | 8 | 8 | 8 | 8 |

### poisson_times

| Solver | n = 31 (s) | n = 63 (s) | n = 127 (s) | n = 255 (s) |
|---|---|---|---|---|
| direct | 0.005 (0.005-0.006) | 0.015 (0.013-0.016) | 0.048 (0.047-0.048) | 0.316 (0.298-0.397) |
| cg | 0.004 (0.003-0.005) | 0.016 (0.015-0.016) | 0.065 (0.060-0.080) | 0.538 (0.517-0.591) |
| cg-jacobi | 0.004 (0.003-0.005) | 0.011 (0.011-0.013) | 0.078 (0.076-0.083) | 0.604 (0.565-0.624) |
| cg-ic0 | 0.006 (0.006-0.007) | 0.017 (0.016-0.017) | 0.105 (0.104-0.167) | 0.811 (0.778-0.952) |
| mg | 0.013 (0.010-0.013) | 0.012 (0.012-0.013) | 0.026 (0.024-0.027) | 0.079 (0.073-0.095) |
| mg-cg | 0.010 (0.009-0.010) | 0.009 (0.009-0.010) | 0.023 (0.019-0.025) | 0.063 (0.059-0.065) |

Fastest median in this run: n = 31: cg-jacobi, n = 63: mg-cg, n = 127: mg-cg, n = 255: mg-cg.

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

| Linear solver | iterations per step | one-time setup (s) | total, 20 steps incl. setup (s) | max difference to direct |
|---|---|---|---|---|
| direct | 1 LU solve | 0.010 (0.010-0.010) | 0.021 (0.020-0.023) | 0.00e+00 |
| cg | 48.6 | 0.000 (0.000-0.000) | 0.054 (0.050-0.063) | 1.89e-13 |
| cg-ic0 | 16.9 | 0.006 (0.005-0.008) | 0.084 (0.067-0.096) | 1.29e-13 |
| mg-cg | 8.0 | 0.001 (0.001-0.001) | 0.146 (0.144-0.156) | 5.57e-14 |

Fastest to slowest median in this run: direct, cg, cg-ic0, mg-cg.

### timing_note

Wall times are the median of 5 runs after one warm-up run, with the fastest and slowest of those runs in brackets.

### environment

Python 3.13.17, NumPy 2.5.3, SciPy 1.18.1, Linux
