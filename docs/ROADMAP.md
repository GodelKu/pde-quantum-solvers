# Roadmap

Planned extensions. Each one has a test in `tests/test_roadmap.py` marked `xfail(strict=True)`:
it fails today, and the suite turns red as soon as an implementation makes it pass while the
marker is still there. Implement, remove the marker, keep CI green. The known answers below are
exact results or were cross-checked by a second method; derive each one before coding it.

1. **Richardson-extrapolated eigenvalues.** Add
   `stationary.richardson_energies(method, potential, x_min, x_max, n_points, n_states)`. Solve
   on `n` and `2n + 1` interior points (the spacing halves exactly) and combine
   `(2^p E_fine − E_coarse) / (2^p − 1)`, with `p` the nominal order of the method. For these
   symmetric stencils the eigenvalue error is a series in even powers of `h`, so the leading term
   cancels: `fd2` becomes fourth order and `numerov` sixth order.
   Known answer: on the harmonic oscillator, `fd2` with 159 and 319 points is below 1e-5 (plain
   `fd2` at 319 points: about 2e-3), and the observed orders are 4 for `fd2` and 6 for `numerov`.

2. **Ground state by imaginary-time propagation.** Add
   `timedep.imaginary_time_ground_state(grid, potential, psi0, dtau, n_steps)` returning
   `(energy, psi)`. Replace `t` by `−iτ` in the split-step propagator,
   `e^{−V dτ/2} e^{−T dτ} e^{−V dτ/2}`, and renormalise after every step; the excited components
   decay like `e^{−(E_n − E_0) τ}`. Report the energy as `<ψ|H|ψ>` with the Fourier Hamiltonian,
   not as the decay rate of the norm, and explain why the first is much more accurate (the fixed
   point is the ground state of a modified Hamiltonian `H + O(dτ²)`).
   Known answer: at `dτ = 0.01` after `τ = 20`, the oscillator gives `E_0 = 1/2` within 1e-9 and
   the Morse potential `E_0 = ω₀/2 − ω₀²/(16 D)` within 1e-6.

3. **A 2D stationary solver and the degeneracies of the oscillator.** Add
   `stationary2d.eigenvalues_2d(potential, length, n, n_states)` on the square
   `[−L/2, L/2]²` with Dirichlet walls. The five-point Hamiltonian is the Kronecker sum
   `(I ⊗ T + T ⊗ I)/2 + V` of 1D second-difference matrices; use
   `scipy.sparse.linalg.eigsh` in shift-invert mode.
   Known answer: the isotropic oscillator has levels 1, 2, 2, 3, 3, 3, 4, 4, 4, 4 (within 0.02 at
   `n = 127`, `L = 16`). The pair (1, 0), (0, 1) is degenerate to round-off by the `x ↔ y`
   symmetry, while (1, 1) differs from (2, 0) and (0, 2) at `O(h²)`; explain why. With
   `ω_y = 2` the levels are `(n_x + 1/2) + 2 (n_y + 1/2)`.

4. **Transfer matrices and resonant tunnelling.** Add
   `tunneling.transfer_matrix_transmission(energies, edges, heights)` for a piecewise-constant
   potential with `V = 0` outside (`heights` has one more entry than `edges`). Match `ψ` and `ψ'`
   at each edge with 2 × 2 matrices; inside a barrier below its top the wavenumber is imaginary.
   Known answers: a single barrier reproduces the closed-form `T(E)` of `potentials.py` to 1e-12.
   Two barriers of height 1 and width 1 around a well of width 3 show two resonances below the
   barrier top, near `E = 0.2350` and `0.9542`, where `T = 1` although one of the barriers alone
   transmits only about 0.22 and 0.65 there. These energies agree to about 1e-6 with a
   finite-difference scattering calculation at `h = 0.005` (the lattice method of
   `transmission_fd`, extended to two barriers).

5. **An L-stable start for Crank–Nicolson.** Add `startup_steps` to `elliptic.heat_solve`: the
   first `startup_steps` Crank–Nicolson steps are each replaced by two backward-Euler steps of
   `dt/2` (the start proposed by Rannacher, Numerische Mathematik 43, 1984). Backward Euler damps
   the stiff modes that Crank–Nicolson only flips; `I/(dt/2) − Δ_h` is the matrix of both the
   half step and the Crank–Nicolson step, so one factorisation serves both.
   Known answer: on the discontinuous hot patch of the README, one step with `dt = 0.001` stays
   non-negative to round-off (plain Crank–Nicolson undershoots to about −0.11), and on smooth data
   the method stays second order.
