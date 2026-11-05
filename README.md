# Radial Migration with Spiral Arm Perturbation

This project simulates stellar radial migration in a galactic disc with a background axisymmetric logarithmic potential and spiral arm peerturbations.The orbit integrator evolves positions and velocities in the **inertial frame**. Radial migration and kinematic heeating are recorded through changes of the angular momentum and radial action of stars respectively.

## Requirements for the 2D build

- A C++17 compiler.
- CMake 3.16 or later.

## Build the 2D executables

From the repository root:

```bash
cd 2D
cmake -B build
cmake --build build
```

This builds two executables in `2D/build/bin/`:

| Executable | Purpose |
|---|---|
| `ShuDFExec` | Generate 30,000 stellar initial conditions using the Shu distribution function. |
| `LogPotExec` | Integrate the orbits from the existing initial-condition file. |

## Generate initial conditions

From `2D/`:

```bash
./build/bin/ShuDFExec
```

This creates a data file `DF_initial_conditions.dat` with the initial condition parameters of 30,000 stars in the current working directory. The recorded parameters are:

```text
R  L  v_R  v_phi  x  y
```

Positions are in kpc, velocities in kpc/Myr, and angular momentum in kpc²/Myr.

## Run a 2D simulation

From `2D/`:

1) For the background axisymmetric logarithmic potential only:

```bash
./build/bin/LogPotExec -no_pert
```

For the logarithmic potential plus the spiral perturbation:

```bash
./build/bin/LogPotExec -pert
```

 The simulation results are stored in `sim_results.dat`.

## Effective potential and Lagrange-points

We calculate the effective potential as

$$
\Phi_{\mathrm{eff}} = \Phi_0 + \Phi_{\mathrm{spiral}} - \frac{1}{2}\Omega_p^2 R^2,
$$

where $\Phi_0$ is the logarithmic background potential. At each snapshot,
Lagrange points satisfy $\nabla\Phi_{\mathrm{eff}}=0$: we solve the radial equation by bisection and calculate the corresponding angles analytically.

To reproduce the animation, run from `2D/`:

```bash
cmake -B build
cmake --build build --target LagrangeSnapshotsExec
./build/bin/LagrangeSnapshotsExec build/lagrange_snapshots
python3 plots/animate_lagrange.py build/lagrange_snapshots
```
The animation is reproduced from snapshots taken every 10Myr.

Python libraries required: Matplotlib, NumPy, and Pillow.


The GIF is saved as `2D/build/lagrange_snapshots/lagrange_motion.gif`.
It shows the instantaneous equilibria in inertial and rotating coordinates
from 1900 to 3100 Myr.

## G 3D results with gnuplot

For the 3D disk the same commands should be used in the terminal form the 3D directory.
