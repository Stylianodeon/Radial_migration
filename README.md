# Radial Migration with Spiral Arm Perturbation

This project simulates stellar radial migration in a galactic disc with an axisymmetric logarithmic potential and rotating spiral-arm perturbations. It supports 2D and 3D configurations.

The orbit integrators evolve positions and velocities in the **inertial frame**. The spiral's rotation is represented by the explicit time dependence of its potential. Rotating-frame positions and effective potentials are used for diagnostics.

## Features

- Logarithmic background potential and spiral perturbations with constant pattern speed.
- Leapfrog integration with midpoint force evaluation.
- Energy, angular-momentum, radial-action, and Jacobi-energy diagnostics.
- Separate 2D initial-condition generator and simulation executables.

## Requirements for the 2D build

- A C++17 compiler.
- CMake 3.16 or later.

## Build the 2D executables

From the repository root:

```bash
cd 2D
cmake -S . -B build
cmake --build build
```

This builds two executables in `2D/build/bin/`:

| Executable | Purpose |
|---|---|
| `ShuDFExec` | Generate 30,000 stellar initial conditions using the Shu distribution function. |
| `LogPotExec` | Integrate the orbits from the existing initial-condition file. |

To build only one target, run the corresponding command from `2D/`:

```bash
cmake --build build --target ShuDFExec
cmake --build build --target LogPotExec
```

## Generate initial conditions

From `2D/`:

```bash
./build/bin/ShuDFExec
```

This writes `DF_initial_conditions.dat` in the current working directory, using a fixed random seed of 42. The six columns are:

```text
R  L  v_R  v_phi  x  y
```

Positions are in kpc, velocities in kpc/Myr, and specific angular momentum in kpc²/Myr. The default 2–20 kpc bounds apply to guiding radius.

The disc parameters are defined in `2D/Log_pot/include/Dehnen_DF.h`.

## Run a 2D simulation

From `2D/`:

1) For the axisymmetric logarithmic potential only:

```bash
./build/bin/LogPotExec -no_pert
```

For the logarithmic potential plus the spiral perturbation:

```bash
./build/bin/LogPotExec -pert
```

 The simulation writes `deltaLz_from_DF.dat`.


## Plot 3D results with gnuplot

With an existing 3D initial-condition file, you can build only the simulation target and save separate control and perturbed results. From the repository root:

```bash
cd 3D
cmake -S . -B build
cmake --build build --target LogPotExec
./build/bin/LogPotExec -no_pert -out control_3D.dat
./build/bin/LogPotExec -pert -out perturbed_3D.dat
```

Without `-out`, the default population output is `deltaLz_from_DF_3D.dat`. The examples below use `perturbed_3D.dat`; replace that filename if needed. The present command-line workflow writes population summaries, not individual orbit trajectories.

### Initial spatial distribution

The initial-condition file has these columns:

```text
1:id  2:R  3:phi  4:z  5:v_R  6:v_phi  7:v_z  8:Lz  9:Rg
10:x  11:y  12:z_cart  13:vx  14:vy  15:vz_cart
```
