# Classical Diatomic Scattering Simulator

A Python simulation for computing classical mechanical scattering angles and trajectories for diatomic molecules interacting via a Morse potential.

## What It Does

The program integrates Newton's equations of motion for two atoms and computes:

1. **Scattering angle θ vs. impact parameter b** — shows how deflection depends on the collision geometry
2. **2D Trajectories** — visualizes the paths of all simulated collisions
3. **Total energy over time** — verifies energy conservation throughout each trajectory

## How It Works

1. **Input**: The user provides molecular parameters (masses, vibrational frequency, anharmonicity, equilibrium distance, dissociation energy) and collision conditions (initial displacement and kinetic energy).
2. **Morse Potential**: Atom–atom interaction is modeled with a Morse potential, with the well depth `De` calculated from the dissociation energy plus the zero-point energy.
3. **Trajectory Integration**: Newton's equations are integrated using a simple finite-difference (Euler) scheme over a fixed number of time steps for each impact parameter `b`.
4. **Scattering Angle**: After finding the turning point from each trajectory, the scattering integral is evaluated numerically using the trapezoidal rule until convergence.

## Dependencies

```
numpy
scipy
matplotlib
```

Install with:
```bash
pip install numpy scipy matplotlib
```

## Usage

Run the script and follow the prompts:

```bash
python scat_ang_CM_for_diatomics.py
```

You will be asked to enter:
- `init_d` — initial atomic separation in Å
- `Ekin_initial` — initial kinetic energy in eV

> **Note:** Molecular constants (`m1`, `m2`, `ve`, `vexe`, `re`, `D0`) are currently hardcoded for F₂ but can be changed directly in the script or uncommented to accept user input.

## Output

Three `matplotlib` plots are displayed after the simulation completes:
- `θ/180°` vs. `b` with the rainbow angle marked
- Collision trajectories in the x–y plane
- Total energy per time step for each trajectory (energy conservation check)

