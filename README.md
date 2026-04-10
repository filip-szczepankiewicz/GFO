# GFO - Geometric Filter Optimization for Rotation Sampling

A MATLAB toolbox for generating optimally isotropic sets of rotations in SO(3), with applications to axi-symmetric and triaxial diffusion MRI and other orientation-sensitive experiments.

Repository based on the paper:

[Beyond directions: Symmetry-aware rotation sets for triaxial diffusion encoding by geometric filter optimization](https://arxiv.org/abs/2601.14970), 

by Sune Nørhøj Jespersen and Filip Szczepankiewicz

---

## Overview

In diffusion MRI and related techniques, experiments are repeated across many orientations (directions or full rotations). The quality of subsequent analyses depends critically on how *uniformly* those orientations sample the relevant space. Poor sampling leads to orientation bias: the measured powder-average signal will vary depending on how the sample happens to be placed in the scanner.

**GFO** addresses this by optimizing a set of N rotations to be as orientation-independent (isotropic) as possible. It supports two families of cost functions:

- **GFO (Geometric Filter Optimization)** — a band-limited, spectral cost on SO(3), based on Chebyshev polynomials and Sobolev-type band weights. This directly minimizes the variance of the frame energy across orientations.
- **ESR (Electrostatic Repulsion)** — a Coulomb-like repulsion energy on SO(3) (or S²), which spreads rotations by maximizing pairwise geodesic distances.

Both families support a **D₂-symmetry** variant (suffix `D2`), appropriate for acquisition schemes where each rotation and its dihedral equivalents represent the same physical experiment.

---

## Repository structure

```
GFO/
├── GFO_generateSet.m          % High-level interface: generate a set of rotation matrices
├── isotropic_SO3.m            % Core SO(3) optimizer (returns quaternions)
├── isotropic_S2.m             % S² optimizer for directions (LTE/ESRS2 case)
├── cv_SO3.m                   % Evaluate isotropy via coefficient of variation
├── GFOSO3_cost_fct.m          % GFO cost function (base)
├── GFOSO3_cost_fct_D2.m       % GFO cost function with D2 symmetry
├── electroSO3_cost_fct.m      % Electrostatic repulsion cost on SO(3)
├── electroSO3_cost_fct_D2.m   % Electrostatic repulsion with D2 symmetry
├── example.m                  % Worked example with CV evaluation
├── +library/                  % Pre-computed and cached rotation sets
│   ├── loadSet.m              % Load a cached set by (n, mode)
│   ├── makeSet.m              % Compute and cache a set by (n, mode)
│   ├── GFOD2_NNNN             % Pre-computed GFO-D2 sets
│   ├── ESRD2_NNNN             % Pre-computed ESR-D2 sets
│   ├── ESRS2_NNNN             % Pre-computed ESR-S2 sets
│   └── StD_NNNN               % Pre-computed spherical t-design sets
└── +util/                     % Internal utility functions
│   ├── sobolev.m              % Sobolev spectral weights
│   ├── hopf_kronecker.m       % Deterministic S³ initializer
│   ├── canon_antipodes.m      % Antipodal canonicalization (q ↔ −q)
│   ├── normalize_rows.m       % Row-wise unit normalization
│   ├── euler_so3_grid.m       % Uniform Euler angle grid on SO(3)
│   ├── hopf_kronecker.m       % Hopf–Kronecker lattice on S³
│   ├── quat_to_euler.m        % Quaternion → ZYZ Euler angles
│   ├── quat_to_R.m            % Quaternion → rotation matrix
│   ├── quat_to_axis_angle.m   % Quaternion → axis–angle
│   ├── Rzyz.m                 % ZYZ Euler angles → rotation matrix
│   ├── randrotations.m        % Haar-random rotations
│   ├── rotVec2Vec.m           % Rotation matrix mapping one vector to another
│   ├── applyRotMatToVec.m     % Apply rotation matrices to vectors
│   ├── readRotMatsFromFile.m  % Read rotation matrices from file
│   └── writeRotMatsToFile.m   % Write rotation matrices to file
└── figures/                   % Scripts for reproducing figures
    ├── figure_gfo_vs_esr_for_LTE.m
    ├── figure_dti_fitting.m
    ├── figure_tensorDist.m
    ├── figure_triangleGlyphs.m
    ├── figure_whatIsTriaxial_v2.m
    ├── figure_colTria.m
    └── btensor_triangle.m
```

---

## Optimization modes

The following modes are available in this framework.

| Mode     | Description                                                                 |
|----------|-----------------------------------------------------------------------------|
| `GFOD2`  | GFO on SO(3) with D2 symmetry. **Recommended general-purpose choice.**      |
| `GFO`    | GFO on SO(3) without additional symmetry constraint.                        |
| `ESRD2`  | ESR on SO(3) with D2 symmetry.                                              |
| `ESR`    | ESR on SO(3) without additional symmetry constraint.                        |
| `ESRS2`  | ESR on S² (directions). Use for axi-symmetric encoding only.                |

GFO = Geometric Filter Optimization;
ESR = Electrostatic Repulsion

---

## Quick start

### Generate a rotation set directly

```matlab
% Generate 15 GFO-D2-optimized rotations
n = 15;
rotMats = GFO_generateSet(n, 'GFOD2');
% rotMats is a 3×3×n array of rotation matrices
```

### Load a pre-computed set from the library

```matlab
% Load a pre-computed GFOD2 set or rotation matrices with n = 15 rotations
R = library.loadSet(15, 'GFOD2');
```

### Optimize and write a set to the library

```matlab
library.makeSet(15, 'GFOD2');  % Skips if already exists
```

---

## Detailed example with CV evaluation

This follows `example.m` and illustrates the full workflow including evaluation.

```matlab
clear

%% Problem definition
nDirs = 15;

% Sobolev band weights for the GFO cost
s = 8; kappa = 7; Lmax = 8;
lvals = 0:Lmax;
Vdeg  = (2*lvals+1) .* util.sobolev(kappa, s, Lmax).^2;
Vdeg  = Vdeg .* (mod(lvals,2) == 0);  % use only even bands

% b-tensor (encoding) and diffusion tensor
B = diag([0 1 2]) / 2;       % e.g. a planar tensor encoding
D = diag([0.1 0.1 2.8]);     % prolate diffusion tensor

%% Optimize rotations
[t, x, y, z, fval, exitflag, output] = isotropic_SO3(nDirs, 'GFOD2', Lmax, Vdeg);

% Convert to other representations
Q         = [t, x, y, z];
[a, b, g] = util.quat_to_euler(t, x, y, z);
rotMats   = util.Rzyz(a, b, g);

%% Evaluate CV (isotropy quality)
Euler  = [a, b, g];
w      = ones(1, nDirs) / nDirs;  % uniform weights

% Use a quasi-uniform quadrature grid up to L = 6 for averaging
EU      = util.euler_so3_grid(6);
U       = util.Rzyz(EU(:,1), EU(:,2), EU(:,3));
cvOpts  = struct('batch', 5e3, 'useParfor', false, 'rots', U);

[CV, ~, sPowder, ~] = cv_SO3(Euler, w, B, D, size(U,3), cvOpts);

fprintf('CV = %.4f%%\n', CV * 100);
fprintf('Mean powder average = %.4f\n', mean(sPowder));
```

---

## Optimization details

### Initialization

Optimization starts from a deterministic **Hopf–Kronecker lattice** on S³ (`util.hopf_kronecker`), which provides a well-distributed, non-random starting point. This makes results reproducible without a fixed random seed.

### Optimizer

`isotropic_SO3` uses MATLAB's `fminunc` with the quasi-Newton algorithm and central finite differences. No analytic gradient is provided; the cost functions are smooth enough for reliable convergence with the default tolerances.

### Post-processing

After optimization, quaternions are row-normalized to enforce unit length and then canonicalized via `util.canon_antipodes` to enforce the SO(3) identification q ↔ −q.

### Hyperparameters

The GFO cost is controlled by three hyperparameters with sensible defaults:

| Parameter | Default SO(3) | Default S² | Description |
|-----------|---------------|------------|-------------|
| `Lmax`    | 8             | 8          | Maximum spherical harmonic band; higher = more precise spectral control |
| `s`       | 8             | 5          | Sobolev smoothness exponent; larger = faster high-frequency decay       |
| `kappa`   | 7             | 4          | Sobolev scale parameter; larger = more weight on higher bands           |

These defaults work well for the range of N relevant to typical diffusion MRI experiments (N ≈ 6–60). They can be tuned if the optimization does not converge or the CV is unexpectedly high.

---

## File format for pre-computed rotation matrices

Pre-computed sets in `+library/` are stored as plain-text files. Each row contains the 9 entries of a 3×3 rotation matrix in row-major order:

```
# Rotation matrices: N entries, row-major flattening of 3x3 -> 1x9
# https://github.com/Neurophysics-CFIN/GFO
# Col:  1    2    3    4    5    6    7    8    9
#       R11  R12  R13  R21  R22  R23  R31  R32  R33
   -0.63277280    0.51221459    0.58071921  ...
   ...
```

Files are named `MODE_NNNN`, e.g. `GFOD2_0015` for a GFOD2 set with N = 15. Read and write via `util.readRotMatsFromFile` / `util.writeRotMatsToFile`.

---

## Dependencies
The code was developed and tested on Matlab 2025b.

- **MATLAB** with the 'Optimization Toolbox' and 'Parallel Computing Toolbox'.
- No third-party dependencies
- Some figure scripts require the [`fix_matlab`](https://github.com/filip-szczepankiewicz/fix_matlab) repository

---

## Citation

If you use this toolbox in your work, please cite:

> *[citation to be added]*

The repository is hosted at: https://github.com/Neurophysics-CFIN/GFO
