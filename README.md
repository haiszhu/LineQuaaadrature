# LineQuaaadrature

Line-quadrature reference code (Fortran + MATLAB MEX on Mac).

## Context
This code was written by Codex and Claude Code as a testbed to reproduce parts of:
- Shidong Jiang and Hai Zhu, *Recursive reduction quadrature for the evaluation of Laplace layer potentials in three dimensions*, arXiv:2411.08342 (2024).  
  https://arxiv.org/abs/2411.08342
- John S. Asvestas, *Line integrals and physical optics. Part I. The transformation of the solid-angle surface integral to a line integral*, Journal of the Optical Society of America A, 2(6): 891-895 (1985).  
  https://doi.org/10.1364/JOSAA.2.000891

The main purpose is to provide extended-precision solid-angle integral tests, to debug whether machine-precision singular / near-singular quadrature is achievable.

## Layout
- `src/` Fortran source modules
- `matlab/` mwrap interface (`LineQuaaadrature.mw`) and generated MEX wrappers
- `test/` Fortran test programs
- `build/` compiled artifacts

## Build
From this folder:

```bash
cd /LineQuaaadrature
make
```

This builds:
- static library: `build/libLineQuaaadrature.a`
- MATLAB MEX: `matlab/LineQuaaadrature_mex.mexmaca64` (on Apple Silicon macOS)

## Run Fortran Tests
### Real64 test
```bash
make test
./build/test_solid_angle
```

### Real128 test (extended precision)
```bash
make test_r128
./build/test_solid_angle_r128
```

Optional thread control:
```bash
OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 ./build/test_solid_angle_r128
```

## Clean
```bash
make clean
```

## Notes
- Compiler settings are in `Makefile` (currently `gfortran-15`, `gcc-15`).
- MEX wrappers are generated from `matlab/LineQuaaadrature.mw` using `mwrap`.

## Reference Results
The following reference outputs were recorded on local runs.

### Extended precision (`./build/test_solid_angle_r128`)
Parameters: `order=14`, `nq=12`, `nvr=432`, `ntri=768`, `ratio=1.0`

```text
=== test_solid_angle_r128 (quad precision) ===
order=14  nq=12  nvr=432  ntri=768  ratio= 1.0
Mesh built in    2.32 s

--- Part A: solid angle sum ---
far exterior  (1.5,0,0): |omega_sum| =   7.080390E-31
near exterior (1.1,0,0): |omega_sum| =   9.781399E-31

--- Part B: DLP accuracy (direct sum + near correction) ---
N_src = 331776
ntarget = 834 exterior points
Direct DLP sum ...
  done in   23.57 s
Near-field correction ...
  done in   33.85 s

--- Result: DLP(1) = 0 outside ---
max |u|  =   2.751315E-28
rms |u|  =   4.642503E-29
```

### Double precision (`./build/test_solid_angle`)
Parameters: `order=14`, `mp=np=8`, `ntri=768`, `ratio=1.0`

```text
=== test_solid_angle ===
order=14  mp=np=8  ntri=768  ratio= 1.0
Mesh built in    0.02 s

--- Part A: solid angle sum ---
far exterior  (1.5,0,0): |omega_sum| =  4.846E-14
near exterior (1.1,0,0): |omega_sum| =  1.312E-14

--- Part B: DLP accuracy (direct sum + near correction) ---
N_src = 80640
ntarget = 6598 exterior points
Direct DLP sum ...
  done in    0.13 s
Near-field correction ...
  done in    0.02 s

--- Result: DLP(1) = 0 outside ---
max |u|  =  3.508E-14
rms |u|  =  1.267E-15
```
