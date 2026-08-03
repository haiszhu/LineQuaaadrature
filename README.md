# LineQuaaadrature (outdated)

Line-quadrature reference code (Fortran + MATLAB MEX on Mac).

## Context
This code was written by Codex and Claude Code as a testbed to reproduce parts of:
- Shidong Jiang and Hai Zhu, *Recursive reduction quadrature for the evaluation of Laplace layer potentials in three dimensions*, arXiv:2411.08342 (2024).  
  https://arxiv.org/abs/2411.08342
- John S. Asvestas, *Line integrals and physical optics. Part I. The transformation of the solid-angle surface integral to a line integral*, Journal of the Optical Society of America A, 2(6): 891-895 (1985).  
  https://doi.org/10.1364/JOSAA.2.000891

One purpose is to provide extended-precision solid-angle integral tests, to debug whether machine-precision singular / near-singular quadrature is achievable.

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
./build/test_solid_angle adaptive
./build/test_solid_angle nearroot
```

```bash
make test
./build/test_solid_angle false
./build/test_solid_angle true
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

### Extended precision, nearroot (`./build/test_solid_angle_r128 nearroot`)
Parameters: `order=16`, `nq=16`, `nvr=768`, `ntri=768`, `ratio=1.0`

```text
=== test_solid_angle_r128 (quad precision) ===
order=16  nq=16  nvr=768  ntri=768  ratio= 1.0
use_nearroot = T
Mesh built in    4.31 s

--- Part A: solid angle sum ---
far exterior  (1.5,0,0): |omega_sum| =   2.789740E-32
near exterior (1.001,0,0): |omega_sum| =   1.224617E-31
near exterior (1.000001,0,0): |omega_sum| =   4.902802E-30
near exterior (1.000000001,0,0): |omega_sum| =   3.646596E-27

--- Part B: DLP accuracy (direct sum + near correction) ---
N_src = 589824
ntarget = 834 exterior points
Direct DLP sum ...
  done in   50.28 s
Near-field correction ...
  done in   16.45 s

--- Result: DLP(1) = 0 outside ---
max |u|  =   2.173480E-32
rms |u|  =   3.423790E-33
argmax: i=630  tx=  5.555555555555555E-02  2.000000000000000E+00  6.944444444444444E-01
u(imax)        =  -2.173480294783766E-32
K_corr(imax)   =  -4.615084860876367E-33
```

### Double precision, nearroot (`./build/test_solid_angle nearroot`)
Parameters: `order=14`, `mp=np=8`, `ntri=768`, `ratio=1.0`

```text
=== test_solid_angle ===
order=14  mp=np=8  ntri=768  ratio= 1.0
use_nearroot = T
Mesh built in    0.01 s

--- Part A: solid angle sum ---
far exterior  (1.5,0,0): |omega_sum| =  8.428E-14
near exterior (1.1,0,0): |omega_sum| =  2.639E-13

--- Part B: DLP accuracy (direct sum + near correction) ---
N_src = 80640
ntarget = 6598 exterior points
Direct DLP sum ...
  done in    0.17 s
Near-field correction ...
  done in    0.01 s

--- Result: DLP(1) = 0 outside ---
max |u|  =  3.419E-14
rms |u|  =  1.260E-15
```
## To do list

* solid angle non-adaptive