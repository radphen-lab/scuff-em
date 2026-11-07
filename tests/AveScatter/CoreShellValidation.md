# Core-shell validation

This records the independent validation performed on 2026-11-07. The historical
temporary before/after comparison reproduces the pre-cleanup result, including
an absorption bug. Those historical files do not establish physical accuracy.

## Geometry and independent reference

The supplied case has a concentric core of radius 0.5 um and relative
permittivity `4+1i`, inside a shell of outer radius 1 um and relative
permittivity `10+0.5i`. The host is vacuum and all relative permeabilities are
one. The SCUFF frequencies are omega = 0.3 and 1, corresponding to vacuum
wavelengths 20.944 and 6.283 um. The outer size parameters are 0.3 and 1.

For a concentric sphere, rotational symmetry makes the orientation-averaged
cross sections equal to the plane-wave coated-sphere Mie solution.
`ValidateCoreShell.py` independently matches the electric and magnetic
Riccati-Bessel fields at both interfaces. It computes extinction, scattering,
and asymmetry from the multipole coefficients. It computes core absorption
from the inward Poynting flux at the inner interface and shell absorption
from the difference between the outer and inner inward fluxes.

The reference passes checks of multipole truncation, the homogeneous-sphere
limit, the quasistatic coated-sphere limit, and zero absorption in a lossless
core or shell. Its total absorption, scattering, and first angular moment
also agree with PyMieScatt's independent `CoreShell_ab` calculation to a
relative tolerance of 1e-11 at both frequencies. See the
[PyMieScatt coated-sphere documentation](https://pymiescatt.readthedocs.io/en/latest/forwardCS.html)
and its [implementation](https://github.com/bsumlin/PyMieScatt/blob/master/PyMieScatt/CoreShell.py).

| omega | Exact shell Cabs | Exact core Cabs | Exact total Cabs | Exact Csca | Exact g |
|---|---|---|---|---|---|
| 0.3 | 0.04783227158 | 0.02069559962 | 0.06852787120 | 0.03854184021 | 0.04214150277 |
| 1.0 | 3.269330286 | 0.7818301448 | 4.051160431 | 14.38040699 | 0.1743410934 |

Cross sections are in um^2. The asymmetry factor g is dimensionless. Column 6
of the BEM output is g times Csca. Compute g from the sum of column 6 divided
by the sum of column 5. For this geometry, row 01 contains shell absorption
and row 02 contains core absorption.

## Confirmed source error and correction

`MatrixTraceFlux` contracts an observable of shape D by S with a response of
shape S by D. It used the response's column count D for the S loop. When the
two interfaces have different numbers of basis functions, the off-diagonal
absorption contractions access entries outside the matrices. The correction
uses the observable's column count S.

The historical shell absorption was -121.518865 um^2 at omega = 0.3 and
345860.321 um^2 at omega = 1. The corrected values on the same original
meshes are 0.0464130304 and 3.45163369 um^2. Core absorption, scattering, and
the scattering moment on this case retain their previous values.

A bounds-checked test of the trace routine rejects the original code and
passes the corrected code for 4 by 2, 2 by 4, and 4 by 4 observables, with
and without the SCUFF matrix scaling. Fresh original-mesh runs with 1 and
4 OpenMP threads give identical printed results. Small-mesh oblate, dimer,
and lossy-sphere outputs also agree with their pre-cleanup outputs.

## Mesh convergence

The acceptance criterion is a relative error of at most 5% for each of shell
absorption, core absorption, total absorption, scattering, g times Csca,
and g, together with finite values, nonnegative material absorption, positive
scattering, and an asymmetry magnitude no greater than one. This is a chosen
BEM discretization tolerance. It is separate from the historical 1e-12
binary-comparison tolerance.

The original meshes contain 436 outer triangles and 218 inner triangles.
The medium meshes contain 820 and 450 triangles. The fine meshes contain
1258 and 708 triangles. The saved results are under `tests/AveScatter/reference/`.

Both fine-mesh frequencies pass the criterion. The largest relative errors
among the six observables decrease with mesh refinement.

| Mesh | Maximum error at omega = 0.3 | Maximum error at omega = 1 | Passes both frequencies |
|---|---|---|---|
| Original, after correction | 7.083% | 15.237% | No |
| Medium | 3.635% | 7.643% | No |
| Fine | 2.368% | 4.928% | Yes |

The corrected fine-mesh totals are

| omega | Total Cabs | Csca | g |
|---|---|---|---|
| 0.3 | 0.0676922502 | 0.0378628380 | 0.04188158347 |
| 1.0 | 4.108903509 | 14.37746210 | 0.1829321352 |

The remaining largest error is in g at omega = 1. At that frequency, total
absorption differs from Mie theory by 1.425% and scattering by 0.0205%.
The saved [fine output](reference/CoreShell.Fine.SIFlux.AVSCAT),
[medium output](reference/CoreShell.Medium.SIFlux.AVSCAT), and
[original-mesh output](reference/CoreShell.SIFlux.AVSCAT) contain the full
per-material and per-observable results.

## Reproducing the validation

The tests use the same `Test*.sh`, `reference/`, and checklist layout as the
Mie and Fresnel tests. From a configured build directory, run

```bash
make -C tests/AveScatter check
```

To run directly from the repository root, with the executable built and
Python 3, NumPy, and SciPy available, use

```bash
SCUFF_AVESCATTER="$PWD/applications/scuff-avescatter/scuff-avescatter" \
  bash tests/AveScatter/TestCoreShell.sh
```

The test uses supplied fine meshes and does not require gmsh. It creates
a temporary result directory and writes `comparison.json` there. It exits
unsuccessfully if a physical or accuracy check fails. Automake reports the
test as skipped when NumPy or SciPy is unavailable.

The validator can also inspect the saved fine result using

```bash
python3 tests/AveScatter/ValidateCoreShell.py \
  tests/AveScatter/reference/CoreShell.Fine.SIFlux.AVSCAT
```

The saved fine result passes. The medium and original-mesh results fail the
5% criterion at frequency 1. The original-mesh result also fails at frequency
0.3. They document convergence and support the separate regression comparison
using `Checklist.SIFlux`.

The historical incorrect absorption values are recorded above. The temporary
binaries, logs, duplicate outputs, and incorrect baseline files were discarded
during the migration. `CheckAveScatter.py` rejects empty files, missing or
extra rows, duplicate keys, NaN, and infinity. Checks of `ValidateCoreShell.py`
confirm that it rejects incomplete rows, negative material absorption,
incorrect scattering, and an impossible asymmetry factor.

The `CoreShell_Avscat` example contains the fine geometry and gmsh inputs.
The `BiSphere_Avscat` example preserves the two-sphere case and transformations.

This benchmark covers one concentric, isotropic, nonmagnetic material pair,
one radius ratio, and two frequencies. Passing it supports the implementation
for that tested regime. Other shapes, shell thicknesses, materials, and
resonances require their own convergence checks and references.
