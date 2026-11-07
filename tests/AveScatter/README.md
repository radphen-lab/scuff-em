# Orientation-averaged scattering tests

This directory follows the Mie and Fresnel test layout. Automake installs the
scripts, input meshes and geometries, checklists, and `reference/` data under
`share/scuff-em/tests/AveScatter`. The two scripts are registered with
`make check` through `tests/Makefile.am`.

## Running

From a configured and built source tree, run

```bash
make -C tests/AveScatter check
```

For an installed executable on PATH, run either script from this directory

```bash
./TestAveScatter.sh
./TestCoreShell.sh
```

Set `SCUFF_AVESCATTER` to select another executable. Set `OMP_NUM_THREADS` to
select the solver thread count, which defaults to 4. The regression script
requires Python 3. The independent Mie test additionally requires NumPy and
SciPy. Automake reports a missing Python dependency as a skipped test.
The supplied meshes remove the need to install gmsh to run the tests.

Each script creates and prints a temporary result directory. Outputs and logs
remain there for inspection, including after a failure. Input files and saved
references remain in the source or installed test directory.

## Regression cases

`TestAveScatter.sh` computes each case and compares the result with its
`reference/*.SIFlux.AVSCAT` file using `Checklist.SIFlux`. Transform,
frequency, and surface label identify each row. Column 3 contains the label
from the geometry's OBJECT or SURFACE declaration, as in scuff-scatter.
All three result columns are checked.
The relative tolerance is 1e-8 with a 1e-10 um^2 absolute floor. The floor
covers cancellation differences in the nearly zero absorption of the
lossless bi-sphere case.

| Case | Coverage |
|---|---|
| Oblate | One lossy dielectric surface at three frequencies |
| BiSphere | Identical-body reuse, inter-surface coupling, and two separations |
| CoreShell | Nested surfaces with unequal matrix-block dimensions |
| LossySphere | Single absorbing sphere at two frequencies |

The references preserve the corrected solver's printed results. The CoreShell
reference uses the original coarse mesh and includes the rectangular trace
correction. Its physical accuracy is assessed separately on the fine mesh.
Repeated runs with `OMP_NUM_THREADS=1` and `OMP_NUM_THREADS=4` check agreement
with the same saved references.

`CheckAveScatter.py` uses the repository's KEY/DATA checklist format and fails
on numerical differences, empty files, missing or extra rows, duplicate keys,
and non-finite values. It is used here because the existing `CheckSCUFFData`
program prints mismatch counts but currently returns success even when data
do not match. The new checker returns a failing process status for mismatches.

## Independent core-shell reference

`TestCoreShell.sh` computes `CoreShellMie.scuffgeo` on the fine mesh and calls
`ValidateCoreShell.py`. It checks shell absorption, core absorption, total
absorption, scattering, the first scattering moment, and g against an
independent coated-sphere Mie solution, with a 5% discretization tolerance.
It also checks physical signs and the bound on g.

The validator checks its reference through homogeneous-sphere, quasistatic,
lossless-material, and multipole-truncation limits. Its saved numerical
validation is documented in [CoreShellValidation.md](CoreShellValidation.md).
The original, medium, and fine BEM outputs in `reference/` document mesh
convergence. The independent test evaluates Mie theory at runtime.

The corresponding input geometries are `CoreShell.scuffgeo`,
`CoreShellMedium.scuffgeo`, and `CoreShellMie.scuffgeo`. Each uses the same
`OmegaFile`, and all three sets of meshes are supplied. These inputs allow
the convergence study to be repeated without regenerating meshes.

To inspect the saved fine output directly, run

```bash
python3 ValidateCoreShell.py reference/CoreShell.Fine.SIFlux.AVSCAT
```

Runnable examples are distributed under `examples/CoreShell_Avscat` and
`examples/BiSphere_Avscat`.
