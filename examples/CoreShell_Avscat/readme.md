# Orientation-averaged scattering by a core-shell sphere

The core radius is 0.5 um with relative permittivity `4+1i`. The shell outer
radius is 1 um with relative permittivity `10+0.5i`. The host is vacuum and
the materials have relative permeability one. `OmegaFile` contains SCUFF
frequencies 1 and 0.3, corresponding to vacuum wavelengths 6.283 and 20.944 um.

From this directory, run

```bash
scuff-avescatter --geometry CoreShell.scuffgeo --OmegaFile OmegaFile --Threads 4
```

The output is `CoreShell.SIFlux.AVSCAT`. Column 4 gives absorption, column 5
gives scattering, and column 6 gives g times scattering, in um^2. Column 3
contains the surface label from the geometry file. The `OuterSphere` row
contains shell absorption and the `InnerSphere` row contains core absorption.
Sum each column over the two surfaces for the particle totals. The particle asymmetry
factor is the total of column 6 divided by the total of column 5.

The supplied meshes have 1258 outer and 708 inner triangles. For these inputs,
the six observables checked by `tests/AveScatter/TestCoreShell.sh` agree with
coated-sphere Mie theory within 5%. The largest error is 4.93% in g at frequency
1. See `tests/AveScatter/CoreShellValidation.md` for the convergence study and
the independent reference calculation.

The supplied meshes make this example runnable without gmsh. To regenerate
them using gmsh, run

```bash
gmsh -2 Sphere_R0P5.geo -o Sphere_R0P5.msh
gmsh -2 Sphere_R1P0.geo -o Sphere_R1P0.msh
```

Changing the meshing parameters changes the discretization error and requires
another convergence check.
