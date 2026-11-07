# Orientation-averaged scattering by two dielectric spheres

Each sphere has radius 1 um and relative permittivity 10. The host is vacuum.
The lower sphere is centered at the origin and the upper sphere is initially
centered at z = 3 um. `BiSphere.trans` specifies center separations of 3 and
4 um, corresponding to surface gaps of 1 and 2 um. `OmegaFile` contains SCUFF
frequencies 1 and 0.3.

From this directory, run

```bash
scuff-avescatter --geometry BiSphere.scuffgeo --OmegaFile OmegaFile \
  --TransFile BiSphere.trans --Threads 4
```

The output is `BiSphere.SIFlux.AVSCAT`, with two surface rows for each
frequency and transformation. Columns 4, 5, and 6 give absorption, scattering,
and g times scattering in um^2. Sum the two surface rows to obtain cluster
totals. The cluster asymmetry factor is the total of column 6 divided by the
total of column 5. Individual scattering rows include interference terms.

This example exercises identical-object block reuse and transformed
interactions. `tests/AveScatter/TestAveScatter.sh` checks the numerical output
against a saved regression reference. The reference records solver behavior
for the supplied mesh. It has no independent two-sphere analytical benchmark.
The permittivity is real, so physical absorption is zero. Small computed
absorption values reflect discretization and numerical integration error.

`Sphere_327.msh` is supplied with the example. `Sphere.geo` is its original
gmsh geometry source. Regenerating a mesh with a different gmsh version may
change the triangulation and its regression values.
