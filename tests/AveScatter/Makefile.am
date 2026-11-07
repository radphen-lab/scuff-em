pkgdatadir = $(datadir)/scuff-em/tests/AveScatter
pkgdata_DATA = OmegaFile                     \
               Oblate.OmegaFile              \
               CoreShell.scuffgeo            \
               CoreShellMie.scuffgeo         \
               CoreShellMedium.scuffgeo      \
               BiSphere.scuffgeo             \
               BiSphere.trans                \
               LossySphere.scuffgeo          \
               Oblate.scuffgeo               \
               Sphere_327.msh                \
               Sphere_R0P5.msh               \
               Sphere_R1P0.msh               \
               Sphere_R0P5_Fine.msh          \
               Sphere_R1P0_Fine.msh          \
               Sphere_R0P5_Medium.msh        \
               Sphere_R1P0_Medium.msh        \
               oblate.msh                    \
               Checklist.SIFlux              \
               README.md                     \
               CoreShellValidation.md

referencedir = $(pkgdatadir)/reference
reference_DATA = reference/CoreShell.SIFlux.AVSCAT        \
                 reference/CoreShell.Medium.SIFlux.AVSCAT \
                 reference/CoreShell.Fine.SIFlux.AVSCAT   \
                 reference/BiSphere.SIFlux.AVSCAT         \
                 reference/LossySphere.SIFlux.AVSCAT      \
                 reference/Oblate.SIFlux.AVSCAT

pkgdata_SCRIPTS = TestAveScatter.sh TestCoreShell.sh CheckAveScatter.py ValidateCoreShell.py
EXTRA_DIST = $(pkgdata_DATA) $(reference_DATA) $(pkgdata_SCRIPTS)
TESTS = TestAveScatter.sh TestCoreShell.sh

SCUFF_AVESCATTER = $(abs_top_builddir)/applications/scuff-avescatter/scuff-avescatter
TESTS_ENVIRONMENT = srcdir='$(abs_srcdir)' SCUFF_AVESCATTER='$(SCUFF_AVESCATTER)'
