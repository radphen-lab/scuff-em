# SCUFF-EM with scuff-avescatter

Orientation- and polarization-averaged light scattering by randomly oriented
particles and particle clusters.

This repository is a research fork of **SCUFF-EM**, the open-source
boundary-element electromagnetic solver developed by M. T. Homer Reid and
the SCUFF-EM contributors, including Steven G. Johnson and Johannes Feist.
The underlying solver, supporting libraries, geometry and material framework,
and most applications come from their work. This fork focuses on the
`scuff-avescatter` application, its examples, numerical validation, and
implementation notes. We gratefully acknowledge the upstream authors and
contributors for making this work possible.
[Upstream repository](https://github.com/HomerReid/scuff-em),
[SCUFF-EM documentation](https://homerreid.github.io/scuff-em-documentation/),
[original authors and contributors](AUTHORS).

## Research paper

`scuff-avescatter` implements the orientation- and polarization-averaging
formulation presented in

**F. V. Ramirez-Cuevas, K. L. Gurunatha, I. P. Parkin, and I. Papakonstantinou.**
*Universal Theory of Light Scattering of Randomly Oriented Particles:
A Fluctuational-Electrodynamics Approach for Light Transport Modeling in
Disordered Nanostructures.* **ACS Photonics 9**(2), 672–681 (2022).
[Published paper](https://doi.org/10.1021/acsphotonics.1c01710),
[arXiv preprint](https://arxiv.org/abs/2110.14799v3).

The paper derives averaged scattering quantities for radiative-transfer
models of disordered materials. This application evaluates those quantities
using SCUFF-EM's surface-current BEM response and trace formulas. The average
over incident directions and polarizations is evaluated analytically for each
fixed geometry and frequency. The implementation draws on the matrix and
fluctuating-surface-current infrastructure of `libscuff` and `scuff-neq`.

The current averaged-scattering workflow targets finite particles and
clusters in a homogeneous, nonabsorbing host. The correspondence between the
paper, matrix conventions, and source code is documented in
[the theory and implementation notes](applications/scuff-avescatter/THEORY.md).

## What the application computes

For a geometry described by triangular surface meshes, `scuff-avescatter`
computes

- The orientation- and polarization-averaged absorption cross section.
- The averaged scattering cross section.
- The first angular moment of scattering, used to obtain the asymmetry factor.

The repository also contains the upstream applications for scattering under
specified illumination (`scuff-scatter`), thermal radiation (`scuff-neq`),
electrostatics (`scuff-static`), and other electromagnetic calculations.
Their documentation is available in the
[upstream application reference](https://homerreid.github.io/scuff-em-documentation/reference/TopLevel/).

## Installation

The commands below describe an Ubuntu/Debian command-line build, including
Ubuntu running under WSL on Windows. This repository contains the SCUFF-EM
source tree and builds `scuff-avescatter` with the other applications.
For platform-specific configuration and optional interfaces, consult the
[upstream installation guide](https://homerreid.github.io/scuff-em-documentation/reference/Installing/).

### Build dependencies

Install the compilers, GNU Autotools, and numerical libraries from an Ubuntu
or Debian terminal.

```bash
sudo apt-get update
sudo apt-get install git build-essential gfortran autoconf automake libtool \
  pkg-config libopenblas-dev liblapack-dev libreadline-dev
```

Python tools support validation and plotting, and gmsh supports mesh
generation and visualization. Install these optional tools with

```bash
sudo apt-get install python3 python3-numpy python3-scipy python3-matplotlib gmsh
```

The supplied examples and tests include meshes, so they can be run using
those meshes directly.

### Clone, build, and install

```bash
git clone https://github.com/radphen-lab/scuff-em.git
cd scuff-em

sh autogen.sh --prefix="$HOME/.local" --without-python --without-hdf5
make -j4
make install

export PATH="$HOME/.local/bin:$PATH"
```

Choose the parallel build job count to suit your machine. The installation
prefix can be changed with `--prefix=/your/installation/path`. Add the PATH
setting to your shell startup file to make it available in future sessions.
Applications are installed under `bin`, with examples and test data under
`share/scuff-em` inside the chosen prefix.

This build profile selects the command-line applications and text output.
`--without-python` skips the optional SCUFF-EM Python bindings, and
`--without-hdf5` skips HDF5 support. The validation and plotting scripts use
the system Python 3 installation separately from those bindings. OpenMP
support is enabled by default. Further configuration options are described
in the upstream installation guide linked above.

## Running a calculation

The main inputs are a `.scuffgeo` geometry file, its surface meshes and
material definitions, and a list of frequencies. Geometry and material
descriptions follow SCUFF-EM conventions.
[Geometry reference](https://homerreid.github.io/scuff-em-documentation/reference/Geometries/),
[material reference](https://homerreid.github.io/scuff-em-documentation/reference/Materials/).

From the repository root, enter the supplied core-shell example and run

```bash
cd examples/CoreShell_Avscat
scuff-avescatter --geometry CoreShell.scuffgeo --OmegaFile OmegaFile \
  --Threads 4 --FileBase CoreShell
```

The result is written to `CoreShell.SIFlux.AVSCAT`. Useful options include

| Option | Purpose |
|---|---|
| `--geometry file.scuffgeo` | Specify the geometry and material assignments |
| `--Omega value` | Evaluate a frequency supplied on the command line |
| `--OmegaFile file` | Read a list of frequencies |
| `--FileBase prefix` | Set the output filename prefix |
| `--TransFile file` | Evaluate listed geometric configurations |
| `--Threads N` | Set the solver's OpenMP thread count |

In vacuum, the frequency convention is `omega = 2*pi/lambda` when wavelength
is expressed in micrometres. For example, `--Omega 1` corresponds to a vacuum
wavelength of approximately 6.283 um. Output files are appended to across
runs. Use a new `--FileBase` for a separate calculation.

### Output columns

For the finite-geometry calculations described here, each output row contains

| Column | Meaning |
|---|---|
| 1 | Geometric transformation tag |
| 2 | Frequency in SCUFF units |
| 3 | Surface or object label from the geometry file |
| 4 | Averaged absorption cross section, in um^2 |
| 5 | Averaged scattering cross section, in um^2 |
| 6 | First scattering moment, often written `g*Csca`, in um^2 |

For the supplied multi-object examples, sum the surface contributions at each
frequency and transformation to obtain particle or cluster totals. The
asymmetry factor is the sum of column 6 divided by the sum of column 5,
where total scattering is nonzero. In the core-shell example, `OuterSphere`
identifies shell absorption and `InnerSphere` identifies core absorption.
The example READMEs explain the interpretation of their surface rows.

## Examples

| Example | Description |
|---|---|
| [CoreShell_Avscat](examples/CoreShell_Avscat/readme.md) | A concentric dielectric core-shell sphere with supplied fine meshes |
| [BiSphere_Avscat](examples/BiSphere_Avscat/readme.md) | Two dielectric spheres evaluated at two separations |
| [Spheroid_Avscat](examples/Spheroid_Avscat/readme.md) | Oblate and prolate spheroids with archived comparison data and a plotting script |

The `examples/` tree also contains the original SCUFF-EM examples for the
other applications.

## Tests and numerical accuracy

From the configured and built repository root, run the averaging-application
tests with

```bash
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 make -C tests/AveScatter check
```

`TestAveScatter.sh` checks oblate, bi-sphere, core-shell, and lossy-sphere
results against saved regression references. `TestCoreShell.sh` compares
separate shell and core absorption, total absorption, scattering, the first
scattering moment, and the asymmetry factor against an independent
coated-sphere Mie calculation. The scripts create temporary output
directories and print their locations. Automake reports a skipped test when
its required Python dependencies are unavailable.

For the supplied fine core-shell meshes at the two test frequencies, the
checked quantities agree with Mie theory within the selected 5% tolerance.
The largest recorded discrepancy is about 4.93% in the asymmetry factor.
These results document the tested geometry and discretization. Mesh
convergence is part of assessing accuracy for additional frequencies,
material contrasts, and geometries.
[Test instructions](tests/AveScatter/README.md),
[core-shell validation and convergence study](tests/AveScatter/CoreShellValidation.md).

## Repository layout

| Directory | Contents |
|---|---|
| `applications/scuff-avescatter/` | Averaged-scattering application and implementation notes |
| `applications/` | SCUFF-EM command-line applications |
| `libs/` | SCUFF-EM's BEM solver and supporting libraries |
| `examples/` | Example geometries, meshes, and run instructions |
| `tests/` | Numerical tests and reference data |
| `doc/` | SCUFF-EM documentation and technical material |

## Citation and acknowledgments

For research using the orientation-averaged scattering formulation, please
cite the Ramirez-Cuevas et al. paper given above. Please also acknowledge the
SCUFF-EM software and the relevant upstream methods. One SCUFF-EM reference is

**M. T. Homer Reid and S. G. Johnson.** *Efficient Computation of Power,
Force, and Torque in BEM Scattering Calculations.* **IEEE Transactions on
Antennas and Propagation 63**(8), 3588–3598 (2015).
[Published paper](https://doi.org/10.1109/TAP.2015.2438393),
[SCUFF-EM citation guidance](https://homerreid.github.io/scuff-em-documentation/reference/TopLevel/#citing-scuff-em).

The upstream project overview and contributor acknowledgments are available
in [README](README) and [AUTHORS](AUTHORS).

## License

SCUFF-EM is distributed under the GNU General Public License. The repository
retains the upstream copyright and licensing notices. See [LICENSE](LICENSE),
[COPYING](COPYING), [COPYRIGHT](COPYRIGHT), and the notices in individual
source files.
