# Orientation-averaged light scattering in SCUFF-EM

## Project scope

SCUFF-EM is a C++ surface-integral-equation boundary element method (BEM) suite. Its core, `libs/libscuff`, represents electric and magnetic surface currents in the RWG basis on triangular meshes. It also provides the matrix assembly and Green function integrals used by the command-line tools. This fork builds `scuff-avescatter` alongside `scuff-scatter`, `scuff-neq`, and the other applications. The root README describes the averaging application, while `doc/docs/reference/TopLevel.md` describes the full suite.

`scuff-avescatter` is the user's primary research focus and implements their 2022 paper. The application computes direction- and polarization-averaged absorption and scattering cross sections and the first angular moment of scattering. The usual inputs are a fixed `.scuffgeo` geometry with surface meshes and material assignments, and frequencies given by `--Omega` or `--OmegaFile`.

## From the analytical average to a BEM trace

For a fixed geometry and frequency, the scattered response is linear in the incident field, while the observables are quadratic. Average the outer product of the incident field over all plane-wave directions and both transverse polarizations before applying the scattering operator. In the paper's electric-field notation, Eq. 2 gives an incident correlation proportional to `Asym(G0)`, where `G0` is the exterior Green operator and `Asym(A) = (A - A†)/(2i)`. The extra incident-direction component in the asymmetry calculation gives a correlation proportional to `Sym(∂j G0)` in Eq. 3, where `Sym(A) = (A + A†)/2`. The derivative follows from differentiating the plane-wave phase. Main-text Eqs. 4a–c contract these kernels with the scatterer's response. Supporting Information Eqs. S8 and S11 derive the two correlations.

In the six-component surface-current BEM formulation, the `ik` prefactor of Eq. S12 changes the decomposition. The exterior source matrix in Eqs. S13a–b uses `Sym(G0)` in the BEM convention. Equation S13c uses `Asym(∂j G0)` on both sides of the response and sums over the three Cartesian directions. The absorption observable uses the interior material operator and the scattering observable uses the exterior radiation operator. The optical derivation extends to a linear, homogeneous, nonabsorbing host using the host wavenumber. For signs and prefactors, distinguish the paper's operators, the surface-current BEM convention, and SCUFF's stored matrix convention.

Rodríguez et al. (2013), §II.C–D, is the source of the fluctuating-surface-current trace architecture. For a projected source correlation `S` and a BEM matrix `M`, write `W = M⁻¹`. The response propagates the source correlation to `D = W S W†`, and the observable is a trace contraction with `D`. In thermal transfer the source comes from fluctuations inside a material body: Rodríguez Eq. 16 gives `⟨ss†⟩ = −(2/π) sym G₁`, and the single-body emissivity of Eq. 21 is `Tr[(sym G₁) W† (sym G₀) W]` with `W = (G₁ + G₀)⁻¹`. In `scuff-avescatter` the roles are swapped: the *exterior* block is the source, which encodes isotropic incident illumination, while the interior block (absorption) or the exterior block again (scattering) is the observable. The code computes optical cross sections without thermal occupation factors or frequency integration. Equation 5 relates the spectral absorption factor to the average cross section by `π²/k²` in the lossless host.

Abraham Ekeroth et al. (2017) gives a related covariance-and-trace formulation in a volume discrete-dipole basis, including directional emissivity and absorptivity. It provides a comparison for the physical idea. The direct implementation lineage is the surface-current route through `scuff-neq` and the energy/momentum transfer routines.

### Paper-to-code convention map

Let `Tₑ` and `Tᵢ` be SCUFF's exterior and interior self blocks after `UndoSCUFFMatrixTransformation` (scale `EE` by `ZVAC`, `NN` by `−1/ZVAC`, `EM` by `−1`). The code's BEM matrix is the Rodríguez sum `M = Tₑ + Tᵢ (+ U blocks)`. The SI writes `W = [Gin − G0]⁻¹`, which matches `M⁻¹` if the paper's operators correspond to `G0 ↔ −Tₑ` and `Gin ↔ Tᵢ`. Under that map, every kernel in S13a–c appears either squared (`G0 … G0`, `∂G0 … ∂G0`) or as the product `(−G0)(−Gin)`. The code's choice of `+Tₑ` and `+Tᵢ` for both factors therefore reproduces S13a–c exactly, with no hidden sign flip.

## Prefactor audit

This checks the arithmetic of the source against the SI. The independent core-shell benchmark described below also checks the normalized observables numerically. The source matrix is scaled by `RYTOVPF = −4/π`, every trace is multiplied by `−1/2`, and the result is divided by `kPI2 = |khost/π|²`.

| Output | SI target | Code chain | Net |
|---|---|---|---|
| `⟨Cabs⟩` | `(2π/k²) Tr[Sym G0 W† Sym(−Gin) W]` (S13a) | `(−½)(−4/π)·(π²/k²)·Re Tr[Tᵢ · W Sym(Tₑ) W†]` | `2π/k²` ✓ |
| `⟨Csca⟩` | `(2π/k²) Tr[Sym G0 W† Sym G0 W]` (S13b) | same, with `Tₑ` as observable | `2π/k²` ✓ |
| `⟨µsca Csca⟩` | `(2π/k⁴) Σₗ Tr[Asym ∂ₗG0 W† Asym ∂ₗG0 W]` (S13c) | `(−½)(−4/π)·(π²/k²)·|cR|²·Σ Tr[A W A W†]`, with `A` built from EMTPFT force kernels | `2π/k⁴` ✓ |

For the asymmetry row, the force kernels in `GetForceIntegrals` carry `TENTHIRDS·MuR·∂G`, while the power kernel in libscuff `EMTPFT.cc` carries `Omega·MuR·G`. Each derivative factor is therefore `(10/3)/ω` times the power-normalized `∂G`. Two factors give `(10/3)²/ω²`, and `cR = 1/(TENTHIRDS·n)` contributes `(0.3/n)²`. The product is `1/(ωn)² = 1/k²` in SCUFF units, where ω = 1 means `k = 1 µm⁻¹`. This `1/k²` turns the cross-section `2π/k²` into the SI's `2π/k⁴`. The `|cR|²` factor is a unit conversion from force-normalized to power-normalized kernels. For a real host index `n` it is exact. `kPI2` and `|cR|²` use moduli, so an absorbing host would silently use `|n|²`. The theory itself assumes a nonabsorbing host.

The `(−½)(−4/π) = 2/π` combination is the normalization of the paper's Eq. 5, `Φabs = (2/π)Tr{…}` with `⟨Cabs⟩ = (π²/k²)Φabs`. It is four times the `1/(2π)` in Rodríguez Eq. 21 (`−¼` trace × `−2/π` source), a convention difference already present in `scuff-neq`'s `ComputeDRMatrix`.

## Source map

| Where | Role |
|---|---|
| `scuff-avescatter.cc`, `main` | Reads geometry and frequencies, creates `SNEQData`, and calls `NeqScattering` once per requested frequency. The active method is EMT, and it is the default. `--OPFT`, `--DSI*`, and `--OmegaKBFile` are parsed but only log "not supported". `--Omega` accepts at most 10 values (`MAXFREQ`), while `--OmegaFile` has no such limit. `--Threads` sets the OpenMP thread count. |
| `NeqScattering.cc`, `NeqScattering` | Assembles `TExt` and `TInt` once per frequency, reusing blocks for identical "mate" surfaces. For each `--TransFile` configuration it assembles the `U` blocks and stamps `M = TExt + TInt + U`. It copies every block that touches region 0 into `T0`, undoes the SCUFF transformation, LU-factorizes once, and calls `AverageScatteringFlux`. The orientation average itself is analytic. |
| `NeqScattering.cc`, `HermitianPart` / `AdvancedDRMatrix` | Builds `DR = W · (−4/π)·Sym(T0) · W†` with `LUSolve → Adjoint → LUSolve`. This mirrors `scuff-neq/WriteFlux.cc`, `ComputeDRMatrix`, except that the source is the exterior block `T0` instead of one body's `TInt`. |
| `NeqScattering.cc`, `ScaTraceFlux` / `AbsTraceFlux` / `MatrixTraceFlux` | `MatrixTraceFlux` computes `Re Tr[G̃ · DR]`, applying the SCUFF un-scaling to `G̃` on the fly. Because `DR` is Hermitian, `Re Tr[G·DR] = Tr[Sym(G)·DR]`, which is why the observable needs no explicit `Sym`. `ScaTraceFlux` sums the exterior blocks `T0(d, s)` for all sources `s` touching region 0. `AbsTraceFlux` sums every block of the region-`r` Green operator, meaning the self blocks of all surfaces bounding `r` plus their `U` couplings in both orders. |
| `NeqScattering.cc`, `gradT0_Matrix` / `StoreMatrix` / `GetForceIntegrals` | A trimmed copy of libscuff `EMTPFT.cc`, including its `SCUFF_EMTPFTI_EHDERIVATIVES2` branch and its region-`Sign` logic. It overwrites `T0` with an antisymmetric matrix of exterior force kernels for one axis. Entries are `KK = +i·Q_KK`, `KN = NK = −Q_KN−NK`, `NN = −i·Q_NN`, and the transposed entries take the opposite sign. |
| `NeqScattering.cc`, `iAntiHermitianPart` | Despite the name, it returns `Asym(A) = (A − A†)/(2i)` after un-scaling. On a single surface it uses the closed forms valid for an antisymmetric matrix: `−i·Re` for KK and NN, `Im` for KN and NK, and zero for the same edge. Across different surfaces it uses the general formula. |
| `NeqScattering.cc`, `AverageScatteringFlux` | Order of operations: first the `Sym(T0)` source with per-surface `Csca` and `Cabs`; then, for x, y, z, rebuild `T0 ← ∂T0`, build `DR`, and contract with `Asym(∂T0)`. It applies the normalizations in the prefactor audit. |
| `CreateSNEQData.cc`, `WriteSIFluxFilePreamble` | Defines the columns and sets the output suffix `.SIFlux.AVSCAT`. Columns 4–6 are `Cabs`, `Csca`, and `µsca·Csca`, all in µm². One row is written per surface. Column 3 uses `G->Surfaces[nsd]->Label`, matching the surface-label output in `scuff-scatter/OutputModules.cc`. |

One calculation means one BEM assembly and factorization per frequency and configuration. It uses four matrix right-hand-side solves of size `N` (the scattering kernel plus three derivative kernels), each requiring two `LUSolve` calls. There is no orientation or polarization sampling loop.

### Per-surface output and clusters (SI §1.4–1.5)

The per-row values are the SI Eq. S26 decomposition. For destination surface `n`, column 5 is `Σₘ Tr[S_nm Sym G0_mn]`, combining self and interference scattering. Column 6 is the analogous contribution to `µ·C`. Column 4 is the absorption in the interior region of surface `n`. For a core–shell geometry this is the per-material absorption, because the shell row covers the shell material and the core row covers the core. Cluster totals (Eq. S25) come from summing rows. The total `⟨µsca⟩ = Σ col6 / Σ col5`. The SI's `⟨µⁿsca⟩` (Eq. S26c) is a different, `N`-scaled per-object quantity: it is not `col6/col5` of one row.

### Cost and memory

`M`, `DRMatrix`, `T0`, and the temporary `asymdG` are dense `N×N` complex matrices, where `N = TotalBFs` (two per interior edge). Together they take about `64 N²` bytes, plus the `TExt`, `TInt`, and `U` blocks. Runtime is one `O(N³)` LU plus eight `N`-column `O(N³)` triangular solves. `gradT0_Matrix` runs once per axis, and each pass calls `GetGCMatrixElements` with all three force components but keeps only one, so it does three times the necessary force-integral work. Assembling the three derivative matrices in one pass would remove that cost at the price of more memory.

## References and retrieval

| Reference | Stable source | Zotero item / BibTeX key |
|---|---|---|
| Ramirez-Cuevas, Gurunatha, Parkin, and Papakonstantinou, *ACS Photonics* 9, 672–681 (2022) | [Journal article](https://doi.org/10.1021/acsphotonics.1c01710), [arXiv with supplement](https://arxiv.org/abs/2110.14799v3) | `YRFHPYPY` (in "non-reciprocal photonics"), duplicate `7X55H62L` / `Ramirez-Cuevas2022` |
| Rodriguez, Reid, and Johnson, *Phys. Rev. B* 88, 054305 (2013) | [Journal article](https://doi.org/10.1103/PhysRevB.88.054305), especially Eqs. 7–21 | `6IJHNA9I` (in "non-reciprocal photonics"), duplicate `J5L9WGT5` / `Rodriguez2013b` |
| Abraham Ekeroth, García-Martín, and Cuevas, *Phys. Rev. B* 95, 235428 (2017) | [Journal article](https://doi.org/10.1103/PhysRevB.95.235428), [public full text](https://arxiv.org/abs/1702.04273v2) | `I9ESJW5R` / `AbrahamEkeroth2017` |

The lowercase "non-reciprocal photonics" Zotero collection, distinct from "Non-reciprocal photonics", holds the working copies. Item `YRFHPYPY` has attachment `NLISFXH2`, the published article, and attachment `GKDBKI77`, the Supporting Information, despite its Abraham Ekeroth filename. Item `6IJHNA9I` has attachment `3RZ63QNN`, the Rodriguez article. The same collection also contains Reid's "libscuff Implementation and Technical Details" (`35JQAAUT`, attachment `XE3I334A`), useful for SCUFF's matrix conventions. It also has Krüger et al. 2012, "Trace formulas for nonequilibrium Casimir interactions…" (`8KMEUSRV`, attachment `UYTN234I`), and Homer Reid's MIT thesis (`FN4U7BJG`). The older duplicates have attachments `X6ZV2B7D` (article), `QZEDQFKA` (SI), and `I84B39JR` (Rodriguez). From WSL, attachments are at `/mnt/c/Users/frami/Zotero/storage/<attachment key>/`. Copy `zotero.sqlite` before querying it, because Zotero locks it while running. The SI PDF's equations extract legibly with `pdftotext -layout`, but symbols such as `𝕋` and `𝔾` are lost in the main article's extraction.

## Validation status

- The paper compares against T-matrix results for a silver sphere dimer and oblate/prolate spheroids. The repository includes spheroid meshes, frequency lists, and historical outputs under `examples/Spheroid_Avscat/`. Its 50 frequencies cover λ = 0.25–10 µm (`λ = 2π/ω` in SCUFF units).
- The archived outputs (`*.AVSCAT.EMTPFT`, 2021) were compared with `*_exact.dat` without rerunning the code. In the size range the SI calls reliable (`Lmax/λ < 2.37`; 31 oblate and 28 prolate points), the maximum relative errors are below 1.7% for `Cabs`, 2.8% for `Csca`, and 1.0% for `µsca`. Across all 50 points the median errors are about 1.5–2.3%. The maximum errors reach 50–60% for `Cabs` near λ ≈ 0.35 µm (`Lmax/λ` ≈ 7–8), consistent with the SI's mesh-resolution explanation.
- The archived files label column 6 "Rad. Pressure Cross Section", but the values are `µsca·Csca`: dividing by column 5 matches the exact asymmetry parameter to within 1%. The current source labels it correctly.
- The spheroid comparison script prefers `.SIFlux.AVSCAT` and falls back to the historical `.AVSCAT.EMTPFT` files. It reads only the numerical frequency and observable columns, so it accepts both surface labels and historical numeric surface identifiers.
- The default top-level test build registers Mie, Fresnel, and AveScatter tests. `tests/AveScatter/TestAveScatter.sh` compares oblate, bi-sphere, core-shell, and lossy-sphere runs against corrected references using a strict checklist checker. `tests/AveScatter/TestCoreShell.sh` compares shell and core absorption, scattering, and asymmetry against an independent coated-sphere Mie solution in `ValidateCoreShell.py`. Its geometry, accuracy criterion, convergence data, and reproduction commands are recorded in `tests/AveScatter/CoreShellValidation.md`. The independent Mie test requires NumPy and SciPy and is reported as skipped if those dependencies are unavailable. Runnable examples are under `examples/CoreShell_Avscat/` and `examples/BiSphere_Avscat/`.
- On 2026-11-07, independent core-shell validation exposed an existing dimension error in `MatrixTraceFlux`. It used `DRMatrix->NC` for the observable's source dimension, which should be `GDest->NC`. Rectangular off-diagonal blocks could therefore access entries outside the matrices and produce unphysical shell absorption. The corrected source has been rebuilt. The historical incorrect values are recorded in the validation report. The migrated test references use the corrected results.
- After that correction, a three-mesh study of the concentric core-shell case (radii 0.5/1 um, core epsilon `4+1i`, shell epsilon `10+0.5i`, vacuum host) shows decreasing errors against Mie theory. The fine mesh has 1258 outer and 708 inner triangles. At SCUFF frequencies 0.3 and 1, shell absorption, core absorption, total absorption, scattering, the first scattering moment, and g all agree within the selected 5% tolerance. The largest remaining discrepancy is 4.93% in g at frequency 1. This validates that tested regime, with the reported discretization error.

## Implementation notes and open questions

- When stamping `T0`, an off-diagonal `U` block is kept if surface `nsp` touches region 0, without checking that region 0 is the region `ns` and `nsp` share. For closed surfaces this is safe because `RWGGeometry` topologically sorts nested surfaces with the container first. Geometries with open surfaces or explicit `REGION` assignments could insert an interior-region coupling into the exterior source kernel.
- PEC surfaces have `RegionIndices[1] = −1`, so their absorption is 0 and `TInt` is never assembled for them.
- The sign of `µsca` does not depend on the sign convention of `∂G`, because the derivative kernel appears twice.
