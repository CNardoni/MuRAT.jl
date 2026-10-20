# MuRAT.jl

Julia implementation of **MuRAT 4.0** (Multi-Resolution seismic Attenuation
Tomography, MATLAB, https://github.com/LucaDeSiena/MuRAT, version of
2026-10-06): peak-delay, coda attenuation (Qc) and coda-normalization (Q)
tomography from SAC data.

The code follows the MATLAB functions step by step, uses the same input fields,
fills the same result structure (`Murat.input`, `Murat.rays`, `Murat.PD`,
`Murat.Qc`, `Murat.Q`) and writes the same output tables. On the Mount St. Helens
example it reproduces the MATLAB results (see *Validation*). Every function has a
docstring that names its MATLAB counterpart (`?murat_qc` in the Julia REPL).

## Installation

In the Julia REPL, press `]` to enter the package mode, then:

```
pkg> add https://github.com/CNardoni/MuRAT.jl
pkg> test MuRAT
```

`add` installs MuRAT.jl and all the packages it needs for the computation, the
text tables and the ParaView files.

### Optional features

MuRAT.jl runs the complete analysis (measurements, inversions, tests, text
tables, ParaView files) without any other package. Four packages add optional
features; install only those you need:

| Package | You need it if | Without it |
|---------|----------------|------------|
| CairoMakie | you want the figures (PNG) | no figures (message at the end of the run) |
| MAT | you set `SaveMatFiles = 1` (`.mat` files for MATLAB) | no `.mat` files (message at the end of the run) |
| NCDatasets | your input uses `namev = "LITHO1.0.nc"` | the run stops at the start with an error asking for it |
| GeophysicalModelGenerator | you want the topography on the vertical sections (`topo` keyword) | figures without topography |

Install them with, e.g., `pkg> add CairoMakie MAT`, and load them with `using`
before the run; MuRAT.jl then enables the feature automatically (Julia package
extensions):

```julia
using MuRAT, CairoMakie, MAT     # figures and .mat files
murat = run_murat(input)
```

## Running

A run needs an input file (a `MuratInput`, see *Input*) and is started from the
folder that contains the data folder and `velocity_models/`: the paths in the
input are relative to it, as in MATLAB, and the results go to the `label` folder
there.

```julia
using MuRAT
using CairoMakie, MAT                 # optional: figures, .mat files

cd("/path/to/my/project")             # contains the SAC folder and velocity_models/
input = include("my_input.jl")        # your MuratInput
murat = run_murat(input)
```

**Test data sets.** `examples/` contains the three test data sets of MuRAT
(`sac_MSH`, `sac_Toba`, `sac_Romania`, with their velocity models in
`examples/velocity_models/`; copied from the MuRAT repository), their inputs and
a run script for each. A test run works from any folder; the results go to the
folder you are in:

```julia
using MuRAT
cd("/path/where/results/should/go")
include(joinpath(pkgdir(MuRAT), "examples", "run_MSH.jl"))   # or run_Toba.jl, run_Romania.jl
```

`pkgdir(MuRAT)` gives the folder where MuRAT.jl is installed. The run scripts load
CairoMakie (and NCDatasets for Romania), which must be installed (see *Optional
features*).

**Parallel processing.** The traces are processed with Julia threads (MATLAB:
`parfor`). Start Julia with them, e.g. `julia -t 6`; `Threads.nthreads()` shows
how many are available. In VS Code, set the number of threads with the setting
`julia.NumThreads`. Julia cannot add threads to a running session, so the input
field `workers` (MATLAB: size of the parallel pool) is kept only for
compatibility with the MATLAB input files.

The stages can also be run one by one, like the sections of `MuRAT.m`:

```julia
m = murat_checks(input)      # folders, SAC headers, grid, velocity models
murat_data(m)                # measurements, rays, kernels, data selection
murat_inversion(m)           # peak delay, Qc, Q, checkerboard and spike tests
murat_summary(m)             # average Qc per frequency, Qc-frequency fit
murat_paraview(m)            # ParaView files
murat_matfiles(m)            # .mat files (needs `using MAT`)
murat_plot(m; topo = Topo)   # figures (needs `using CairoMakie`)
```

A finished run can be reloaded and plotted again without recomputing it:

```julia
m = murat_load("label/temp/Murat_inverse.jld2")
murat_plot(m; topo = Topo)
```

Runs saved by an earlier version of the package also load: input fields added
since then take their default values.

## Input

`MuratInput` has the fields of `Murat.input` in the MATLAB input files (see
`examples/input_MSH.jl`, `input_Toba.jl` and `input_Romania.jl`, converted from
the MATLAB test inputs). To convert an input
file, copy the values and change:

| MATLAB                              | Julia                               |
|-------------------------------------|-------------------------------------|
| `Murat.input.end = [...]`           | `ending = [...]` (`end` is reserved)|
| `[]`                                | `nothing`                           |
| `'text'`                            | `"text"`                            |
| `[3 6 12 18]`                       | `[3, 6, 12, 18]`                    |
| `Murat.input.originTime = 'o'`      | `originTime = "o"` (header name)    |

Fields that exist only in the Julia version (all optional):

| Field | Default | Meaning |
|-------|---------|---------|
| `velocityFolder` | `"velocity_models"` | folder of the velocity model file |
| `kernelType` | `"differential"` | `"differential"`: kernel at the coda start minus kernel at the coda end (as MuRAT 4.0); `"middle"`: one kernel at the middle of the coda window |
| `plotInterpolation` | `3` | refinement factor of the result and checkerboard maps (1 = flat cells on the grid) |
| `SaveMatFiles` | `0` | `1` writes the `.mat` files listed in *Outputs* |

`workers` is read but only as the number of threads wanted (see *Running*).

### Topography

The vertical sections of all maps and ray figures can show the topography
profile. Pass it to `run_murat` (or `murat_plot`), e.g. from
GeophysicalModelGenerator over the grid limits:

```julia
using GeophysicalModelGenerator
lonlim = sort([input.origin[2], input.ending[2]]) .+ [-0.05, 0.05]
latlim = sort([input.origin[1], input.ending[1]]) .+ [-0.05, 0.05]
Topo = import_topo(lon = lonlim, lat = latlim, file = "@earth_relief_15s")
murat = run_murat(input; topo = Topo)
```

`import_topo` downloads the topography (an internet connection is needed).

## Outputs

Inside the `label` folder:

- `TXT/peakdelay_<f>_Hz.txt`, `TXT/Qc_<f>_Hz.txt`, `TXT/Q_<f>_Hz.txt`: same
  columns as MATLAB (model, hit count, checkerboard and spike input and output);
- `Results/MuRAT_<f>_Hz.vts`: all models and tests for ParaView;
- `Tests/LCurve/*.txt`: L-curve values (residual norm, solution norm,
  regularization parameter);
- `Tests/DataHeaders.txt`: table of the SAC headers;
- `temp/Murat_checks.jld2`, `temp/Murat_forward.jld2`: the run after the checks
  and after the data stage; `temp/Murat_inverse.jld2`: only the inversion results
  (compressed; `murat_load` on the inverse file adds them to the forward run, so
  keep the two files together);
- figures (PNG) in `Rays`, `Tests`, `Results`, `Checkerboard`, `Spike` and
  `Kernels`, with `using CairoMakie`;
- with `SaveMatFiles = 1` and `using MAT`:
  - `matfiles/*.mat`: per frequency `PD_inter`, `PD_not_inter`, `Qc_inter`,
    `Qc_not_inter`, `Q`, `hitmap`, `Qc_check_in`, `Qc_check_out`, with the file
    and variable names of `Murat_plot.m`;
  - `Murat_summary.mat`: the `Murat` structure (variable `Murat`) without its
    largest arrays (inversion matrices, ray paths), whose names are listed in
    `Murat.skippedFields`.

In Julia the results are in the same places as in the MATLAB structure
(`murat.Qc.inverseQc`, `murat.Qc.modvQc`, `murat.PD.peakDelay`, `murat.Q.modvQ`,
`murat.rays.locationsDeg`, ...). The derived settings that MATLAB adds to
`Murat.input` (grid, velocity models, file list, headers) are in `murat.grid`.

## Differences from MATLAB MuRAT 4.0

### Computation

- **Inversion methods**: only `"Tikhonov"` (FISTA with non-negativity) is
  implemented. `"Particle"`, `"Annealing"` and `"Global"` (MATLAB Global
  Optimization Toolbox) stop with an error; note that `Murat_inputMSH.m` uses
  `"Particle"`, the Julia example uses `"Tikhonov"`.
- **LITHO1.0** velocity models are read with NCDatasets (`using NCDatasets` before
  the run); the profile matches MATLAB `Murat_readLithos1` to single precision.
- **Kernels**: `kernelType = "differential"` is the MATLAB 4.0 computation; the
  `"middle"` option has no MATLAB 4.0 counterpart. With `"differential"` only the
  two kernels needed are computed (no kernel figure inside the data loop).
- **Random numbers**: the FISTA solver starts its step-size estimate from random
  numbers (seed 0, `murat_inversion(m; seed = 0)`), which cannot reproduce
  MATLAB's; Q maps differ by a few percent (see *Validation*).
- **3D `griddata`** uses Qhull's Delaunay triangulation; on a regular grid the
  triangulation is not unique, so rescaled velocity models differ slightly from
  MATLAB at some nodes.
- **Parallel processing**: Julia threads instead of a MATLAB parallel pool
  (`workers`, see *Running*).

### Files and folders

- An existing `label` folder is moved to `label_old` (MATLAB deletes or renames
  it). Intermediate results are compressed JLD2 files in `label/temp` instead of
  `.mat` files in `./temp`, and `Murat_inverse.jld2` contains only the inversion
  results instead of repeating the whole structure. To redo only the inversion,
  e.g. with other damping values, load `Murat_forward.jld2` (or the inverse file),
  change `m.input` and call `murat_inversion(m)`, `murat_summary(m)`.
- `DataHeaders` is a tab-separated text file instead of an Excel file.
- `.mat` files only with `SaveMatFiles = 1`. `matfiles/` also contains the
  `hitmap` and `Qc_check_in/out` files (not saved by MATLAB 4.0), and
  `Murat_summary.mat` replaces the full `Murat.mat` (several GB for large data
  sets). No `.fig` files.
- No `Results/Parameter` folder (the parameter-separation map is disabled in
  MATLAB and not implemented).

### Figures

Made with Makie; same content as MATLAB, different layout:

- **Maps** (results, checkerboards, spikes, hit maps, kernels): three 2D panels
  instead of 3D slices — horizontal section on the left with the aspect of the
  real distances, west-east and south-north sections on the right with vertical
  exaggeration ×2 and **depth** (km, positive down) instead of altitude.
- **Section positions**: results and checkerboards are taken at the exact
  `sections` positions and interpolated `plotInterpolation` times finer (MATLAB:
  fixed `divi = 5` and Gaussian smoothing). Hit maps and spikes are flat cells on
  the grid planes nearest to `sections`; kernel sections pass through the kernel
  maximum.
- **Topography** profile on all vertical sections when `topo` is given.
- **Hit map**: `hot` colour map saturated at 50 rays (as `Murat_hitmap.m`).
- **Qc test**: only the measurements within mean ± 1 std, with their 500-point
  moving mean (MATLAB 4.0 plots all measurements and a Picard plot below).
- **Peak-delay test**: logarithmic axes labelled in seconds.
- **Rays**: map view and projections on the west-east and south-north planes;
  no coastline or location panels.
- **L-curves** are saved for every inversion (PNG and TXT).

### MATLAB behaviour kept on purpose

To reproduce the MATLAB results, these details are kept as they are: MATLAB's
colon operator (`a:d:b`) and `round` (halves away from zero); the ray-cell
assignment in `Murat_segments` (sub-segments go to the cell after the one
containing their midpoint, and only ray points with positive x and y are used);
the units of the second branch of `Murat_traveling`; the L-curve without the
residual outside the range of the matrix (MATLAB passes the full `U` of the
SVD); `atan(x/y)` instead of `atan2` when converting grid points to lon/lat; the
declustering grid (latitude and longitude steps follow the original code); only
the first column of the first Paasschens sample is corrected; the Qc
uncertainty is stored as `1/|R|` and the "uncertainties below threshold" message
of `Murat_dataWarning` tests it against `fitTresholdLinear` as MATLAB does.

## Validation

MuRAT.jl was compared with MATLAB MuRAT 4.0 (R2025b) on the three test data sets
distributed with MuRAT (`sac_MSH`, `sac_Toba`, `sac_Romania`, with the input
files in `examples/`), using `inversionMethod = 'Tikhonov'` in both:

| | Mount St. Helens | Toba | Romania |
|---|---|---|---|
| Data | 448 traces, 4 frequencies | 863 recordings (3 components), 2 frequencies | 1059 traces, 4 frequencies |
| Velocity model | 3D (`modvMSH.txt`) | 1D (`iasp91.txt`) | 1D LITHO1.0 |
| Qc measurement | NonLinear | NonLinear, `fitTresholdLinear = 0` | Linearized |
| Velocity grid | median rel. diff. 1e-12 (max 6e-3, triangulation) | ≤ 2e-13 | ≤ 7e-8 (single precision) |
| Peak delays, travel times | identical / ≤ 1e-11 | identical | identical |
| Qc⁻¹, Qc kernels, peak-delay matrix | ≤ 3e-5 | ≤ 1e-6 | ≤ 5e-6 |
| Energy ratios | median 1e-8 | median 3e-5 | median 2e-7 |
| Selected data (all masks) | identical | identical | identical |
| Final maps, correlation (PD / Qc / Q) | 0.988-0.9997 / 1.0000 / ≥ 0.997 | 1.0000 / 1.0000 / 1.0000 | 1.0000 / 1.0000 / 1.0000 |
| Misfit reductions | same | same | same |

The data stage is compared on MATLAB's velocity models (only the data loop
differs), the final maps on Julia's own models: the small map differences for
Mount St. Helens come from the 3D velocity model, whose Delaunay triangulation of
a regular grid is not unique (MATLAB and Qhull choose differently).

Time (Apple M-series): Mount St. Helens MATLAB 75 s (serial), Julia 15-25 s
with 6 threads. On a data set of about 180,000 traces and 2 frequencies, a
complete Julia run with figures takes about 14 minutes with 6 threads.

## MATLAB → Julia functions

| MATLAB | Julia (file) |
|--------|--------------|
| `Murat_checks`, `Murat_modv1D`, `Murat_modv3D`, `Murat_rescale` (`fillNaN3`), `Murat_readLithos1` | `murat_checks`, `murat_modv1d`, `murat_modv3d`, `rescale_field` (`fill_nan3`), `read_litho1` (checks.jl, interpolation.jl, ext/MuratNCExt.jl) |
| `fget_sac`, `Murat_testData` | `read_sac`, `test_data` (io_sac.jl) |
| `Murat_dataParallelized`, `Murat_data` | `murat_data` (data.jl) |
| `Murat_envelope`, `Murat_location`, `Murat_picking`, `Murat_originTime`, `Murat_peakDelay(Check)`, `Murat_codaCheck`, `Murat_Qc`, `Murat_body` | `murat_envelope`, `murat_location`, `murat_picking`, `murat_origin_time`, `murat_peak_delay(_check)`, `murat_coda_check`, `murat_qc`, `murat_body` (measurements.jl) |
| `Murat_rays`, `Murat_tracing`, `Murat_bending`, `Murat_traveling`, `Murat_velocity`, `Murat_veldicity`, `Murat_cornering`, `Murat_segments` | same names in snake_case (rays.jl) |
| `Murat_kernels`, `Murat_paasschensFunction`, `Murat_codaMatrix` (`fixKernelLocal`) | `murat_kernels`, `paasschens`, `convolution_catalog`, `murat_coda_matrix_differential` (`fix_kernel!`), `murat_coda_matrix` for `kernelType = "middle"` (kernels.jl) |
| `Murat_selection`, `Murat_dataWarning`, `Murat_components`, `Murat_retainPeakDelay`, `Murat_retainQc`, `Murat_declustering` | same names in snake_case (selection.jl) |
| `Murat_inversion`, `Murat_inversionQc`, `Murat_inversionQ`, `Murat_tikhonovQc/Q`, `Murat_tikhonov_nonneg_fista`, `Murat_lsqlinQmean`, `Murat_smoothing`, `Murat_inputTesting`, `l_curve_tikh_svd` | `murat_inversion`, `murat_inversion_qc`, `murat_inversion_q`, `tikhonov_solve`, `tikhonov_nonneg_fista`, `murat_lsqlin_qmean`, `murat_smoothing`, `murat_input_testing`, `l_curve_tikh_svd` (inversion.jl) |
| `Murat_plot` and `Murat_image*` | `murat_summary`, `murat_paraview` (output.jl), `murat_matfiles` (ext/MuratMATExt.jl), `murat_plot` (ext/MuratMakieExt.jl) |
| `Murat_fold`, `Murat_unfold`, `Murat_unfoldXYZ`, `Murat_DDcoordinates`, `Murat_checkerboard3D`, `envelope(...,'rms')`, `fillmissing(...,'nearest')`, `interp1(...,'pchip')`, `distance`, `reckon`, `fit(...,'Robust','on')` | `fold_field`, `unfold_mesh`, `unfold_xyz`, `dd_coordinates`, `checkerboard3d`, `rms_envelope`, `fill_nearest!`, `pchip`, `geodistance`, `georeckon`, `robust_line_fit` (utils.jl, geodesy.jl, interpolation.jl) |

## Credits and license

MuRAT.jl is a Julia implementation of MuRAT
(https://github.com/LucaDeSiena/MuRAT), Copyright (c) 2024 Luca De Siena and
Authors (see MuRAT's `AUTHORS.md`), released under the MIT License. The MuRAT
copyright and permission notice are included in [`LICENSE.md`](LICENSE.md), as
the MIT License requires for copies and derived works. MuRAT.jl itself is
released under the same MIT License, Copyright (c) 2026 Chiara Nardoni and
Luca De Siena.

The test data sets and velocity models in `examples/` are copied from the MuRAT
repository.

If you use MuRAT.jl, please cite MuRAT (see its `CITATION.cff` and Zenodo DOI)
together with this Julia implementation.

## How MuRAT.jl was developed

MuRAT.jl was ported from the MATLAB code of MuRAT 4.0 with the help of Claude
(Anthropic), an AI coding assistant, under the direction of the authors. The
authors reviewed the code, decided on the changes from MATLAB (documented in
*Differences from MATLAB MuRAT 4.0*) and validated the results against MATLAB on
the three MuRAT test data sets (see *Validation*). The authors are responsible
for the code.
