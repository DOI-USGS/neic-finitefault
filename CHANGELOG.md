# Change Log
- [Change Log](#change-log)
- [1.X.X](#1xx)
  - [Changed](#changed)
  - [Removed](#removed)
- [1.1.0](#110)
  - [Added](#added)
  - [Changed](#changed-1)
  - [Fixed](#fixed)
  - [Removed](#removed-1)
- [1.0.0 ](#100-)
  - [Added](#added-1)
  - [Changed](#changed-2)
  - [Fixed](#fixed-1)
  - [Removed](#removed-2)
- [0.1.0 (Provisional Release)](#010-provisional-release)
  - [Added](#added-2)
  - [Changed](#changed-3)

# 1.X.X

## Changed
- Simplified Dockerfile for ease of use

## Removed
- Removed automated publishing of Docker images from pipeline

# [1.1.0](https://code.usgs.gov/ghsc/neic/algorithms/neic-finitefault/-/releases/1.1.0)

## Added 
- Option for static offset ramp added for imagery data
- New `ShakeRupture` class added to write polygons for ShakeMap. Includes polygons for multiple segments.
- Waveform plots now include units
- Added option to reformat synthetic output as SAC files (#175)
- New `misfit_details.txt` included in output
- Added `velocity_model` option to `wasp manage update-inputs`
  

## Changed
- InSAR data is all grouped into a general imagery category rather than divided into ascending/descending.
- Range of allowed values for imagery data ramp coefficients increased
- References to "GPS" updated to "GNSS"
- CLI now `ffm` instead of `wasp`
- `remove-response` command now `skip-remove-response`
- `downloads` command now `generate-downloads`
- Updated Spanish text to English: "Solucion" now "Solution"
- `wasp process greens` now recalculates local GF bank (#185)
- cGPS filename convention updated to `[LXY][ENZ].sac` instead of `[LXY]*.sac`
- Surface wave sorted into directories named `RAYLEIGH` and `LOVE` rather than all in `LONG`. `SH` updated to `BHT`
- Waveform file names updated from `final_*` to `processed_*`
- Made rake vector scale dependent on fault segment width if using autoscale to ensure rake vectors scale correctly for multi-segment models where segments are of different widths
- Labels use "AsymmetricCosine" instead of "Asymetriccosine"


## Fixed
- Previously .fsp format did not accommodate multiple segments with different Dx/Dz values. New code adds Dx, Dz metadata into each segment header in case of different subfault sizes for multisegment models.
- Fixed moment rate units to be Nm/s
- Shift match plotting fixed for negative shifts (#195)
- Fixed issue with static grid size for Mw higher than 8.4 in `write_Coulomb_file`.
- Python text inconsistency on different systems due to long decimals
- vel_model.txt not overwritten during plotting
- Fixed accidental `final_final_` surface wave file naming
- `static_to_fsp` labels datetime as "UTC"
  
## Removed
- Scale removed from map plot


# [1.0.0 ](https://code.usgs.gov/ghsc/neic/algorithms/neic-finitefault/-/releases/1.0.0)
https://doi.org/10.5066/P1EKKUNW

## Added 
- Manual shift  match option added for `shift_match` for #156
- cutde added as dependency for plotting
- Checks and raised errors added for incorrectly entered stations/components arguments in `_modify_by_dict` 
- Jupyter notebook added as dependency
- Added option to add multiple ascending/descending InSAR files in `wasp manage fill-dicts`


## Changed
- Installation process replaced Miniconda with Miniforge
- `wasp manage update-inputs` `directory` argument now required
- Installation process updated to utilize conda rather than installing from source
- `wasp manage fill-dicts` InSAR options changed from `insar_ascending: pathlib.Path`, `insar_ascending_ramp: InsarRamp`, `insar_descending: pathlib.Path`, `insar_descending_ramp: InsarRamp` to `insar_ascending: List[str]` and `insar_descending: List[str]` in format `<path>:<ramp>`
- Defaults for `filling_data_dicts` method arguments (`insar_asc` and `insar_desc`) updated from `List[None]` to `None`

## Fixed
- Fixed description of map limits parameter, and changed KML plot so that map limits parameter is consistent with other plots.
- SRF format skipped if solution is static-only in `wasp plot neic`, since SRF format is incompatible with static solutions.
- Fixed strong motion and cgs labels for `execute_plot`
- `wasp model run-multiple` script executes plots in series instead of parallel to avoid issues with PyGMT
- Install script process for adding configuration line to the low.in file now checks if the line already exists before adding
- `wasp model run-multiple` `solutions_folders` option fixed to be type `List[pathlib.Path]`
- InSAR plot modified to better fit aspect ratio regions
- Observed cGPS waveforms are extended by a constant value (the static displacement). However, the cGPS Green's functions, whose size is always 1024, are truncated to zero about 80 points before the end of the waveform. In the wavelet domain, this is an issue because the FFT of the observed waveform for cGPS observations (not truncated to zero) then does not match the FFT of the cGPS synthetics (truncated to zero), leading to a mismatch of wavelet coefficients. The cGPS waveforms were extended to the static displacement, until reaching the point at 1024-80. From then on, the observed waveform is truncated to zero, for consistency with the Green's functions. We also update the cgps_waves.json wavelet_weight values to encourage better fitting of low frequencies.

## Removed
- Removed code related to the no longer used `Event_mult.in` file (#134)
- okada_wrapper dependency removed
- Removed unused files in `fortran_code/` directory

# [0.1.0 (Provisional Release)](https://code.usgs.gov/ghsc/neic/algorithms/neic-finitefault/-/releases/0.1.0)

## Added
- All code transferred from GitHub
- Added enhanced tests
- Added Typer backed command line interface

## Changed
- Command names have slight variations from GitHub code