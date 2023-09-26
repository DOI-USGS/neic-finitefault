# NEIC Finite Fault Process

The information in this document outlines the steps taken to create a finite fault model for the NEIC.

- [NEIC Finite Fault Process](#neic-finite-fault-process)
  - [1. Getting Data](#1-getting-data)
  - [2. Inversion](#2-inversion)
  - [3. Plot the results](#3-plot-the-results)
  - [4. Make Changes](#4-make-changes)
    - [Removing Bad Channels](#removing-bad-channels)
    - [Shift Timing for Autocorrelation](#shift-timing-for-autocorrelation)
    - [Update Plane Orientation](#update-plane-orientation)
  - [5. Officially Incorporate Changes](#5-officially-incorporate-changes)
    - [Update Input Files](#update-input-files)
    - [Recalculate Green's Functions](#recalculate-greens-functions)
  - [6. Rerun modelling](#6-rerun-modelling)
  - [7. Rerun plots](#7-rerun-plots)

## 1. Getting Data

Put all data into a folder that you plan to work within:

1. Put waveform and displacement data in a data folder
   1. Teleseismic and strong motion data can both be retrieved from [IRIS fdsnws](https://service.iris.edu/fdsnws/). The data acquisition script uses obspy to get this data: `wasp manage acquire --help`. Other data types (insar, cgps, gps, etc) are generally, manually acquired from a contributor.
2. Get the moment tensor information and save it in a file in the format `<eventid>_cmt_CMT`. QuakeML data with tensor information can be retrieved from a USGS event's origin product page. For example: https://earthquake.usgs.gov/earthquakes/eventpage/us20003k7a/origin/detail

## 2. Inversion

Execute the full process and inversion script: `wasp model run --help`. For example, running an event ("us20003k7a") with all data from a folder "/home/user/us20003k7a_product" might look something like `wasp model run /home/user/us20003k7a_product usgs_model -g /home/user/us20003k7a_product/us20003k7a_cmt_cmt -d /home/user/us20003k7a_product/data -t cgps -t gps -t insar -t strong -t surf -t body -ina /home/user/us20003k7a_product/insar_ascending.txt -v /home/user/us20003k7a_product/vel_model.txt`

## 3. Plot the results

For each nodal plane folder (NPX), run the plotting routine. For example: `wasp plot neic /home/user/us20003k7a_product/20150916225432/ffm.0/NP1 /home/user/us20003k7a_product/us20003k7a_cmt_cmt -a -t cgps -t gps -t insar -t strong -t surf -t body -d -p --pub -e us20003k7a`

## 4. Make Changes

Copy the nodal plane folder (NPX) to a new location (NP3) so that you preserve change history and begin running commands on the new folder.

### Removing Bad Channels

Bad channels can be removed or downweighted using the provided management script: `wasp manage modify-dicts --help`. For example to downweight a channels "HN1" and "HNZ" at station "ABC" **and** to delete the stations channel "HN2": `wasp manage modify-dicts /home/user/us20003k7a_product/20150916225432/ffm.0/NP3 downweight -t strong -sc "ABC:HN1,HNZ"` **and** `wasp manage modify-dicts /home/user/us20003k7a_product/20150916225432/ffm.0/NP3 delete -t strong -sc "ABC:HN2"`

### Shift Timing for Autocorrelation

To shift timing use the shift matching script: `wasp process shift-match --help`. For example `wasp process-data shift-match /home/user/us20003k7a_product/20150916225432/ffm.0/NP3 cgps /home/user/us20003k7a_product/us20003k7a_cmt_cmt`

### Update Plane Orientation

To update the planes orientation, subfault size, etc update the segments.json file.

## 5. Officially Incorporate Changes

### Update Input Files

In order for the changes from [4. Make Changes](#4-make-changes) to be accepted the input files script needs to be run: `wasp manage update-inputs --help`

### Recalculate Green's Functions

Rerun the greens function for all the data sets you updated: `wasp process greens /home/user/us20003k7a_product/20150916225432/ffm.0/NP3 /home/user/us20003k7a_product/us20003k7a_cmt_cmt -t cgps -t tele`

## 6. Rerun modelling

Rerun the modelling fortran script with the datatypes specified: `/home/user/neic-finitefault/fortran_code/bin_inversion_gfortran_f95/run_moddelling body surf strong cgps gps insar`

## 7. Rerun plots

Rerun writing plots with the same commands from [3. Plot the results](#3-plot-the-results).
