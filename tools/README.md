!!! Info

    This file describes the contents of the `#!bash kite/tools`-folder. Note that the main post-processing tool is located in the `#!bash kite/build`-folder.

# `kite/tools`-folder
This folder contains bespoke tools:

## 1. `process_arpes.py`
to process the output of a spectral function calculation.

## 2. `process_single_shot.py`
to extract the output data from a single-shot calculation of the longitudinal DC-conductivity.
This tool is a simple python script that outputs a data file with three columns read out of
the processed HDF5, namely:

| Column 1       | Column 2                       | Column 3         |
|----------------|--------------------------------|------------------|
| Fermi Energies | `preserve_disorder` flag (0/1) | The conductivity |

KITE-tools extracts the same result:

``` bash
./build/KITE-tools output.h5 --SingleShot
```
