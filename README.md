# Streaming Calculation of Local Seismic Attributes

MATLAB implementation of the methods developed for my M.Sc. thesis on the **streaming calculation of local seismic attributes**.

## Overview

Local seismic attributes provide time-varying characteristics of seismic signals, such as local frequency, through localized time-frequency analysis.

In this project, a streaming approach is developed for the calculation of local seismic attributes. Unlike conventional approaches that may require repeated or computationally expensive calculations, the streaming formulation updates the desired attribute as new data samples become available.

The repository contains the MATLAB implementation of the main algorithms, demonstration scripts, and supporting experiments used during the development and evaluation of the proposed approach.

## Main Methods

The repository includes implementations for:

* Instantaneous frequency
* Local frequency
* Streaming local frequency
* Local time-frequency mapping
* Streaming local time-frequency mapping
* Common-frequency local and streaming mappings

The main frequency estimation methods are implemented through:

```text
inst_freq.m
loc_freq.m
stream_lf.m
```

where:

* `inst_freq` calculates instantaneous frequency.
* `loc_freq` calculates local frequency using regularized local estimation.
* `stream_lf` calculates streaming local frequency.

## Repository Structure

```text
streaming-local-seismic-attributes/
│
├── README.md
├── LICENSE
├── .gitignore
│
├── src/
│   ├── inst_freq.m
│   ├── loc_freq.m
│   ├── stream_lf.m
│   ├── ltfm.m
│   ├── sltfm.m
│   ├── lcf_slices.m
│   └── scf_slices.m
│
├── examples/
│   ├── frequency_comparison_demo.m
│   └── streaming_map_demo.m
│
├── experiments/
│   ├── lambda_effect.m
│   └── regularization_window_effect.m
│
└── data/
```

### `src/`

Contains the core functions used by the examples and experiments.

### `examples/`

Contains demonstration scripts for the main methods.

#### `frequency_comparison_demo.m`

Calculates and compares instantaneous, local, and streaming local frequency for two synthetic seismic signals using:

```text
inst_freq
loc_freq
stream_lf
```

#### `streaming_map_demo.m`

Calculates and displays:

* Local time-frequency map
* Streaming local time-frequency map
* Short-Time Fourier Transform (STFT)

for a synthetic signal using:

```text
ltfm
sltfm
```

The STFT is calculated using MATLAB's built-in functionality.

### `experiments/`

Contains supporting experiments used to investigate the behavior of the methods.

#### `lambda_effect.m`

Investigates the effect of the regularization parameter `lambda` on streaming local frequency by calculating the attribute for four different values of `lambda`.

#### `regularization_window_effect.m`

Investigates the effect of the regularization window size (`win_size`) on local frequency by calculating the attribute for four different window sizes.

## Data Availability

The repository does not contain the proprietary or non-distributable field dataset used in the thesis.

Consequently, the demonstration script associated with the field-data experiment is not included in the public repository.

The functions required for common-frequency local and streaming mapping are provided independently:

```text
lcf_slices.m
scf_slices.m
```

These functions can be used with appropriately available seismic data.

Synthetic examples are included where possible to demonstrate the main algorithms without requiring restricted datasets.

## Requirements

* MATLAB R2022a
* Signal Processing Toolbox

The code was developed and tested in **MATLAB R2022a**. The `spectrogram` function used in the examples requires the **Signal Processing Toolbox**.

Compatibility with other MATLAB versions has not been systematically evaluated.

## Usage

Clone or download the repository and add the `src` directory to the MATLAB path.

For example:

```matlab
addpath('src')
```

The demonstration scripts can then be executed from the `examples` directory.

For example:

```matlab
frequency_comparison_demo
```

or:

```matlab
streaming_map_demo
```

Before running an individual script, make sure that its required input data and MATLAB dependencies are available.

## Citation

If you use this code in academic work, please cite the corresponding thesis or publication associated with this repository.

Citation information will be added here when the associated academic publication is available.

## Author

**Mohammad Hasan Safari Araghi**

M.Sc. in Geophysics (Seismic Exploration)

Institute of Geophysics, University of Tehran

GitHub: [mh-safari](https://github.com/mh-safari)

## License

The license information for this repository will be provided in the `LICENSE` file.
