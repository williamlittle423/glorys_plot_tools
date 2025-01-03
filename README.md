# GLORYS Plot Tools

![Version](https://img.shields.io/badge/version-1.0-blue.svg)

## Overview

GLORYS Plot Tools is a Python 3 package for generating plots and performing duct analysis on Copernicus Marine products, including Global Ocean Physics Reanalysis and Analysis and Forecast datasets. Key features include:

- **Depth Profiles**: Create static or animated profiles for temperature, salinity, density, and soundspeed.
- **Duct Calculations**: Compute subsurface duct properties within specified regions and timeframes.
- **Duct Prevalence**: Visualize the prevalence of ducts over selected areas and periods.
- **Variable Animations**: Generate animations of various GLORYS variables.
- **Plot Variables**: Plot specific GLORYS variables at single time points or averaged over periods.

## Prerequisites

Ensure the following dependencies are installed:

- **Python 3.10+**
- **HDF5**
- **Rust**

### Installation Methods

- **Python**:
  - **Conda**:
    ```bash
    conda install python=3.10
    ```
  - **Ubuntu**:
    ```bash
    sudo apt-get update
    sudo apt-get install python3.10 python3-pip
    ```
  - **Official Installer**: [Download Python](https://www.python.org/downloads/)
  
- **HDF5**:
  - **Conda**:
    ```bash
    conda install hdf5
    ```
  - **Ubuntu**:
    ```bash
    sudo apt-get install libhdf5-dev
    ```
  
- **Rust**: [Install Rust](https://www.rust-lang.org/tools/install)

*Ensure all dependencies are added to your system’s PATH.*

## Installation

1. Download the `.whl` file from the `dist/` directory of the GLORYS Plot Tools project.
2. Navigate to the directory containing the `.whl` file:
    ```bash
    cd path/to/directory
    ```
3. Install the package using pip:
    ```bash
    pip install glorys_plot_tools-1.0.0-py3-none-any.whl
    ```

## Usage

Import the desired functions from the package and execute them with the necessary parameters.

### Example: Plot Depth Profiles

```python
from glorys_plot_tools.depth_profiles.plot_profiles import plot_profiles

plot_profiles(
    lat=42.6,
    lon=-60.2,
    start_year=2020,
    start_month=6,
    start_day=26,
    end_year=2020,
    end_month=8,
    end_day=26,
    max_depth=0.5,
    animation=1,
    include_duct_depth=1,
    fps=2,
    static=1,
)
```