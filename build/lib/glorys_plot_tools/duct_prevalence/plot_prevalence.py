from glorys_plot_tools import glorys_download
from glorys_plot_tools.duct_calculations.calculate_duct_data import calculate_duct_properties
from glorys_plot_tools.duct_prevalence.plotting_util import plot_probability_with_contours
from glorys_plot_tools.duct_prevalence.calculate_prevalence import create_dataset, calculate_total_prevalence
import glob
import netCDF4 as nc


def plot_prevalence(
    start_year: int,
    start_month: int,
    start_day: int,
    end_year: int,
    end_month: int,
    end_day: int,
    lat_min: float,
    lat_max: float,
    lon_min: float,
    lon_max: float,
    plot_title = "",
    gulf_stream_contour=1,
    shelf_slope_contour=1,
    depth_contour=1,
    exising_fp=None,
    dataset_id="cmems_mod_glo_phy_my_0.083deg_P1D-m",
    dataset_version="202311"
    ):

    if exising_fp is not None:
        file_path = exising_fp
    else:
        # Download the dataset
        glorys_download.download_data(
            dataset_id=dataset_id,
            dataset_version=dataset_version,
            start_year=start_year,
            start_month=start_month,
            start_day=start_day,
            end_year=end_year,
            end_month=end_month,
            end_day=end_day,
            minimum_latitude=lat_min,
            maximum_latitude=lat_max,
            minimum_longitude=lon_min,
            maximum_longitude=lon_max,
            maximum_depth=400.0
        )

        # GLORYS date format: 1993-01-01T00:00:00",
        start_month = f'0{start_month}' if start_month < 10 else str(start_month)
        start_day = f'0{start_day}' if start_day < 10 else str(start_day)

        end_month = f'0{end_month}' if end_month < 10 else str(end_month)
        end_day = f'0{end_day}' if end_day < 10 else str(end_day)

        # Find the dataset of the following format: 'cmems*.nc
        file_pattern = 'cmems*.nc'

        matching_files = glob.glob(file_pattern)

        if not matching_files:
            raise FileNotFoundError(f"Error finding downloaded file: {file_pattern}")

        file_path = matching_files[0]

    stem = file_path.rsplit('.', 1)[0]
    parts = stem.rsplit('_', 4)
    suffix = '_'.join(parts[-4:])

    duct_fp = 'duct_data_' + suffix + '.nc'

    # Calculate duct data using multiprocessing
    calculate_duct_properties(file_path, duct_fp)

    # Calculate total prevalence
    calculate_total_prevalence(duct_fp, file_path, suffix)

    prevalence_fp = f"duct-prevalence_{suffix}.nc"

    # Plot the prevalence
    plot_probability_with_contours(file_path, prevalence_fp, suffix, depth_contour, shelf_slope_contour, gulf_stream_contour, plot_title)
