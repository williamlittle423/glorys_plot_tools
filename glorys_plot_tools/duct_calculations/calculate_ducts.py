from glorys_plot_tools.duct_calculations.calculate_duct_data import calculate_duct_properties
from glorys_plot_tools import glorys_download
import glob

def calculate_ducts(
    lat_min,
    lat_max,
    lon_min,
    lon_max,
    start_year,
    start_month,
    start_day,
    end_year,
    end_month,
    end_day,
    dataset_id="cmems_mod_glo_phy_my_0.083deg_P1D-m",
    dataset_version="202311",
    existing_fp=None,
    ducts_save_fp=None
    ):

    if existing_fp is not None:
        glorys_fp = existing_fp
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
            maximum_depth=350
        )

        # Find the dataset of the following format: 'cmems*.nc
        file_pattern = 'cmems*.nc'

        matching_files = glob.glob(file_pattern)

        if not matching_files:
            raise FileNotFoundError(f"Error finding downloaded file: {file_pattern}")

        glorys_fp = matching_files[0]

    # Extract date and location suffix for naming duct file
    stem = glorys_fp.rsplit('.', 1)[0]
    parts = stem.rsplit('_', 4)
    suffix = '_'.join(parts[-4:])

    duct_fp = 'duct_data_' + suffix + '.nc'

    # Calculate ducts
    calculate_duct_properties(
        glorys_fp,
        ducts_save_fp=duct_fp
    )
