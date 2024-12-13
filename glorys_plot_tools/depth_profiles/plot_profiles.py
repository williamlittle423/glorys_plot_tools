import netCDF4 as nc 
from glorys_plot_tools import glorys_download
from glorys_plot_tools.depth_profiles.plot_depth_animation import plot_depth_animation
from glorys_plot_tools.depth_profiles.plot_static_profiles import plot_static_profile
from glorys_plot_tools.duct_calculations.calculate_duct_data import calculate_duct_properties
import glob

def plot_profiles(
    lat,
    lon,
    start_year,
    start_month,
    start_day,
    end_year,
    end_month,
    end_day,
    max_depth,
    animation=1,
    animation_plot_title=None,
    animation_save_path=None,
    include_duct_depth=1,
    static_plot_title=None,
    static_save_path=None,
    static_include_mean=1,
    static_include_daily=1,
    fps=2,
    static=1,
):
    """
    Plot the static and/or animated profiles of the temperature, salinity, and density at a given location and time.
    """

    # Download the dataset
    glorys_download.download_data(
        start_year=start_year,
        start_month=start_month,
        start_day=start_day,
        end_year=end_year,
        end_month=end_month,
        end_day=end_day,
        minimum_latitude=lat,
        maximum_latitude=lat,
        minimum_longitude=lon,
        maximum_longitude=lon,
        maximum_depth=max_depth
    )
    
    # Parameters validation
    if start_month < 1 or start_month > 12:
        raise ValueError(f"Invalid start month: {start_month}")
    if end_month < 1 or end_month > 12:
        raise ValueError(f"Invalid end month: {end_month}")
    if start_day < 1 or start_day > 31:
        raise ValueError(f"Invalid start day: {start_day}")
    if end_day < 1 or end_day > 31:
        raise ValueError(f"Invalid end day: {end_day}")
    if animation == 0 and static == 0:
        raise ValueError("At least one of 'animation' or 'static' must be set to 1.")

    # GLORYS date format: 1993-01-01T00:00:00",
    start_month = f'0{start_month}' if start_month < 10 else str(start_month)
    start_day = f'0{start_day}' if start_day < 10 else str(start_day)
    
    end_month = f'0{start_month}' if end_month < 10 else str(end_month)
    end_day = f'0{end_day}' if end_day < 10 else str(end_day)

    # Find the dataset of the following format: 'cmems*.nc
    file_pattern = 'cmems*.nc'

    matching_files = glob.glob(file_pattern)

    if not matching_files:
        raise FileNotFoundError(f"Error finding downloaded file: {file_pattern}")

    file_path = matching_files[0]

    glorys_dataset = nc.Dataset(file_path, 'r')
    len_t = len(glorys_dataset['time'])
    glorys_dataset.close()

    stem = file_path.rsplit('.', 1)[0]
    parts = stem.rsplit('_', 4)
    suffix = '_'.join(parts[-4:])
    
    data = nc.Dataset(file_path, 'r')
    
    len_t = len(data['time'])
    
    data.close()
    
    if include_duct_depth == 1:
        print('Duct depth included. Performing duct analysis on data...')
        # Perform duct analysis on the data
        duct_fp = 'duct_data_' + suffix + '.nc'
        calculate_duct_properties(file_path, duct_fp, len_t)


    if animation == 1:
        # Create an animation of the profiles
        plot_depth_animation(file_path, duct_fp, suffix, animation_plot_title, animation_save_path, fps=fps, plot_duct_depth=include_duct_depth)
    
    if static == 1:
        # Create a static plot of the profiles
        plot_static_profile(file_path, duct_fp, static_save_path, static_plot_title, include_duct_depth, suffix, static_include_mean, static_include_daily) 

    return None
