from glorys_plot_tools.var_animations.plot_var_animation import plot_var_animation
from glorys_plot_tools import glorys_download
import glob

def plot_animation(
    var,
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
    plot_title=None,
    save_path=None,
    fps=2,
    dpi=1,
    depth_m=0,
    overlay_lonl=None,
    overlay_lonr=None,
    overlay_latb=None,
    overlay_latt=None,
    existing_fp=None,
    cmocean_cmap='thermal',
    fig_width=None,
    fig_height=None,
):
    
    # Parameters validation
    if start_month < 1 or start_month > 12:
        raise ValueError(f"Invalid start month: {start_month}")
    if end_month < 1 or end_month > 12:
        raise ValueError(f"Invalid end month: {end_month}")
    if start_day < 1 or start_day > 31:
        raise ValueError(f"Invalid start day: {start_day}")
    if end_day < 1 or end_day > 31:
        raise ValueError(f"Invalid end day: {end_day}")
    
    if depth_m == 0:
        min_depth = 0.49402499198913574
        max_depth = 0.49402499198913574
    
    if not existing_fp:
        # Download the dataset
        glorys_download.download_data(
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
            minimum_depth=depth_m,
            maximum_depth=depth_m
        )
        
        # Find the dataset of the following format: 'cmems*.nc
        file_pattern = 'cmems*.nc'

        matching_files = glob.glob(file_pattern)

        if not matching_files:
            raise FileNotFoundError(f"Error finding downloaded file: {file_pattern}")

        glorys_fp = matching_files[0]
    else:
        glorys_fp = existing_fp
    
    stem = glorys_fp.rsplit('.', 1)[0]
    parts = stem.rsplit('_', 4)
    suffix = '_'.join(parts[-4:])
    
    # Turn the overlay region into a tuple
    if overlay_lonl is not None and overlay_lonr is not None and overlay_latb is not None and overlay_latt is not None:
        overlay_region = (overlay_lonr, overlay_lonl, overlay_latb, overlay_latt)
    else:
        overlay_region = None
    
    plot_var_animation(glorys_fp, suffix, var, fig_width, fig_height, fps, dpi, depth_m, cmocean_cmap, overlay_region)
