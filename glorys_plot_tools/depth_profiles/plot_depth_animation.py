import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import netCDF4 as nc
from matplotlib import animation
from matplotlib.animation import FuncAnimation, PillowWriter
import cartopy.crs as ccrs
from glorys_plot_tools.density import fill_density
from glorys_plot_tools.soundspeed import sound_speed
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from datetime import datetime, timedelta

def increment_latter_date(input_str, days=1):
    """
    Increments the day of the latter date in the input string.

    Parameters:
    - input_str (str): The input string formatted as '62.00W_42.00N_0.49-266.04m_1996-06-06-1996-08-07'
    - days (int): Number of days to increment. Default is 1.

    Returns:
    - str: The updated string with the incremented date.
    """
    try:
        # Split the string by underscores
        parts = input_str.split('_')
        
        if len(parts) != 4:
            raise ValueError("Input string does not have the expected number of parts separated by underscores.")
        
        # The last part contains two dates separated by a hyphen
        date_part = parts[3]
        date_strings = date_part.split('-')
        
        if len(date_strings) != 6:
            raise ValueError("Date part does not have the expected format with two dates.")
        
        # Extract the two dates
        start_date_str = '-'.join(date_strings[:3])  # '1996-06-06'
        
        # Parse the start date
        start_date = datetime.strptime(start_date_str, "%Y-%m-%d")
        
        # Increment the day
        new_start_date = start_date + timedelta(days=days)
        
        # Format the new start date
        new_start_date_str = new_start_date.strftime("%Y-%m-%d")
        
        # Reconstruct the date part
        new_date_part = f"{new_start_date_str}"
        
        # Replace the old date part with the new one
        parts[3] = new_date_part
        
        # Reconstruct the full string
        new_input_str = '_'.join(parts)
        
        return new_input_str
    
    except Exception as e:
        print(f"Error: {e}")
        return None

def plot_depth_profiles_animated(filename, duct_filename, xmin_xmax, daily_lists, suffix, save_path, fps, plot_duct_depth, num_t):
    with nc.Dataset(filename, 'r') as ds:
        depths = ds.variables['depth'][0:30]
        temperature = np.mean(ds.variables['thetao'][:, 0:30, 0, 0], axis=0)
        salinity = np.mean(ds.variables['so'][:, 0:30, 0, 0], axis=0)
        density = np.mean(ds.variables['density'][:, 0:30, 0, 0], axis=0)
        sound_speed = np.mean(ds.variables['c'][:, 0:30, 0, 0], axis=0)

        # Duct data
        if plot_duct_depth:
            with nc.Dataset(duct_filename, 'r') as duct_data:
                # Duct depth from mean of duct depths that are not zero
                duct_depth = np.mean(duct_data['duct_depth'][:, 0, 0][duct_data['duct_depth'][:, 0, 0] != 0])

    # Create figure and axes with increased figure size
    fig, axes = plt.subplots(1, 4, figsize=(18, 13), sharey=True)

    fig.subplots_adjust(
        left=0.1,      # Left margin
        right=0.95,    # Right margin
        bottom=0.1,    # Bottom margin
        top=0.9,      # Top margin
        hspace=0.1,     # Horizontal space between subplots
        wspace=0.3      # Vertical space between subplots
    )

    def update(time_idx):
        for ax in axes:
            ax.clear()

        # Temperature profile
        axes[0].plot(daily_lists['thetao'][time_idx, :], depths, label='Temperature', color='grey', alpha=0.5)
        axes[0].plot(temperature, depths, label='Temperature', color='blue')
        if plot_duct_depth:
            axes[0].axhline(y=duct_depth, color='black', linestyle='--', label=f'Duct Depth ({duct_depth:.1f} m)')
        axes[0].set_xlabel('Temperature (°C)', fontsize=14, labelpad=10)
        axes[0].set_ylabel('Depth (m)', fontsize=14, labelpad=10)
        axes[0].invert_yaxis()
        axes[0].set_xlim([xmin_xmax['thetao'][0], xmin_xmax['thetao'][1]])
        axes[0].tick_params(axis='both', labelsize=14)

        # Salinity profile
        axes[1].plot(daily_lists['so'][time_idx, :], depths, color='grey', alpha=0.5)
        axes[1].plot(salinity, depths, label='Salinity', color='orange')
        if plot_duct_depth:
            axes[1].axhline(y=duct_depth, color='black', linestyle='--', label=f'Duct Depth ({duct_depth:.1f} m)')
        axes[1].set_xlabel('Salinity (PSU)', fontsize=14, labelpad=10)
        axes[1].set_xlim([xmin_xmax['so'][0], xmin_xmax['so'][1]])
        axes[1].tick_params(axis='x', labelsize=14)

        # Density profile
        axes[2].plot(daily_lists['density'][time_idx, :], depths, color='grey', alpha=0.5)
        axes[2].plot(density, depths, label='Density', color='green')
        if plot_duct_depth:
            axes[2].axhline(y=duct_depth, color='black', linestyle='--', label=f'Duct Depth ({duct_depth:.1f} m)')
        axes[2].set_xlabel('Density (kg/m³)', fontsize=14, labelpad=10)
        axes[2].set_xlim([xmin_xmax['density'][0], xmin_xmax['density'][1]])
        axes[2].tick_params(axis='x', labelsize=14)

        # Sound Speed profile
        axes[3].plot(daily_lists['c'][time_idx], depths, color='grey', alpha=0.5)
        axes[3].plot(sound_speed, depths, color='red')
        if plot_duct_depth:
            axes[3].axhline(y=duct_depth, color='black', linestyle='--', label=f'Mean Duct Depth ({duct_depth:.1f} m)')
        axes[3].set_xlabel('Sound Speed (m/s)', fontsize=14, labelpad=10)
        axes[3].set_xlim([xmin_xmax['c'][0], xmin_xmax['c'][1]])
        axes[3].tick_params(axis='x', labelsize=14)
        # Place the legend in the top right
        axes[3].legend(loc='upper right', fontsize=16)

        days = time_idx
        new_date_str = increment_latter_date(suffix, days=days)
        title = f'Glorys Depth Profiles {new_date_str}'
        fig.suptitle(title, fontsize=20, y=0.95)

    ani = animation.FuncAnimation(fig, update, frames=num_t, repeat=True)

    if save_path:
        ani.save(save_path, writer='pillow', fps=fps, dpi=200)
        print(f'Animation saved to {save_path} successfully!')

    plt.show()


def find_var_xmin_xmax(fp, var_name):
    """
    Finds the global vmin and vmax for a duct variable.

    Parameters:
        duct_fp (str): Path to the duct data file.
        var_name (str): Name of the duct variable.

    Returns:
        float: Minimum value for the color scaling.
        float: Maximum value for the color scaling.
    """
    with nc.Dataset(fp, 'r') as ds:
        # Extract all data for the variable across all days
        var_all_days = ds[var_name][:, 0:30, 0, 0]
        # Compute global min and max, excluding NaNs
        vmin = np.nanmin(var_all_days)
        vmax = np.nanmax(var_all_days)
    return vmin, vmax

def find_global_x_lim(glorys_fp):
    var_names = ['thetao', 'so', 'density', 'c']
    
    xmin_xmax = {}
    
    # loop through all variables
    for var in var_names:
        xmin = []
        xmax = []
        # loop through all duct files
        # find vmin and vmax for each variable
        xmin_, xmax_ = find_var_xmin_xmax(glorys_fp, var)
        xmin.append(xmin_)
        xmax.append(xmax_)
        print(f"Variable: {var} - Min: {min(xmin)}, Max: {max(xmax)}")
        
        xmin_xmax[var] = (min(xmin), max(xmax))
    
    return xmin_xmax

def plot_depth_animation(glorys_fp, duct_fp, suffix, title, save_path, fps, plot_duct_depth):
    
    fill_density(glorys_fp)
                    
    dataset = nc.Dataset(glorys_fp, 'r')
    
    latitudes = dataset['latitude'][:]
    longitudes = dataset['longitude'][:]
        
    if save_path == None:
        save_path = f'animated_profile_{suffix}.gif'
    
    if title == None:
        title = f'GLORYS Depth Profiles {suffix}'
            
    xmin_xmax = find_global_x_lim(glorys_fp)
    
    print('NUMBER OF TIME STEPS: ', len(dataset['time']))

    daily_temperature = dataset['thetao'][:, :, 0, 0]
    daily_salinity = dataset['so'][:, :, 0, 0]
    daily_density = dataset['density'][:, :, 0, 0]
    daily_sound_speed = []
    for i in range(len(daily_temperature)):
        daily_sound_speed.append([sound_speed(d, s, t) for d, s, t in zip(dataset['depth'][:], daily_salinity[i, :], daily_temperature[i, :])])
    daily_lists = {'thetao': daily_temperature, 'so': daily_salinity, 'density': daily_density, 'c': daily_sound_speed}
    dataset.close()
    plot_depth_profiles_animated(glorys_fp, duct_fp, xmin_xmax, daily_lists, suffix, save_path, fps, plot_duct_depth, len(daily_temperature))
