import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import netCDF4 as nc
from matplotlib.animation import FuncAnimation, PillowWriter
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from datetime import datetime, timedelta
import cmocean

var_to_cbar_label = {
    'thetao': 'Temperature [C°]',
    'so': 'Salinity [PS]',
    'uo': 'Zonal Velocity [m/s]',
    'vo': 'Meridional Velocity [m/s]',
    'zos': 'Sea Surface Height [m]'
}

surface_var_to_label = {
    'thetao': 'SST',
    'so': 'Surface Salinity',
    'uo': 'Surface Zonal Velocity',
    'vo': 'Surface Meridional Velocity',
    'zos': 'Sea Surface Height'
}

non_surface_var_to_label = {
    'thetao': 'Temperature',
    'so': 'Salinity',
    'uo': 'Zonal Velocity',
    'vo': 'Meridional Velocity',
    'zos': 'Sea Surface Height'
}

cmap_string_to_map = {
    'thermal': cmocean.cm.thermal,
    'haline': cmocean.cm.haline,
    'speed': cmocean.cm.speed,
    'deep': cmocean.cm.deep,
    'dense': cmocean.cm.dense,
    'matter': cmocean.cm.matter,
    'turbid': cmocean.cm.turbid,
    'amp': cmocean.cm.amp,
    'balance': cmocean.cm.balance,
    'curl': cmocean.cm.curl,
    'delta': cmocean.cm.delta,
    'diff': cmocean.cm.diff,
    'tarn': cmocean.cm.tarn,
    'gray': cmocean.cm.gray,
    'ice': cmocean.cm.ice,
    'phase': cmocean.cm.phase,
    'topo': cmocean.cm.topo,
    'solar': cmocean.cm.solar,
    'oxy': cmocean.cm.oxy,
    'rain': cmocean.cm.rain,
    'jet': 'jet'
}

def increment_latter_date(input_str, days):
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

def plot_frame(day_idx, var, var_data, latitudes, longitudes, depth, depth_idx, suffix, cmocean_cmap, overlay_region=None):
    
    if overlay_region is not None:
        plot_overlay_region = True
        FVlonr, FVlonl, FVlatb, FVlatt = overlay_region
    else:
        plot_overlay_region = False
        
    plt.clf()  # Clear the previous frame
    ax = plt.axes(projection=ccrs.Mercator(min_latitude=latitudes.min(), max_latitude=latitudes.max()))  # Mercator projection for the map

    cmap = cmap_string_to_map[cmocean_cmap]
    cmap.set_bad(color='white')  # Set NaN values to appear as white - helpful for debugging

    img = ax.imshow(var_data[day_idx], cmap=cmap, transform=ccrs.PlateCarree(),
                    extent=[longitudes.min(), longitudes.max(), latitudes.min(), latitudes.max()],
                    origin='lower')
    
    ax.set_extent([longitudes.min(), longitudes.max(), latitudes.min(), latitudes.max()], crs=ccrs.PlateCarree())

    var_min = np.nanmin(var_data)
    var_max = np.nanmax(var_data)
    
    img.set_clim(var_min, var_max)

    # Create the colorbar without the label
    cbar = plt.colorbar(img, ax=ax, orientation='vertical')
    
    # Set the colorbar label with desired fontsize
    cbar.set_label(var_to_cbar_label[var], fontsize=18)
    
    # Set the tick labels fontsize
    cbar.ax.tick_params(labelsize=18)
    
    updated_date = increment_latter_date(suffix, days=day_idx)

    if depth_idx == 0:
        ax.set_title(f'GLORYS {surface_var_to_label[var]} {updated_date}', fontsize=20)
    else:
        ax.set_title(f'GLORYS {non_surface_var_to_label[var]} at depth {depth} {updated_date}', fontsize=20)
    
    ax.coastlines()
    ax.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='lightgrey')

    # Add the rectangle
    if plot_overlay_region:
        rect = mpatches.Rectangle((FVlonl, FVlatb), FVlonr - FVlonl, FVlatt - FVlatb,
                              linewidth=1, edgecolor='black', facecolor='none', transform=ccrs.PlateCarree())
        ax.add_patch(rect)

    gl = ax.gridlines(draw_labels=True, crs=ccrs.PlateCarree(), linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 18, 'color': 'black'}
    gl.ylabel_style = {'size': 18, 'color': 'black'}
    
    plt.tight_layout()

def plot_var_animation(glorys_fp, suffix, var, fig_width, fig_height, fps, dpi, depth, cmocean_cmap, overlay_region):
    glorys_data = nc.Dataset(glorys_fp)
    latitudes = glorys_data.variables['latitude'][:]
    longitudes = glorys_data.variables['longitude'][:]
    depths = glorys_data.variables['depth'][:]
    
    # Find the depth index closest to the desired depth
    depth_idx = np.abs(depths - depth).argmin()
    
    # Check if the variable has a depth dimension
    if len(glorys_data.variables[var].shape) == 3:
        var_data = glorys_data.variables[var][:, :, :]
    else:
        var_data = glorys_data.variables[var][:, depth_idx, :, :]
    
    # Compute the spatial ranges
    lat_range = latitudes.max() - latitudes.min()
    lon_range = longitudes.max() - longitudes.min()
    
    # Adjust figure size based on the aspect ratio of the data if dimensions are not provided
    if fig_width == None:
        aspect_ratio = lat_range / lon_range if lon_range != 0 else 1.0
        fig_width = 16
        fig_height = fig_width * aspect_ratio
    
    fig = plt.figure(figsize=(fig_width, fig_height))
    
    anim = FuncAnimation(
        fig, 
        lambda day_idx: plot_frame(day_idx, var, var_data, latitudes, longitudes, depth, depth_idx, suffix, cmocean_cmap, overlay_region), 
        frames=var_data.shape[0]
    )

    anim.save(f'GLORYS_{var}_animated_{suffix}.gif', writer=PillowWriter(fps=fps), dpi=dpi)
    plt.show()

    glorys_data.close()