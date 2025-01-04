import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import netCDF4 as nc
from matplotlib.animation import FuncAnimation, PillowWriter
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
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

def increment_latter_date(day_idx, start_date):
    """
    Increment the given date by a specified number of days.

    Parameters:
    day_idx (int): Number of days to increment (can be negative for decrement).
    start_date (str): Starting date in 'YYYY-MM-DD' format.

    Returns:
    str: Resulting date in 'YYYY-MM-DD' format.
    """
    try:
        # Parse the input date
        date_obj = datetime.strptime(start_date, '%Y-%m-%d')
        
        # Increment the date
        new_date = date_obj + timedelta(days=day_idx)
        
        # Return the new date in the same format
        return new_date.strftime('%Y-%m-%d')
    except ValueError as e:
        # Handle invalid date format
        return f"Invalid date format: {e}"


def plot_var_animation(glorys_fp, suffix, var, fig_width, fig_height, fps, dpi, depth, cmocean_cmap, overlay_region, start_date):
    """
    Creates an animation of the specified variable in the GLORYS dataset, 
    with a fixed colorbar (and colorbar size) throughout all frames.
    """

    # --- 1) Load the data ---
    glorys_data = nc.Dataset(glorys_fp)
    latitudes = glorys_data.variables['latitude'][:]
    longitudes = glorys_data.variables['longitude'][:]
    depths = glorys_data.variables['depth'][:]
    
    # Find the depth index closest to the desired depth
    depth_idx = np.abs(depths - depth).argmin()
    
    # Check if the variable has a depth dimension
    if len(glorys_data.variables[var].shape) == 3:
        # shape: (time, lat, lon)
        var_data = glorys_data.variables[var][:, :, :]
    else:
        # shape: (time, depth, lat, lon)
        var_data = glorys_data.variables[var][:, depth_idx, :, :]
    
    # Compute the global min and max for the entire var_data
    var_min = np.nanmin(var_data)
    var_max = np.nanmax(var_data)

    # --- 2) Figure setup & colorbar creation ---
    # Determine figure size if not provided
    lat_range = latitudes.max() - latitudes.min()
    lon_range = longitudes.max() - longitudes.min()
    if fig_width is None:
        aspect_ratio = lat_range / lon_range if lon_range != 0 else 1.0
        fig_width = 16
        fig_height = fig_width * aspect_ratio
    
    fig = plt.figure(figsize=(fig_width, fig_height))
    ax = plt.axes(projection=ccrs.Mercator(min_latitude=latitudes.min(), max_latitude=latitudes.max()))
    ax.set_extent([longitudes.min(), longitudes.max(), latitudes.min(), latitudes.max()], 
                  crs=ccrs.PlateCarree())
    
    # Prepare colormap
    cmap = cmap_string_to_map[cmocean_cmap]
    
    # Create a single imshow object, initially showing day 0
    # We set vmin and vmax here so the color scale is fixed
    img = ax.imshow(
        var_data[0],
        transform=ccrs.PlateCarree(),
        extent=[longitudes.min(), longitudes.max(), latitudes.min(), latitudes.max()],
        origin='lower',
        cmap=cmap,
        vmin=var_min,
        vmax=var_max
    )

    # Create colorbar only once
    cbar = plt.colorbar(img, ax=ax, orientation='vertical')
    cbar.set_label(var_to_cbar_label[var], fontsize=18)
    cbar.ax.tick_params(labelsize=18)

    # (Optional) Add coastlines, land, etc. (only once)
    ax.coastlines()
    ax.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='lightgrey')

    # Add gridlines, label styles
    gl = ax.gridlines(draw_labels=True, crs=ccrs.PlateCarree(),
                      linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 18, 'color': 'black'}
    gl.ylabel_style = {'size': 18, 'color': 'black'}

    # If overlay region is provided, create the rectangle once
    if overlay_region is not None:
        FVlonr, FVlonl, FVlatb, FVlatt = overlay_region
        rect = mpatches.Rectangle((FVlonl, FVlatb), FVlonr - FVlonl, FVlatt - FVlatb,
                                  linewidth=1, edgecolor='black', facecolor='none', 
                                  transform=ccrs.PlateCarree())
        ax.add_patch(rect)

    # --- 3) Animation update function ---
    def update_frame(day_idx):
        """Updates the imshow data and the title for each frame."""
        # Update the data on the existing 'img'
        img.set_data(var_data[day_idx])
        
        # The color limits remain unchanged (vmin=var_min, vmax=var_max)
        # so we do NOT call set_clim each time unless you want to reassert them:
        # img.set_clim(var_min, var_max)
        
        # Update the title
        updated_date = increment_latter_date(day_idx, start_date)
        if depth_idx == 0:
            ax.set_title(f'GLORYS {surface_var_to_label[var]} {updated_date}', fontsize=20)
        else:
            ax.set_title(f'GLORYS {non_surface_var_to_label[var]} at depth {depth} {updated_date}', fontsize=20)
        
        # Return the artist we updated
        return [img]

    # --- 4) Create the animation ---
    anim = FuncAnimation(
        fig,
        update_frame,
        frames=var_data.shape[0],  # number of time steps
        interval=1000/fps,        # in milliseconds
        blit=False                # or True if you only return updated artists
    )

    # Save the animation
    anim.save(f'GLORYS_{var}_animated_{suffix}.gif', writer=PillowWriter(fps=fps), dpi=dpi)
    plt.show()

    # Close the dataset
    glorys_data.close()
