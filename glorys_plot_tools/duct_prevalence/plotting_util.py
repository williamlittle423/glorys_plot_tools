import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean
import matplotlib.patches as mpatches
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

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

def plot_glorys_variable(
    var,
    file_path,
    lat_min,
    lat_max,
    lon_min,
    lon_max,
    depth_m=0,
    title=None,
    save_path=None,
    cmap="thermal",
    depth_contour=1,
    ssf_contour=1,
    gulf_stream_contour=1,
):
    """
    Plot a GLORYS variable with the same formatting as the `plot_probability_with_contours` function.

    Parameters:
    - var (str): Variable to plot (e.g., 'thetao', 'so', etc.).
    - file_path (str): Path to the GLORYS NetCDF file.
    - lat_min, lat_max, lon_min, lon_max (float): Geographic boundaries for the plot.
    - depth_m (float): Depth in meters to plot. Defaults to 0.
    - title (str): Title for the plot.
    - save_path (str): Path to save the plot image. If None, displays the plot.
    - cmap (str): Colormap to use. Defaults to "thermal".
    - depth_contour, ssf_contour, gulf_stream_contour (int): Flags to include specific contours.
    """
    ds = nc.Dataset(file_path)
    latitudes = ds.variables['latitude'][:]
    longitudes = ds.variables['longitude'][:]
    depths = ds.variables['depth'][:]
    depth_idx = np.abs(depths - depth_m).argmin()
    data = ds.variables[var][0, depth_idx, :, :]

    # Subset lat/lon
    lat_mask = (latitudes >= lat_min) & (latitudes <= lat_max)
    lon_mask = (longitudes >= lon_min) & (longitudes <= lon_max)
    data = data[lat_mask, :][:, lon_mask]
    latitudes = latitudes[lat_mask]
    longitudes = longitudes[lon_mask]

    # Plot setup
    plt.figure(figsize=(15, 10))
    ax = plt.axes(projection=ccrs.Mercator(min_latitude=latitudes.min(), max_latitude=latitudes.max()))
    cmap = cmap_string_to_map.get(cmap, cmocean.cm.thermal)
    cmap.set_bad(color='white')

    # Main data visualization
    img = ax.imshow(data, cmap=cmap, transform=ccrs.PlateCarree(),
                    extent=[longitudes.min(), longitudes.max(), latitudes.min(), latitudes.max()],
                    origin='lower')
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

    # Add coastlines and land features
    ax.coastlines()
    ax.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='lightgrey')

    # Add optional contours
    if depth_contour:
        ax.contour(longitudes, latitudes, data, levels=[-200], colors='cyan', linewidths=0.5, transform=ccrs.PlateCarree())
    if gulf_stream_contour:
        ax.contour(longitudes, latitudes, data, levels=[15], colors='yellow', linewidths=0.5, transform=ccrs.PlateCarree())
    if ssf_contour:
        ax.plot(longitudes, np.mean(latitudes), color='red', linewidth=1, transform=ccrs.PlateCarree(), label='Shelf-Slope Front (SSF)')

    # Add gridlines
    gl = ax.gridlines(draw_labels=True, crs=ccrs.PlateCarree(), linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 12, 'color': 'black'}
    gl.ylabel_style = {'size': 12, 'color': 'black'}

    # Add colorbar and title
    title = title if title else f"{var} at Depth {depth_m}m"
    plt.title(title, fontsize=12)
    cbar = plt.colorbar(img, ax=ax, label=f"{var} [{ds.variables[var].units}]")
    cbar.ax.tick_params(labelsize=10)

    # Save or display plot
    if save_path:
        plt.savefig(save_path, dpi=200)
    else:
        plt.show()

    ds.close()
