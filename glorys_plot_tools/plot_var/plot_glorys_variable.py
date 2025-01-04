import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import glob
import cartopy.feature as cfeature
import glorys_plot_tools.glorys_download as glorys_download
from datetime import datetime, timedelta
import scipy.io
import pkg_resources
import cmocean
import matplotlib
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

def plot_glorys_variable(
    var,
    start_year,
    start_month,
    start_day,
    lat_min,
    lat_max,
    lon_min,
    lon_max,
    end_year=None,
    end_month=None,
    end_day=None,
    depth_m=0,
    title=None,
    existing_fp=None,
    dataset_id="cmems_mod_glo_phy_my_0.083deg_P1D-m",
    dataset_version="202311",
    depth_contour=1,
    ssf_contour=1,
    gulf_stream_contour=1,
    depth_contour_color='cyan',
    ssf_contour_color='red',
    gulf_stream_contour_color='yellow',
    save_path=None,
    cmap="jet",
    show_plot=True
):
    """
    Plot a GLORYS variable at a specific time point or averaged over a specified period,
    with formatting similar to the provided "plot_probability_with_contours" script.
    """

    if existing_fp is not None:
        glorys_fp = existing_fp
    else:
        # Download the dataset
        if end_year is None:
            end_year = start_year
            end_month = start_month
            end_day = start_day

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
            minimum_depth=depth_m,
            maximum_depth=depth_m,
        )

        # Find the dataset file
        file_pattern = "cmems*.nc"
        matching_files = glob.glob(file_pattern)

        if not matching_files:
            raise FileNotFoundError(f"Error finding downloaded file: {file_pattern}")

        glorys_fp = matching_files[0]

    # Extract date and location suffix for naming duct file
    stem = glorys_fp.rsplit('.', 1)[0]
    parts = stem.rsplit('_', 4)
    suffix = '_'.join(parts[-4:])

    # Open the dataset
    ds = nc.Dataset(glorys_fp)
    latitudes = ds.variables["latitude"][:]
    longitudes = ds.variables["longitude"][:]
    depths = ds.variables["depth"][:]
    depth_idx = np.abs(depths - depth_m).argmin()

    """
    Extract SSF data from .mat file
    """
    ssf_file_path = pkg_resources.resource_filename('glorys_plot_tools', 'contours/SSFANNUAL_1973_2017.mat')

    # Load the .mat file
    mat_file = scipy.io.loadmat(ssf_file_path)

    # Extract the SSFANNUAL structure
    SSFANNUAL = mat_file['SSFANNUAL']

    # Extract the fields
    time = SSFANNUAL['time'][0, 0].flatten()
    lat = SSFANNUAL['lat'][0, 0]
    lon = SSFANNUAL['lon'][0, 0].flatten()

    # Remove NaN values from lat and corresponding lon
    # This is required for the plot to work
    valid_idx = ~np.isnan(lat).any(axis=1)
    lat = lat[valid_idx, :]
    time = time[valid_idx]

    start_idx = 249
    mean_lat = np.nanmean(lat[start_idx::12], axis=0)

    """
    Extract Depth Contour data from netCDF4 file
    """
    depth_contour_fp = pkg_resources.resource_filename('glorys_plot_tools', 'contours/contour_200m_depth.nc')

    gebco_dataset = nc.Dataset(depth_contour_fp)

    # Gebco data
    gebco_latitudes = gebco_dataset['lat']
    gebco_longitudes = gebco_dataset['lon']
    gebco_altitudes = gebco_dataset['elevation'][:]
    gebco_lat_grid, gebco_lon_grid = np.meshgrid(gebco_latitudes, gebco_longitudes, indexing='ij')

    # Extract data and take the mean over time if an end date was specified
    if end_year is not None:
        if ds.variables[var].ndim == 4:
            data = ds.variables[var][:, depth_idx, :, :]
        else:
            data = ds.variables[var][:, :, :]
        if gulf_stream_contour:
            # check if temps is in ds.variables
            if 'thetao' not in ds.variables:
                raise ValueError("Gulf stream contour requires 'thetao' variable in dataset. Unable to perform this with duct datasets.")
            temps = ds.variables['thetao'][:]
            temps = np.mean(temps, axis=0)
        data = np.mean(data, axis=0)
    else:
        if ds.variables[var].ndim == 4:
            data = ds.variables[var][0, depth_idx, :, :]
        else:
            data = ds.variables[var][0, :, :]
        if gulf_stream_contour:
            if 'thetao' not in ds.variables:
                raise ValueError("Gulf stream contour requires 'thetao' variable in dataset. Unable to perform this with duct datasets.")
            temps = ds.variables['thetao'][0, :, :, :]

    if gulf_stream_contour:
        lat_grid, lon_grid = np.meshgrid(latitudes, longitudes, indexing='ij')

    # Subset lat/lon
    lat_mask = (latitudes >= lat_min) & (latitudes <= lat_max)
    lon_mask = (longitudes >= lon_min) & (longitudes <= lon_max)
    data = data[lat_mask, :][:, lon_mask]
    latitudes = latitudes[lat_mask]
    longitudes = longitudes[lon_mask]

    # Convert cmap string to a colormap object and set bad values as white
    cmap_obj = cmap_string_to_map.get(cmap, cmocean.cm.thermal)
    if isinstance(cmap_obj, str):
        cmap_obj = plt.get_cmap(cmap_obj)
    cmap_obj.set_bad(color='white')

    # Prepare the figure and axis with the desired formatting
    plt.figure(figsize=(15, 10))
    # Use a Mercator projection to be consistent with the styling script
    ax = plt.axes(projection=ccrs.Mercator(min_latitude=lat_min, max_latitude=lat_max))
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

    # Plot the variable data
    img = ax.pcolormesh(longitudes, latitudes, data, cmap=cmap_obj, transform=ccrs.PlateCarree())

    # Add coastlines, land, and gridlines
    ax.coastlines()
    ax.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='lightgrey')

    gl = ax.gridlines(draw_labels=True, crs=ccrs.PlateCarree(), linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 12, 'color': 'black'}
    gl.ylabel_style = {'size': 12, 'color': 'black'}

    contour_handles = []
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'

    # Reshape temps from (1, N_lat, N_lon) to (N_lat, N_lon)
    var
    if gulf_stream_contour:
        temps = temps[0, :, :]

    # Optionally plot contours
    if depth_contour:
        gebco_depth_contour = ax.contour(gebco_lon_grid, gebco_lat_grid, gebco_altitudes, levels=[-200], colors=depth_contour_color, linewidths=0.5, transform=ccrs.PlateCarree())
        contour_handles.append(mpatches.Patch(color=depth_contour_color, label='200m Depth'))
    if gulf_stream_contour:
        glorys_temp_contour = ax.contour(lon_grid, lat_grid, temps, levels=[15], colors=gulf_stream_contour_color, linewidths=0.5, transform=ccrs.PlateCarree())
        contour_handles.append(mpatches.Patch(color=gulf_stream_contour_color, label='Gulf Stream'))
    if ssf_contour:
        ax.plot(lon, mean_lat, color=ssf_contour_color, linewidth=1, transform=ccrs.PlateCarree(), label='Shelf Slope Front (SSF)')
        contour_handles.append(mpatches.Patch(color=ssf_contour_color, label='Shelf-Slope Front (SSF)'))
        
    if title is not None:
        ax.set_title(title, fontsize=14)
    else:
        if depth_m < 0.5:
            ax.set_title(f'{surface_var_to_label[var]} {suffix}', fontsize=14)
        else:
            ax.set_title(f'{non_surface_var_to_label[var]} {suffix}', fontsize=14)

    cbar = plt.colorbar(img, ax=ax, orientation="vertical", pad=0.05)
    cbar.ax.tick_params(labelsize=10)
    cbar.set_label(f"{var_to_cbar_label[var]}", fontsize=12)

    # Add legend if contours are present
    if contour_handles:
        plt.legend(handles=contour_handles, loc='upper right')

    if save_path == None:
        save_path = f'{var}_{suffix}.png'
    
    plt.savefig(save_path, dpi=200)
    if show_plot:
        plt.show()

    ds.close()
