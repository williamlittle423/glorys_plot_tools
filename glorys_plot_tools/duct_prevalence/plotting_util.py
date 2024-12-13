import netCDF4 as nc
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pkg_resources
import xarray as xr
import matplotlib
import cmocean
import matplotlib.patches as mpatches
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

# Enum to convert month index to month name
month_enum = {
    0: 'January',
    1: 'February',
    2: 'March',
    3: 'April',
    4: 'May',
    5: 'June',
    6: 'July',
    7: 'August',
    8: 'September',
    9: 'October',
    10: 'November',
    11: 'December',
}

def plot_probability_with_contours(glorys_fp, prevalence_fp, suffix, depth_contour, ssf_contour, gulf_stream_contour):

    glorys_dataset = nc.Dataset(glorys_fp, 'r')

    save_path = F'duct_prevalence_{suffix}.png'

    data = nc.Dataset(prevalence_fp, 'r')
    data = np.array(data['prevalence'][:])

    # Load the .mat file
    mat_file_path = pkg_resources.resource_filename('glorys_plot_utils', 'contours/SSFANNUAL_1973_2017.mat')
    mat_file = scipy.io.loadmat(mat_file_path)

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
    
    # Handle masked and NaN values
    if np.ma.is_masked(data):
        data = np.ma.filled(data, fill_value=np.nan)
    data = np.nan_to_num(data, nan=0.0)  # Replace NaN with 0 or another appropriate fill value
    
    gebco_file_path = pkg_resources.resource_filename('glorys_plot_utils', 'contours/contour_200m_depth.nc')
    gebco_dataset = xr.open_dataset(gebco_file_path)  # gebco elevation data
    
    temp_file_path = pkg_resources.resource_filename('glorys_plot_utils', 'contours/temperature_data.nc')
    temp_dataset = xr.open_dataset(temp_file_path)  # temperature data for gulf stream contour

    temps = temp_dataset['temperature'].isel(time=9)

    # Temperature data for contour (GLORYS data)
    latitudes = np.array(glorys_dataset['latitude'])
    longitudes = np.array(glorys_dataset['longitude'])
    lat_grid, lon_grid = np.meshgrid(latitudes, longitudes, indexing='ij')

    # Gebco data
    gebco_latitudes = gebco_dataset['lat']
    gebco_longitudes = gebco_dataset['lon']
    gebco_altitudes = gebco_dataset['elevation'].values
    gebco_lat_grid, gebco_lon_grid = np.meshgrid(gebco_latitudes, gebco_longitudes, indexing='ij')

    plt.figure(figsize=(15, 10))
    ax = plt.axes(projection=ccrs.Mercator(min_latitude=latitudes.min(), max_latitude=latitudes.max()))  # Use Mercator projection for the map
    cmap = cmocean.cm.haline
    cmap.set_bad(color='white')  # Set NaN values to appear as white

    print('Shape of data: ', data.shape)

    img = ax.imshow(data, cmap=cmap, transform=ccrs.PlateCarree(),
                    extent=[longitudes.min(), longitudes.max(), latitudes.min(), latitudes.max()],
                    origin='lower')
    
    ax.set_extent([longitudes.min(), longitudes.max(), latitudes.min(), latitudes.max()], crs=ccrs.PlateCarree())

    # Add coastlines and land features
    ax.coastlines()
    ax.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='lightgrey')

    # Add contour for 200m depth from gebco data
    if depth_contour == 1:
        matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
        gebco_depth_contour = ax.contour(gebco_lon_grid, gebco_lat_grid, gebco_altitudes, levels=[-200], colors='cyan', linewidths=0.5, transform=ccrs.PlateCarree())
        ax.clabel(gebco_depth_contour, inline=1, fontsize=0)

    # Add gulf stream contour
    if gulf_stream_contour == 1:
        glorys_temp_contour = ax.contour(lon_grid, lat_grid, temps, levels=[15], colors='yellow', linewidths=0.5, transform=ccrs.PlateCarree())
        ax.clabel(glorys_temp_contour, inline=1, fontsize=0)

    # Add shelf slope front contour
    if ssf_contour == 1:
        ax.plot(lon, mean_lat, color='red', linewidth=1, transform=ccrs.PlateCarree(), label='Shelf Slope Front (SSF)')

    # Legend handles
    contour_handles = [mpatches.Patch(color='cyan', label='200m Depth')] + \
                      [mpatches.Patch(color='yellow', label='Gulf Stream')] + \
                      [mpatches.Patch(color='red', label='Shelf-Slope Front (SSF)')]
    
    # Gridlines and labels
    gl = ax.gridlines(draw_labels=True, crs=ccrs.PlateCarree(), linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 12, 'color': 'black'}
    gl.ylabel_style = {'size': 12, 'color': 'black'}

    # Color bar and title
    plt.title(f'Subsurface Duct Prevalence (%) {suffix}', fontsize=12)

    cbar = plt.colorbar(img, ax=ax, label='Probability (%)', orientation='vertical')
    cbar.ax.tick_params(labelsize=10)
    plt.legend(handles=contour_handles, loc='upper right')

    if save_path:
        plt.savefig(save_path, dpi=200)
    plt.show()

if __name__ == '__main__':
    plot_probability_with_contours('cmems_mod_glo_phy_my_0.083deg_P1D-m_multi-vars_80.00W-42.00W_30.00N-53.00N_0.49m_1993-06-01-1993-06-05.nc',
                                    'duct-prevalence_80.00W-42.00W_30.00N-53.00N_0.49m_1993-06-01-1993-06-05.nc',
                                    '1993-06-01-1993-06-05', 1, 1, 1)