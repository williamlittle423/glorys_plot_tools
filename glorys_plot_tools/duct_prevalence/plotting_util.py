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

def plot_probability_with_contours(
    glorys_fp,
    prevalence_fp,
    suffix,
    depth_contour,
    ssf_contour,
    gulf_stream_contour,
    plot_title
):
    """
    Plots duct prevalence (%) on the GLORYS grid, optionally adding:
    1) A 200m isobath from GEBCO.
    2) A Gulf Stream contour (T=15C).
    3) A Shelf-Slope Front line (from SSFANNUAL data).

    The Gulf Stream temperature dataset is typically 277x457, 
    but here we interpolate it onto the potentially smaller GLORYS grid 
    so that the shapes match when calling ax.contour().
    """

    # ---------------------------
    # 1) Read GLORYS and Prevalence Data
    # ---------------------------
    glorys_dataset = nc.Dataset(glorys_fp, 'r')
    latitudes = glorys_dataset['latitude'][:]     # shape: (N_lat,)
    longitudes = glorys_dataset['longitude'][:]   # shape: (N_lon,)

    prevalence_data = nc.Dataset(prevalence_fp, 'r')
    data = np.array(prevalence_data['prevalence'][:])  # shape: (N_lat, N_lon)
    data = np.ma.filled(data, fill_value=0.0)          # Replace any masked with 0

    # ---------------------------
    # 2) Read SSF .mat Data (for shelf slope front)
    # ---------------------------
    mat_file_path = pkg_resources.resource_filename('glorys_plot_tools', 'contours/SSFANNUAL_1973_2017.mat')
    mat_file = scipy.io.loadmat(mat_file_path)
    SSFANNUAL = mat_file['SSFANNUAL']
    time = SSFANNUAL['time'][0, 0].flatten()
    lat = SSFANNUAL['lat'][0, 0]  # shape: (some_n, m)
    lon = SSFANNUAL['lon'][0, 0].flatten()

    # Remove rows that have NaN in lat
    valid_idx = ~np.isnan(lat).any(axis=1)
    lat = lat[valid_idx, :]
    time = time[valid_idx]
    start_idx = 249  # e.g. picking a start time index
    mean_lat = np.nanmean(lat[start_idx::12], axis=0)

    # ---------------------------
    # 3) Read and Subset GEBCO Data (bathymetry)
    # ---------------------------
    gebco_file_path = pkg_resources.resource_filename('glorys_plot_tools', 'contours/contour_200m_depth.nc')
    gebco_dataset = xr.open_dataset(gebco_file_path)
    gebco_latitudes = gebco_dataset['lat'].values        # shape: (G_lat,)
    gebco_longitudes = gebco_dataset['lon'].values       # shape: (G_lon,)
    gebco_altitudes = gebco_dataset['elevation'].values  # shape: (G_lat, G_lon)

    # Subset to GLORYS bounding box
    lat_mask = (gebco_latitudes >= latitudes.min()) & (gebco_latitudes <= latitudes.max())
    lon_mask = (gebco_longitudes >= longitudes.min()) & (gebco_longitudes <= longitudes.max())

    sub_gebco_lat = gebco_latitudes[lat_mask]
    sub_gebco_lon = gebco_longitudes[lon_mask]
    sub_gebco_alt = gebco_altitudes[lat_mask, :][:, lon_mask]

    # Build a meshgrid for the subset
    sub_gebco_lat_grid, sub_gebco_lon_grid = np.meshgrid(sub_gebco_lat, sub_gebco_lon, indexing='ij')

    # ---------------------------
    # 4) Read and Interpolate Temperature Data (Gulf Stream)
    #    from 277x457 -> match GLORYS shape
    # ---------------------------
    temp_file_path = pkg_resources.resource_filename('glorys_plot_tools', 'contours/temperature_data.nc')
    temp_dataset = xr.open_dataset(temp_file_path)
    # Suppose 'temperature' dims = ('latitude', 'longitude') = (277, 457)
    temps = temp_dataset['temperature'].isel(time=9)  # shape: (277, 457)

    # Interpolate onto the GLORYS lat/lon (N_lat, N_lon).
    # Make sure the dimension names in 'temps' match ('latitude','longitude') or similar.
    # If your temp data uses different dimension names (e.g., 'lat','lon'), adjust accordingly:
    lat_da = xr.DataArray(latitudes, dims=['latitude'])
    lon_da = xr.DataArray(longitudes, dims=['longitude'])

    # This interpolation step returns shape (N_lat, N_lon)
    temps_on_glorys = temps.interp(latitude=lat_da, longitude=lon_da, method='nearest')

    # Build a meshgrid for the GLORYS lat/lon (N_lat, N_lon)
    lat_grid, lon_grid = np.meshgrid(latitudes, longitudes, indexing='ij')  # shape: (N_lat, N_lon)

    # ---------------------------
    # 5) Plot Setup
    # ---------------------------
    plt.figure(figsize=(15, 10))
    ax = plt.axes(projection=ccrs.Mercator(min_latitude=latitudes.min(), max_latitude=latitudes.max()))
    cmap = cmocean.cm.haline
    cmap.set_bad(color='white')

    # Plot the prevalence data (N_lat x N_lon)
    img = ax.imshow(
        data,
        cmap=cmap,
        transform=ccrs.PlateCarree(),
        extent=[longitudes.min(), longitudes.max(), latitudes.min(), latitudes.max()],
        origin='lower'
    )
    ax.set_extent([longitudes.min(), longitudes.max(), latitudes.min(), latitudes.max()],
                  crs=ccrs.PlateCarree())

    # Add coastlines and land
    ax.coastlines()
    ax.add_feature(cfeature.LAND, zorder=1, edgecolor='black', facecolor='lightgrey')

    # ---------------------------
    # 6) Optional: 200m Depth Contour
    # ---------------------------
    if depth_contour == 1:
        matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
        # Subset grid is shape (sub_gebco_lat.size, sub_gebco_lon.size)
        gebco_depth_contour = ax.contour(
            sub_gebco_lon_grid,
            sub_gebco_lat_grid,
            sub_gebco_alt,
            levels=[-200],
            colors='cyan',
            linewidths=0.5,
            transform=ccrs.PlateCarree()
        )
        ax.clabel(gebco_depth_contour, inline=1, fontsize=0)

    # ---------------------------
    # 7) Optional: Gulf Stream Contour (T = 15°C)
    #     Use the regridded temps_on_glorys (shape: N_lat x N_lon)
    # ---------------------------
    if gulf_stream_contour == 1:
        glorys_temp_contour = ax.contour(
            lon_grid,             # shape (N_lat, N_lon)
            lat_grid,             # shape (N_lat, N_lon)
            temps_on_glorys,      # shape (N_lat, N_lon)
            levels=[15],
            colors='yellow',
            linewidths=0.5,
            transform=ccrs.PlateCarree()
        )
        ax.clabel(glorys_temp_contour, inline=1, fontsize=0)

    # ---------------------------
    # 8) Optional: Shelf-Slope Front (SSF)
    # ---------------------------
    if ssf_contour == 1:
        ax.plot(
            lon,
            mean_lat,
            color='red',
            linewidth=1,
            transform=ccrs.PlateCarree(),
            label='Shelf Slope Front (SSF)'
        )

    # Prepare legend handles based on requested overlays
    contour_handles = []
    if depth_contour == 1:
        contour_handles.append(mpatches.Patch(color='cyan', label='200m Depth'))
    if gulf_stream_contour == 1:
        contour_handles.append(mpatches.Patch(color='yellow', label='Gulf Stream'))
    if ssf_contour == 1:
        contour_handles.append(mpatches.Patch(color='red', label='Shelf-Slope Front (SSF)'))

    # ---------------------------
    # 9) Gridlines, Colorbar, Title, Legend
    # ---------------------------
    gl = ax.gridlines(
        draw_labels=True,
        crs=ccrs.PlateCarree(),
        linewidth=0.5,
        color='gray',
        alpha=0.5,
        linestyle='--'
    )
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 12, 'color': 'black'}
    gl.ylabel_style = {'size': 12, 'color': 'black'}

    cbar = plt.colorbar(img, ax=ax, label='Probability (%)', orientation='vertical')
    cbar.ax.tick_params(labelsize=10)

    if plot_title:
        plt.title(plot_title, fontsize=12)
    else:
        plt.title(f'Subsurface Duct Prevalence (%) {suffix}', fontsize=12)

    if contour_handles:
        plt.legend(handles=contour_handles, loc='upper right')

    # Save the figure
    save_path = f'duct_prevalence_{suffix}.png'
    plt.savefig(save_path, dpi=200)
    plt.show()

    # ---------------------------
    # 10) Close Datasets
    # ---------------------------
    glorys_dataset.close()
    prevalence_data.close()
    gebco_dataset.close()
    temp_dataset.close()
