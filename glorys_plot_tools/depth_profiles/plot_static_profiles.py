import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import xarray as xr
from glorys_plot_tools.density import fill_density
from glorys_plot_tools.soundspeed import sound_speed

# Function to plot four depth profiles on a single figure with mean monthly line in color and daily lines in increasing opacity grey on each graph
def plot_mean_depth_profiles_w_daily(glorys_fp, duct_filename, daily_lists, save_path, suffix, plot_duct, static_include_mean, static_include_daily, title=None):
    
    with nc.Dataset(glorys_fp, 'r') as ds:
        depths = ds.variables['depth'][0:30]
        temperature = np.mean(ds.variables['thetao'][:, 0:30, 0, 0], axis=0)
        salinity = np.mean(ds.variables['so'][:, 0:30, 0, 0], axis=0)      
        density = np.mean(ds.variables['density'][:, 0:30, 0, 0], axis=0)
        sound_speed = np.mean(ds.variables['c'][:, 0:30, 0, 0], axis=0)
        
        T_length = len(ds['time'])

        # Duct data
        if plot_duct:
            with nc.Dataset(duct_filename, 'r') as duct_data:
                # Duct depth from mean of duct depths that are not zero
                duct_depth = np.mean(duct_data['duct_depth'][:, 0, 0][duct_data['duct_depth'][:, 0, 0] != 0])

    fig, axes = plt.subplots(1, 4, figsize=(18, 12), sharey=True)
    fig.subplots_adjust(hspace=0.07, wspace=0.05, top=1.2, right=0.85)
    fig.tight_layout(pad=3.0)

    # Used to make daily lines with increasing opacity
    alphas = np.linspace(0.1, 0.5, T_length)

    # Plot temperature profile
    if static_include_daily == 1:
        for i in range(T_length):
            axes[0].plot(daily_lists['thetao'][i, :], depths, color='grey', alpha=alphas[i])
    if static_include_mean == 1:
        axes[0].plot(temperature, depths, label='Temperature')
    if plot_duct == 1:
        axes[0].axhline(y=duct_depth, color='black', linestyle='--', label=f'Mean Duct Depth ({duct_depth:.1f} m)')
    axes[0].set_xlabel('Temperature (°C)', fontsize=14)
    axes[0].set_ylabel('Depth (m)', fontsize=14)
    axes[0].invert_yaxis()
    axes[0].tick_params(axis='x', labelsize=14)   
    axes[0].set_xlim([np.nanmin(daily_lists['thetao']), np.nanmax(daily_lists['thetao'])])
    axes[0].tick_params(axis='y', labelsize=14)   

    # Plot salinity profile
    if static_include_daily == 1:
        for i in range(T_length):
            axes[1].plot(daily_lists['so'][i, :], depths, color='grey', alpha=alphas[i])
    if static_include_mean == 1:
        axes[1].plot(salinity, depths, label='Salinity', color='orange')
    if plot_duct == 1:
        axes[1].axhline(y=duct_depth, color='black', linestyle='--', label=f'Mean Duct Depth ({duct_depth:.1f} m)')
    axes[1].set_xlabel('Salinity (PSU)', fontsize=14)
    axes[1].tick_params(axis='x', labelsize=14)   
    axes[1].set_xlim([np.nanmin(daily_lists['so']), np.nanmax(daily_lists['so'])])

    # Plot density profile
    if static_include_daily == 1:
        for i in range(T_length):
            axes[2].plot(daily_lists['density'][i, :], depths, color='grey', alpha=alphas[i])
    if static_include_mean == 1:
        axes[2].plot(density, depths, label='Density', color='green')
    if plot_duct == 1:
        axes[2].axhline(y=duct_depth, color='black', linestyle='--', label=f'Mean Duct Depth ({duct_depth:.1f} m)')
    axes[2].set_xlabel('Density (kg/m³)', fontsize=14)
    axes[2].tick_params(axis='x', labelsize=14)   
    axes[2].set_xlim([np.nanmin(daily_lists['density']), np.nanmax(daily_lists['density'])])

    # Plot sound speed profile
    if static_include_daily == 1:
        for i in range(T_length):
            axes[3].plot(daily_lists['c'][i], depths, color='grey', alpha=alphas[i])
    if static_include_mean == 1:
        axes[3].plot(sound_speed, depths, color='red')
    if plot_duct == 1:
        axes[3].axhline(y=duct_depth, color='black', linestyle='--', label=f'Mean Duct Depth ({duct_depth:.1f} m)')
    axes[3].set_xlabel('Sound Speed (m/s)', fontsize=14)
    axes[3].tick_params(axis='x', labelsize=14)   
    axes[3].set_xlim([np.nanmin(daily_lists['c']), np.nanmax(daily_lists['c'])])
    axes[3].legend(fontsize=14)
    
    if title:
        title = title + f' {suffix}'
    else:
        title = f'Glorys Depth Profiles {suffix}'
    
    fig.suptitle(title, fontsize=18, y=1.0)

    if save_path:
        plt.savefig(save_path, dpi=230)
    else:
        plt.savefig(f'static_depth_profiles_{suffix}.png', dpi=230)
    
    plt.show()

def plot_static_profile(glorys_fp, duct_fp, save_path, title, plot_duct, suffix, static_include_mean, static_include_daily):

    # Density is not a variable in the GLORYS data, so we need to calculate it
    fill_density(glorys_fp)

    glorys_data = nc.Dataset(glorys_fp, 'r')
    daily_temperature = glorys_data['thetao'][:, 0:30, 0, 0]
    daily_salinity = glorys_data['so'][:, 0:30, 0, 0]
    daily_density = glorys_data['density'][:, 0:30, 0, 0]

    daily_sound_speed = []
    for i in range(glorys_data['thetao'].shape[0]):
        daily_sound_speed.append([sound_speed(d, s, t) for d, s, t in zip(glorys_data['depth'][:], daily_salinity[i, :], daily_temperature[i, :])])
    daily_lists = {'thetao': daily_temperature, 'so': daily_salinity, 'density': daily_density, 'c': daily_sound_speed}
    so = np.mean(glorys_data['so'][:, :, 0, 0], axis=0)
    thetao = np.mean(glorys_data['thetao'][:, :, 0, 0], axis=0)
    c = [sound_speed(d, s, t) for d, s, t in zip(glorys_data['depth'][:], so, thetao)]
    plot_mean_depth_profiles_w_daily(glorys_fp, duct_fp, daily_lists, save_path, suffix, plot_duct, static_include_mean, static_include_daily, title=title)
