import gsw
import netCDF4 as nc
import numpy as np

# Function to calculate density
def calculate_density(SP, CT, depth, lat, lon):
    p = gsw.p_from_z(-depth, lat)  # Convert depth to pressure
    SA = gsw.SA_from_SP(SP, p, lon, lat)  # Absolute Salinity
    return gsw.rho(SA, CT, p)  # Density

# Function to compute density for given coordinates and month
def compute_density_for_profile(SP, thetao, depths, lat, lon):
    density = np.zeros_like(depths)
    for k, depth in enumerate(depths):
        SP_value = SP[k]
        thetao_value = thetao[k]
        if not np.isnan(SP_value) and not np.isnan(thetao_value):
            CT_value = gsw.CT_from_pt(gsw.SA_from_SP(SP_value, gsw.p_from_z(-depth, lat), lon, lat), thetao_value)
            density[k] = calculate_density(SP_value, CT_value, depth, lat, lon)
        else:
            density[k] = np.nan
    return density


def fill_density(filename):

    # Fill the density variable in the climatology file
    with nc.Dataset(filename, 'a') as ds:
        # Check if the density variable exists, and create it if it does not
        if 'density' not in ds.variables:
            # Create the density variable
            density_var = ds.createVariable('density', 'f4', ('time', 'depth', 'latitude', 'longitude'))
            density_var.units = 'kg/m^3'  # Set appropriate units
            density_var.long_name = 'Sea Water Density'

        # Loop through time, latitudes, and longitudes - it is only one coordinate right now
        for t_idx in range(ds.variables['time'].shape[0]):
            lat = ds.variables['latitude'][0]
            lon = ds.variables['longitude'][0]
            SP = ds.variables['so'][t_idx, :, 0, 0]
            thetao = ds.variables['thetao'][t_idx, :, 0, 0]
            depths = ds.variables['depth'][:]                                                                                                                         
            density = compute_density_for_profile(SP, thetao, depths, lat, lon)
            ds.variables['density'][t_idx, :, 0, 0] = density
        
        print(f"Density filled for lat: {ds.variables['latitude'][0]:.1f} and lon{ds.variables['longitude'][0]:.1f}")