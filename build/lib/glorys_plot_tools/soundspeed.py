from tqdm import tqdm
import os
import time
import netCDF4 as nc
import numpy as np
from multiprocessing import Pool, cpu_count

def sound_speed(D, S, T):
    """
    Mackenzie (1981) equation for sound speed in seawater
    
    Return:
        c = sound speed in metres

    Parameters:
        T = temperature in degrees Celsius
        S = salinity in parts per thousand
        D = depth in metres

    Range of validity: 2 < T < 30, 25 < S < 40, 0 < D < 8000
    """
    """if (T > 30.0 or T < 0.0):
        raise Exception("Temperature out of range for Mackenzie (1981) sound speed equation: {}".format(T))
        return None
    if (S > 40.0 or S < 25.0):
        raise Exception("Salinity out of range for Mackenzie (1981) sound speed equation: {}".format(S))
        return None
    if (D > 8000.0 or D < 0.0):
        raise Exception("Depth out of range for Mackenzie (1981) sound speed equation: {}".format(D))
        return None"""
    c = ( 
          1448.96 + 4.591*T - 5.304e-2 * T**2 + 2.374e-4 * T**3 
        + 1.340 * (S - 35.0) + 1.630E-2 * D + 1.675e-7 * D**2 
        - 1.025e-2 * T * (S-35.0) - 7.139e-13 * T * D**3
        )
    return c

def compute_sound_speed(thetao, so, depths):
    """
    Vectorized computation of sound speed. 
    thetao, so should be 3D arrays of shape (depth, lat, lon).
    depths is a 1D array corresponding to the depth dimension.
    """
    # Broadcast depths to match (depth, lat, lon)
    D = depths[:, np.newaxis, np.newaxis]

    # Extract variables
    T = thetao
    S = so

    # Compute c using the UNESCO equation (vectorized)
    c = (1448.96 
         + 4.591 * T
         - 5.304e-2 * T**2
         + 2.374e-4 * T**3
         + 1.340 * (S - 35.0)
         + 1.630e-2 * D
         + 1.675e-7 * D**2
         - 1.025e-2 * T * (S - 35.0)
         - 7.139e-13 * T * D**3)

    # Assign NaN where conditions are not met
    return c

def process_time(args):
    """
    Worker function for multiprocessing. 
    args: A tuple of (thetao, so, depths).
    Returns the computed sound speed array for that time.
    """
    thetao, so, depths = args
    return compute_sound_speed(thetao, so, depths)

def populate_c_vals(file_path):
    try:
        # Open the dataset in append mode
        with nc.Dataset(file_path, 'a', format='NETCDF4') as dataset:
            depths = dataset['depth'][:]
            num_t = len(dataset['time'])
            lat_dim = dataset.dimensions['latitude']
            lon_dim = dataset.dimensions['longitude']
            depth_dim = dataset.dimensions['depth']
            time_dim = dataset.dimensions['time']

            # Create the variable 'c' if it doesn't exist
            if 'c' not in dataset.variables:
                dataset.createVariable('c', 'f4', 
                                       (time_dim.name, depth_dim.name, lat_dim.name, lon_dim.name))

            # We will read data for each time step , send them to a pool process, and then write the result.
            tasks = []
            for t in range(num_t):
                thetao = dataset['thetao'][t]  # shape: (depth, lat, lon)
                so = dataset['so'][t]          # shape: (depth, lat, lon)
                tasks.append((thetao, so, depths))

            # Use a multiprocessing pool to process each time in parallel
            with Pool(processes=cpu_count()) as pool:
                for time_idx, result in tqdm(enumerate(pool.imap(process_time, tasks), start=0)):
                    # Write results to the NetCDF file
                    dataset.variables['c'][time_idx, :, :, :] = result

            print("Soundspeed successfully added to the NetCDF file.")

    except OSError as e:
        print(f"Error opening or modifying file: {e}")