import netCDF4 as nc
import numpy as np
from glorys_plot_tools.soundspeed import populate_c_vals
from scipy.signal import find_peaks
import multiprocessing
from tqdm import tqdm
import time
import netCDF4 as nc

def create_duct_file(filename, reference_file):
    """
    Create a netCDF4 file with the given filename by copying the 
    'time', 'latitude', and 'longitude' variables and their attributes 
    from the reference NetCDF file.

    Dimensions:
        - time
        - latitude
        - longitude

    Variables:
        - latitude(latitude)
        - longitude(longitude)
        - time(time)
        - duct(time, latitude, longitude)
        - duct_strength(time, latitude, longitude)
        - duct_width(time, latitude, longitude)
        - duct_depth(time, latitude, longitude)
        - duct_prominence(time, latitude, longitude)
        - cutoff_frequency(time, latitude, longitude)
    """
    try:
        # Open the reference file in read mode
        with nc.Dataset(reference_file, 'r') as ref:
            # Retrieve dimension sizes
            time_dim_size = len(ref.dimensions['time'])
            lat_dim_size = len(ref.dimensions['latitude'])
            lon_dim_size = len(ref.dimensions['longitude'])

            # Retrieve reference variables
            ref_time = ref.variables['time']
            ref_lat = ref.variables['latitude']
            ref_lon = ref.variables['longitude']

            # Open the new file in write mode
            with nc.Dataset(filename, 'w', format='NETCDF4') as new:
                # Create dimensions
                new.createDimension('time', time_dim_size)
                new.createDimension('latitude', lat_dim_size)
                new.createDimension('longitude', lon_dim_size)

                # Function to copy variables along with their attributes
                def copy_variable(var_name, new_file, ref_file, dimensions):
                    ref_var = ref_file.variables[var_name]
                    # Create the variable in the new file with the same data type and dimensions
                    new_var = new_file.createVariable(var_name, ref_var.datatype, dimensions)
                    # Copy variable attributes
                    for attr in ref_var.ncattrs():
                        new_var.setncattr(attr, ref_var.getncattr(attr))
                    # Copy variable data
                    new_var[:] = ref_var[:]

                # Copy 'time', 'latitude', and 'longitude' variables
                copy_variable('latitude', new, ref, ('latitude',))
                copy_variable('longitude', new, ref, ('longitude',))
                copy_variable('time', new, ref, ('time',))

                # Create additional variables
                # Define the dimensions for these variables
                var_dimensions = ('time', 'latitude', 'longitude')

                # Create 'duct' variable (assuming integer type)
                duct_var = new.createVariable('duct', 'i1', var_dimensions)
                duct_var.long_name = "Duct Indicator"
                duct_var.units = "unitless"  # Replace with actual units if available

                # Create 'duct_strength' variable
                duct_strength_var = new.createVariable('duct_strength', 'f4', var_dimensions)
                duct_strength_var.long_name = "Duct Strength"
                duct_strength_var.units = "units"  # Replace with actual units

                # Create 'duct_width' variable
                duct_width_var = new.createVariable('duct_width', 'f4', var_dimensions)
                duct_width_var.long_name = "Duct Width"
                duct_width_var.units = "meters"  # Replace with actual units

                # Create 'duct_depth' variable
                duct_depth_var = new.createVariable('duct_depth', 'f4', var_dimensions)
                duct_depth_var.long_name = "Duct Depth"
                duct_depth_var.units = "meters"  # Replace with actual units

                # Create 'duct_prominence' variable
                duct_prominence_var = new.createVariable('duct_prominence', 'f4', var_dimensions)
                duct_prominence_var.long_name = "Duct Prominence"
                duct_prominence_var.units = "units"  # Replace with actual units

                # Create 'cutoff_frequency' variable
                cutoff_frequency_var = new.createVariable('cutoff_frequency', 'f4', var_dimensions)
                cutoff_frequency_var.long_name = "Cutoff Frequency"
                cutoff_frequency_var.units = "Hz"  # Replace with actual units

                print(f"Successfully created NetCDF file '{filename}' with copied dimensions and variables.")

    except OSError as e:
        print(f"Error opening or modifying file: {e}")
    except KeyError as e:
        print(f"Missing variable or dimension in reference file: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


def duct_properties_calc(data, lat: int, lon: int, time_idx: int, DEPTH_LIMITER=30):
    """ Calculate the duct properties for the given time step"""
    
    depth_vals = np.array(data['depth'][:])
    c_vals = data['c'][time_idx, :, lat, lon]

    if DEPTH_LIMITER > len(c_vals):
        depth_vals = depth_vals[0:len(c_vals)]
    else:
        depth_vals = depth_vals[0:DEPTH_LIMITER]
        c_vals = c_vals[0:DEPTH_LIMITER]

    c_vals = np.array(c_vals)
    inverted_c = -c_vals
    min_peaks, _ = find_peaks(inverted_c)
    max_peaks, _ = find_peaks(c_vals)

    valid_min_peaks = [i for i in min_peaks if 10 < depth_vals[i] < 300]
    valid_max_peaks = [i for i in max_peaks if 10 < depth_vals[i] < 300]

    best_pair = None
    max_prominence = -np.inf

    for min_idx in valid_min_peaks:
        for max_idx in valid_max_peaks:
            if depth_vals[max_idx] > depth_vals[min_idx]:
                prominence = c_vals[max_idx] - c_vals[min_idx]
                if prominence > max_prominence:
                    max_prominence = prominence
                    best_pair = (min_idx, max_idx)

    if best_pair is None:
        return (False, 0, 0, 0, 0, 0)
    else:
        min_idx, max_idx = best_pair
        C_a = c_vals[max_idx]
        C_d = c_vals[min_idx]
        d_a = depth_vals[max_idx]
        d_d = depth_vals[min_idx]
        delta_C = np.abs(C_d - C_a)
        if delta_C < 0.1:
            return (False, 0, 0, 0, 0, 0)

        found_value = False
        loop_idx = min_idx - 1
        delta_z = None

        while loop_idx > 0:
            if c_vals[loop_idx] >= C_a:
                found_value = True
                d_above = depth_vals[loop_idx]
                delta_z = np.abs(d_above - d_a)
                break
            else:
                loop_idx -= 1

        if not found_value:
            for peak_idx in max_peaks:
                if depth_vals[peak_idx] < d_d:
                    found_value = True
                    d_above = depth_vals[peak_idx]
                    delta_z = np.abs(d_above - d_a)
                    break

        if not found_value:
            d_above = 0
            delta_z = np.abs(d_above - d_a)

        duct_strength = delta_C * delta_z
        cutoff_frequency = 3 / (8 * delta_z) * np.sqrt(np.power(C_d, 3) / (2 * delta_C))
        return (True, d_d, duct_strength, cutoff_frequency, delta_z, delta_C)  

def process_chunk(time_idx, lat_idx, lon_idx, glorys_file, duct_file):
    GLORYS_data = nc.Dataset(glorys_file, 'r')
    duct_values = duct_properties_calc(GLORYS_data, lat_idx, lon_idx, time_idx)
    GLORYS_data.close()
    return (time_idx, lat_idx, lon_idx, duct_values)

def add_data_to_duct_file(duct_file, results):
    duct_data = nc.Dataset(duct_file, 'a')
    for result in results:
        time_idx, lat_idx, lon_idx, duct_values = result
        duct, d_d, duct_strength, cutoff_frequency, delta_z, delta_C = duct_values
        duct_data['duct'][time_idx, lat_idx, lon_idx] = duct
        duct_data['duct_depth'][time_idx, lat_idx, lon_idx] = d_d
        duct_data['duct_strength'][time_idx, lat_idx, lon_idx] = duct_strength
        duct_data['cutoff_frequency'][time_idx, lat_idx, lon_idx] = cutoff_frequency
        duct_data['duct_width'][time_idx, lat_idx, lon_idx] = delta_z
        duct_data['duct_prominence'][time_idx, lat_idx, lon_idx] = delta_C

    duct_data.close()

def process_and_save(glorys_file, duct_file, time_idx):
    GLORYS_data = nc.Dataset(glorys_file, 'r')
    latitudes = GLORYS_data['latitude'][:]
    longitudes = GLORYS_data['longitude'][:]
    num_lat = len(latitudes)
    num_lon = len(longitudes)
    GLORYS_data.close()

    args = [(time_idx, lat_idx, lon_idx, glorys_file, duct_file) for lat_idx in range(num_lat) for lon_idx in range(num_lon)]

    with multiprocessing.Pool() as pool:
        results = pool.starmap(process_chunk, args)

    add_data_to_duct_file(duct_file, results)

def calculate_duct_properties(glorys_fp, duct_save_fp):

    glorys_data = nc.Dataset(glorys_fp, 'r')
    len_t = len(glorys_data['time'])
    glorys_data.close()

    # Create the duct file to hold data
    create_duct_file(duct_save_fp, glorys_fp)

    # Fill the existing glorys data with sound speed values
    print('Filling GLORYS dataset with sound speed values...')
    populate_c_vals(glorys_fp)

    # Calculate the duct properties for each time step and store them in the new duct dataset
    print('Beginning duct calculations...')
    for t in tqdm(range(len_t)):
        process_and_save(glorys_fp, duct_save_fp, t)