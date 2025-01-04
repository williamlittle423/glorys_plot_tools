import netCDF4 as nc
import time

def create_dataset(suffix: str, glorys_fp: str) -> str:
    """ 
    Create a dataset with the following:
    Dimensions:
    - lat: len(glorys['latitude'])
    - lon: len(glorys['longitude'])

    Variables:
    - prevalence: 2D variable with dimensions (lat, lon)
    - longitude: 1D variable with dimensions (lon)
    - latitude: 1D variable with dimensions (lat)
    
    Returns:
    - The name of the created dataset
    """

    # Open the reference dataset
    GLORYS_reference = nc.Dataset(glorys_fp)

    # Create a new netCDF dataset for storing prevalence
    dataset = nc.Dataset(f"duct-prevalence_{suffix}.nc", "w")

    # Create dimensions
    dataset.createDimension("lat", len(GLORYS_reference['latitude']))
    dataset.createDimension("lon", len(GLORYS_reference['longitude']))

    # Create variables
    prevalence_var = dataset.createVariable("prevalence", "f4", ("lat", "lon"))
    longitude_var = dataset.createVariable("longitude", "f4", ("lon"))
    latitude_var = dataset.createVariable("latitude", "f4", ("lat"))

    # Copy over latitude and longitude data
    longitude_var[:] = GLORYS_reference['longitude'][:]
    latitude_var[:]   = GLORYS_reference['latitude'][:]

    # Close our newly created dataset (we'll open it again in append mode)
    dataset.close()
    GLORYS_reference.close()
    
    return f"duct-prevalence_{suffix}.nc"


def calculate_prevalence_at_coord(duct_dataset, lat_idx, lon_idx):
    """ 
    Calculate the prevalence of ducts at the given (lat_idx, lon_idx).
    Duct prevalence is defined as (number of times duct=1) / (total time steps) * 100.
    """
    duct_count = 0
    total_count = len(duct_dataset['time'])
    
    for t in range(total_count):
        if duct_dataset['duct'][t, lat_idx, lon_idx] == 1:
            duct_count += 1
    
    prevalence = duct_count / total_count * 100.0
    return prevalence


def calculate_total_prevalence(duct_filename, glorys_fp, suffix):
    """
    Create a prevalence dataset (via `create_dataset`) and fill its 'prevalence'
    variable by looping through each lat/lon in the duct data. 
    No multiprocessing is used here.
    """

    # 1) Create (or open) the prevalence dataset
    prevalence_filename = create_dataset(suffix, glorys_fp)

    # 2) Open duct data for reading and the new prevalence file for appending
    duct_data = nc.Dataset(duct_filename, 'r')
    prevalence_data = nc.Dataset(prevalence_filename, "a")
    
    lat_count = len(duct_data['latitude'])
    lon_count = len(duct_data['longitude'])
    
    start_time = time.time()
    print("Starting prevalence computation...")

    # 3) Nested loop over all lat/lon points
    for lat_idx in range(lat_count):
        for lon_idx in range(lon_count):
            # Calculate prevalence at this lat/lon
            prevalence = calculate_prevalence_at_coord(duct_data, lat_idx, lon_idx)
            # Write to the prevalence variable
            prevalence_data['prevalence'][lat_idx, lon_idx] = prevalence
    
    # 4) Close datasets
    prevalence_data.close()
    duct_data.close()

    end_time = time.time()
    print(f"Prevalence processing complete in {end_time - start_time:.2f} seconds.")
