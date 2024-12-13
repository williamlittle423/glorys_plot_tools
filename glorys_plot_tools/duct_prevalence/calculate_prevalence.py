import netCDF4 as nc
import time
import multiprocessing

def create_dataset(suffix: str, glorys_fp: str) -> str:
    """ Create a dataset with the following:
    Dimensions:
    - lat: len(glorys['latitude'])
    - lon: len(glorys['longitude'])

    Variables:
    - prevalence: 3D variable with dimensions (time, lat, lon)
    - longitude: 1D variable with dimensions (lon)
    - latitude: 1D variable with dimensions (lat)
    
    Returns:
    - The name of the created dataset
    """

    # Create the dataaset
    # TODO: Update this to use the downloaded dataset
    GLORYS_reference = nc.Dataset(glorys_fp)

    dataset = nc.Dataset(f"duct-prevalence_{suffix}.nc", "w")

    dataset.createDimension("lat", len(GLORYS_reference['latitude']))
    dataset.createDimension("lon", len(GLORYS_reference['longitude']))

    dataset.createVariable("prevalence", "f4", ("lat", "lon"))
    longitude = dataset.createVariable("longitude", "f4", ("lon"))
    latitude = dataset.createVariable("latitude", "f4", ("lat"))

    longitude = GLORYS_reference['longitude']
    latitude = GLORYS_reference['latitude']

    dataset.close()
    
    return f"duct-prevalence_{suffix}.nc"

def calculate_prevalence_at_coord(duct_dataset, lat, lon):
    """ Calculate the prevalence of ducts given at latitude and longitude
    """
    duct_count = 0
    total_count = len(duct_dataset['time'])
    
    for t in range(total_count):
        if duct_dataset['duct'][t, lat, lon] == 1:
            duct_count += 1
    
    prevalence = duct_count / total_count * 100
    
    return prevalence

def process_chunk(duct_filename, lat_idx, lon_idx):
    duct_data = nc.Dataset(duct_filename, 'r')
    prevalence = calculate_prevalence_at_coord(duct_data, lat_idx, lon_idx)
    return (lat_idx, lon_idx, prevalence)

def calculate_total_prevalence(duct_filename, glorys_fp, suffix):

    prevalence_filename = create_dataset(suffix, glorys_fp)
    
    duct_data = nc.Dataset(duct_filename, 'r')
    prevalence_data = nc.Dataset(prevalence_filename, "a")
    
    args = [(duct_filename, lat_idx, lon_idx) for lat_idx in range(len(duct_data['latitude'])) for lon_idx in range(len(duct_data['longitude']))]
    
    start_time = time.time()
    print('Initializing multiprocessing pool...')
    with multiprocessing.Pool() as pool:
        results = pool.starmap(process_chunk, args)
        
    for result in results:
        lat_idx, lon_idx, prevalence = result
        prevalence_data['prevalence'][lat_idx, lon_idx] = prevalence
        
    print('Prevalence processing complete in: ', time.time() - start_time, 's')