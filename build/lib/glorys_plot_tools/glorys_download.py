import copernicusmarine

def download_data(
    dataset_id="cmems_mod_glo_phy_my_0.083deg_P1D-m",
    dataset_version="202311",
    variables=["so", "thetao", "uo", "vo", "zos"],
    minimum_longitude=-65,
    maximum_longitude=-59.5-0.15*1.5,
    minimum_latitude=37.53+0.165,
    maximum_latitude=43,
    start_year=1993,
    start_month=1,
    start_day=1,
    end_year=1994,
    end_month=1,
    end_day=1,
    minimum_depth=0.49402499198913574,
    maximum_depth=0.49402499198913574,
    disable_progress_bar=False
):

    # GLORYS date format: 1993-01-01T00:00:00",
    start_month = f'0{start_month}' if start_month < 10 else str(start_month)
    start_day = f'0{start_day}' if start_day < 10 else str(start_day)

    end_month = f'0{end_month}' if end_month < 10 else str(end_month)
    end_day = f'0{end_day}' if end_day < 10 else str(end_day)

    start_datetime = f'{str(start_year)}-{start_month}-{start_day}T00:00:00'
    end_datetime = f'{str(end_year)}-{end_month}-{end_day}T00:00:00'
        
    if minimum_depth < 0.49402499198913574:
        minimum_depth = 0.49402499198913574
        if maximum_depth < minimum_depth:
            maximum_depth = minimum_depth
        print("Minimum depth must be greater than 0.494m. Setting minimum depth to 0.494m.")
    
    copernicusmarine.subset(
        dataset_id=dataset_id,
        dataset_version=dataset_version,
        variables=variables,
        minimum_longitude=minimum_longitude,
        maximum_longitude=maximum_longitude,
        minimum_latitude=minimum_latitude,
        maximum_latitude=maximum_latitude,
        start_datetime=start_datetime,
        end_datetime=end_datetime,
        minimum_depth=minimum_depth,
        maximum_depth=maximum_depth,
        disable_progress_bar=disable_progress_bar,
        username='wlittle',
        password='Finley2019!',
        force_download=True
    )
