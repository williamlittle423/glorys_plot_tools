from glorys_plot_tools import glorys_download

glorys_download.download_data(
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
    end_year=1993,
    end_month=1,
    end_day=2,
    minimum_depth=0.49402499198913574,
    maximum_depth=0.49402499198913574,
    disable_progress_bar=False
)