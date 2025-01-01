from glorys_plot_tools.duct_calculations.calculate_ducts import calculate_ducts

calculate_ducts(
    lat_min=42,
    lat_max=45,
    lon_min=-62,
    lon_max=-59,
    start_year=2019,
    start_month=6,
    start_day=6,
    end_year=2019,
    end_month=6,
    end_day=7,
    dataset_id="cmems_mod_glo_phy_my_0.083deg_P1D-m",
    dataset_version="202311",
    # existing_fp='cmems_mod_glo_phy_my_0.083deg_P1D-m_multi-vars_62.00W-59.00W_42.00N-45.00N_0.49-318.13m_2019-06-06-2019-06-20.nc',
    ducts_save_fp=None
)