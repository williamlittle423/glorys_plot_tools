from glorys_plot_tools.depth_profiles import plot_profiles

plot_profiles.plot_profiles(
    45.0,
    -52.0,
    1996,
    6,
    6,
    1996,
    6,
    10,
    300,
    animation=1,
    include_duct_depth=1,
    static_plot_title=None,
    static_save_path=None,
    static_include_mean=1,
    static_include_daily=1,
    static=1,
    existing_fp='cmems_mod_glo_phy_my_0.083deg_P1D-m_multi-vars_52.00W_45.00N_0.49-266.04m_1996-06-06-1996-06-10.nc'
)