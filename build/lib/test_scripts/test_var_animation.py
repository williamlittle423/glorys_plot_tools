from glorys_plot_tools.var_animations.plot_animation import plot_animation

plot_animation(
    var='thetao',
    lat_min=35,
    lat_max=45,
    lon_min=-62,
    lon_max=-58,
    start_year=2005,
    start_month=6,
    start_day=26,
    end_year=2005,
    end_month=7,
    end_day=26,
    fps=4,
    dpi=100,
    depth_m=0,
    cmocean_cmap='jet'
    #existing_fp='cmems_mod_glo_phy_my_0.083deg_P1D-m_multi-vars_60.00W-42.00W_35.00N-45.00N_0.49m_2002-06-26-2002-07-26.nc'
)