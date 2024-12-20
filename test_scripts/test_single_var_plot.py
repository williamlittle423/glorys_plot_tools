from glorys_plot_tools.plot_var.plot_glorys_variable import plot_glorys_variable

# Test single time plot
plot_glorys_variable(
    var='zos',
    start_year=2019,
    start_month=6,
    start_day=26,
    lat_min=30,
    lat_max=53,
    lon_min=-80,
    lon_max=-42,
    end_year=None,
    end_month=None,
    end_day=None,
    depth_m=0,
    title=None,
    dataset_id="cmems_mod_glo_phy_my_0.083deg_P1D-m",
    dataset_version="202311",
    depth_contour=1,
    ssf_contour=1,
    gulf_stream_contour=1,
    save_path=None,
    cmap="dense",
    show_plot=True
)

# Test time average plot
plot_glorys_variable(
    var='thetao',
    start_year=2019,
    start_month=6,
    start_day=26,
    end_year=2019,
    end_month=7,
    end_day=26,
    lat_min=30,
    lat_max=53,
    lon_min=-80,
    lon_max=-42,
    depth_m=0,
    title=None,
    dataset_id="cmems_mod_glo_phy_my_0.083deg_P1D-m",
    dataset_version="202311",
    depth_contour=1,
    ssf_contour=1,
    gulf_stream_contour=1,
    save_path=None,
    cmap="jet",
    show_plot=True
)