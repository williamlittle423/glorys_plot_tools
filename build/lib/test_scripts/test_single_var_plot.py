from glorys_plot_tools.plot_var.plot_glorys_variable import plot_glorys_variable

# Test single time plot
plot_glorys_variable(
    var='thetao',
    start_year=2019,
    start_month=6,
    start_day=26,
    lat_min=35,
    lat_max=45,
    lon_min=-68,
    lon_max=-57,
    end_year=2019,
    end_month=6,
    end_day=30,
    depth_m=0,
    title='Sea Surface Temperature - June 26-30, 2019',
    dataset_id="cmems_mod_glo_phy_my_0.083deg_P1D-m",
    dataset_version="202311",
    depth_contour=1,
    ssf_contour=1,
    gulf_stream_contour=1,
    save_path=None,
    cmap="dense",
    show_plot=True
)