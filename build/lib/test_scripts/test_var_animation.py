from glorys_plot_tools.var_animations.plot_animation import plot_animation

plot_animation(
    var='thetao',
    lat_min=38,
    lat_max=40,
    lon_min=-62,
    lon_max=-57,
    start_year=2021,
    start_month=6,
    start_day=16,
    end_year=2021,
    end_month=6,
    end_day=30,
    fps=4,
    dpi=150,
    depth_m=0,
    cmocean_cmap='jet',
)
