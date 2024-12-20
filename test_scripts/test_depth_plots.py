from glorys_plot_tools.depth_profiles import plot_profiles

plot_profiles.plot_profiles(
    42.0,
    -62.0,
    1996,
    6,
    6,
    1996,
    8,
    7,
    300,
    animation=0,
    include_duct_depth=1,
    static_plot_title=None,
    static_save_path=None,
    static_include_mean=1,
    static_include_daily=1,
    static=1,
)