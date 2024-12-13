from glorys_plot_utils.duct_prevalence.plotting_util import plot_probability_with_contours

glorys_fp = 'cmems_mod_glo_phy_my_0.083deg_P1D-m_multi-vars_80.00W-42.00W_30.00N-53.00N_0.49-380.21m_1993-06-01-1993-06-05.nc'
prevalence_fp = 'duct-prevalence_80.00W-42.00W_30.00N-53.00N_0.49-380.21m_1993-06-01-1993-06-05.nc'
suffix = '1993-06-01-1993-06-05'
depth_contour = 1
ssf_contour = 1
gulf_stream_contour = 1

plot_probability_with_contours(glorys_fp, prevalence_fp, suffix, depth_contour, ssf_contour, gulf_stream_contour)