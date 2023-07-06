# In this folder, we showcase different types of PSN analyses run
# with PyInKnife. Each analysis is run in a different folder:

# - acPSN -> PyInKnife2 is run to construct only acPSNs
#            on the '../../3-filter_traj/traj_prot_dt1000.xtc'
#            trajectory. The raw data are then aggregated, and
#            the aggregated data are plotted.

# - cmPSN -> PyInKnife' is run to construct only cmPSNs on
#            the '../../3-filter_traj/traj_prot_dt100.xtc'
#            trajectory. The raw data are then aggregated, and
# 	     the aggregated data are plotted.

# - hb -> PyInKnife' is run to construct only hbIINs on the
#         '../../3-filter_traj/traj_prot_dt1000.xtc' trajectory.
#         The raw data are then aggregated, and the aggregated
#         data are plotted.

# - hc -> PyInKnife2 is run to construct only hcIINs on the
#         '../../3-filter_traj/traj_prot_dt100.xtc' trajectory.
#         The raw data are then aggregated, and the aggregated
#         data are plotted.

# - sb -> PyInKnife2 is run to construct only sbIINs on the
#         '../../3-filter_traj/traj_prot_dt100.xtc' trajectory.
#         The raw data are then aggregated, and the aggregated
#         data are plotted.

# These configuration files are used for all analyses:

# - aggregate.yaml -> The configuration file used to aggregate the raw
#                     data.

# - plot_ccs_barplot.yaml -> The configuration file used to plot the
#                            number of nodes found in the most populated
#                            connected components.

# - plot_hubs_barplot.yaml -> The configuration file used to plot the
#                             degree distribution of the hub nodes.