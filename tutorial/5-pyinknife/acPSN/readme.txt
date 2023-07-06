# Requirements:
# - Run the script in '../../4-extract_first_structure.sh'

# In this folder, we:
# - Run the PyInKnife pipeline.
# - Aggregate the raw data.
# - Plot the aggregated data.

# Step 1 - Generate the raw data

# We run using the 'run.yaml' configuration file (-c run.yaml) in the
# current working directory (-d .) with four processes (-n 4)

pyinknife_run -f ../../3-filter_traj/traj_prot_dt1000.xtc -s ../../4-extract_first_structure/first_structure.pdb -r ../../4-extract_first_structure/first_structure.pdb -c run.yaml -d . -n 4

# At the end of the run, we should have one folder containing the data for
# the full trajectory (whose name is defined in the 'run.yaml' file) and
# as many folders as the resamplings performed (whose names are also defined
# in the 'run.yaml' file).

# In this tutorial, we removed the trajectory (.xtc) files containing the
# resampled trajectories because of space limitations in the GitHub
# repository.

# Step 2 - Aggregate the data

# We run using the 'run.yaml' configuration file (-c run.yaml) to get the
# options used for generating the raw data and the '../aggregate.yaml' file
# to get the options that will be used for the aggregation (-ca ../aggregate.yaml).
# We assume that pyinknife_run was run in the current working directory (-d .),
# and we write the files containing the aggregate data in a directory named
# 'aggregate'.
# Only the results for the five most populated connected components will be
# aggregated (less populated connected components will be ignored)
# (--firstccs 5).

# First, we create the 'aggregate' directory. Then, we run.

mkdir aggregate
pyinknife_aggregate -c run.yaml -ca ../aggregate.yaml -d . -od aggregate --firstccs 5 

# Step 3 - plot the results

# First, we plot the degree distribution of the hubs found in the PSNs
# built from the trajectory resamplings (-p hubs).
# We run using the 'run.yaml' configuration file (-c run.yaml) to get the
# options used for generating the raw data, the '../aggregate.yaml' file to
# get the options used for aggregating the data (-ca ../aggregate.yaml), and
# the '../plot_hubs_barplot.yaml' to get the options that will be used to
# generate the plots.
# We generate the plots from the aggregated data contained in the files
# stored inside the 'aggregate' directory (-d aggregate), and we generate
# the plots in a directory named 'plots'.

# First, we create the 'plots' directory. Then, we run.

mkdir plots
pyinknife_plot -c run.yaml -ca ../aggregate.yaml -cp ../plot_hubs_barplot.yaml -p hubs -d aggregate -od plots

# Then, we plot the number of nodes found in each of the five most populated
# connected components in the PSNs built from the trajectory resamplings
# (-p ccs).
# We run using the 'run.yaml' configuration file (-c run.yaml) to get the
# options used for generating the raw data, the '../aggregate.yaml' file to
# get the options used for aggregating the data (-ca ../aggregate.yaml), and
# the '../plot_ccs_barplot.yaml' to get the options that will be used to
# generate the plots.
# We generate the plots from the aggregated data contained in the files
# stored inside the 'aggregate' directory (-d aggregate), and we generate
# the plots in a directory named 'plots'.

pyinknife_plot -c run.yaml -ca ../aggregate.yaml -cp ../plot_ccs_barplot.yaml -p ccs -d aggregate -od plots