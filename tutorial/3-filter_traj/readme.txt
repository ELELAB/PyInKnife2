# Requirements:
# - Run the script in '../1-make_index'

# In this folder, the script '3-filter_traj.sh' is used
# to create two new trajectory files from the full-length
# trajectory. The new trajectories will contain only protein
# atoms but only subsets of the original trajectory's frames.
# To do so, the script reads in both the original trajectory
# file and the index file generated at the '1-make_index'
# step, which contains only the indexes of the protein atoms.

# The Bash line to run the script is:

./3-filter_traj.sh

# Be aware that the full-length trajectory is not
# available in the GitHub repository because it is
# too big to be hosted on GitHub. However, we can
# provide the trajectory and the associated topology
# (TPR) file upon request.
