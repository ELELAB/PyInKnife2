# Requirements:
# - Run the script in '../1-make_index'

# In this folder, the script '2-filter_tpr.sh' is used to
# create a new TPR file from the TPR associated with the
# full-length trajectory. The new TPR file will contain
# only the protein atoms. To do so, the script reads in
# both the original TPR file and the index file generated
# at the '1-make_index' step, which contains only
# the indexes of protein atoms.

# The Bash line to run the script is:

./2-filter_tpr.sh

# Be aware that the full-length trajectory is not
# available in the GitHub repository because it is
# too big to be hosted on GitHub. However, we can
# provide the trajectory and the associated topology
# (TPR) file upon request.