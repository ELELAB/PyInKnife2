# This step cannot be directly reproduced inside this
# tutorial because the trajectory is too big to be uploaded to
# the GitHub repository. However, we can provide the full-length
# trajectory and the corresponding topology (TPR) file upon
# request to the users who want to reproduce the tutorial
# in full.

# You need to load GROMACS to manipulate the trajectory and the
# associated TPR file.
# In this case, we used GROMACS 5.1.2 with MPI support patched
# with PLUMED 2.3b (the PLUMED patch is not needed to run
# the commands we used here). We recommend to use either the
# same version of GROMACS or a higher version to run this
# tutorial.
source /usr/local/gromacs-5.1.2_plumed-2.3b/bin/GMXRC.bash

# Set the path to the input TPR file
input_tpr=../top.tpr

# Run the 'make_ndx' command to extract only the protein atoms
# (group 1, 'Protein') from the TPR file and write an index
# file containing only the indexes of those atoms
gmx_mpi make_ndx -f $input_tpr -o protein.ndx << eof
keep 1 
q
eof
