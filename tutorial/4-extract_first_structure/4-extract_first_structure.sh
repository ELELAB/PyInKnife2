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

# Set the path to the input trajectory file
input_traj=../traj.xtc

# Set the path to the input TPR file
input_tpr=../top.tpr

# Set the path to the input index file
input_ndx=../1-make_index/protein.ndx

# Set the name of the output PDB file
output_pdb=first_structure.pdb

# Run the 'trjconv' command to extract the first structure from
# the full-length trajectory
gmx_mpi trjconv -f $input_traj -s $input_tpr -n $input_ndx -b 0 -e 1 -o $output_pdb
