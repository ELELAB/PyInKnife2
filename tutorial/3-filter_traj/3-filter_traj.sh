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

# Set the path to the input index file
input_ndx=../1-make_index/protein.ndx

# Set the name of the first output trajectory file, where a frame
# from the original trajectory will be written only when
# t MOD 100 = the first time point recorded
output_traj_dt100=traj_prot_dt100.xtc

# Set the name of the first output trajectory file, where a frame
# from the original trajectory will be written only when
# t MOD 1000 = the first time point recorded
output_traj_dt1000=traj_prot_dt1000.xtc

# Run the 'trjcat' command to create the first output trajectory
# file (-dt 100)
gmx_mpi trjcat -f $input_traj -n $input_ndx -o $output_traj_dt100 -dt 100

# Run the 'trjcat' command to create the second output trajectory
# file (-dt 1000)
gmx_mpi trjcat -f $input_traj -n $input_ndx -o $output_traj_dt1000 -dt 1000
