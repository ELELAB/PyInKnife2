This directory runs the PyinKnife pipeline in three scripts:

NB We removed the xtc files from the resampling to avoid issues with space in the repository

#prepare input files
ln -s ../../../4.frames/pdbmovie_1.pdb 


#step 1 use pyinknife_run
#you need to customize  the run.yalm file 
#have only salt bridge mode set to true
#define the cutoffs to scrutinize
#check the compatibility with PyInteraph2 flags and file names

pyinknife_run  -f ../../../3.filt_trjs/traj_prot_dt100.xtc -s pdbmovie_1.pdb -r pdbmovie_1.pdb -c run.yaml -d . -n 4

#at the end of the analyses the outputs are in the different resampling and fulltrj folders with the filtered graph, and two folders for
#connected components and hubs, respectively. Before starting the run verify that the .csv files from the previous step have been created 
ls */hc/*/*csv
#if they are not in the folders it means that there has been a problem at the previous step to solve first
pyinknife_aggregate -c run.yaml -ca aggregate.yaml -d . -od aggregate --firstccs 5 

#plot the results 
#connected components
pyinknife_plot -c run.yaml -ca aggregate.yaml -cp plot_ccs_barplot.yaml -p ccs -d aggregate -od plots
#hubs
pyinknife_plot -c run.yaml -ca aggregate.yaml -cp plot_hubs_barplot.yaml -p hubs -d aggregate -od plots
