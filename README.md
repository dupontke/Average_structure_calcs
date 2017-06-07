# Average_structure_calcs
Scripts to average the structure of a substrate from a MD trajectory 

To run the average structure calculations, change the NPROD= __ to the number of production runs you want to calculate the average structure for in the calc_avg_structure.job file. Then use the following command:
./calc_avg_structure.job 

Once the .job command is complete you need to run the weighted average script in order to get the average structure of the entire equilibrated trajectory. 

Use the following command:
./weighted_average.py calc_avg_structure.output pdb_file system
