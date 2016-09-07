#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# ----------------------------------------
# USAGE:
# ./calc_avg_structure.py pdb_file trajectory_location 1st_trajectory_# last_trajectory_# system_descripter avg_struct_filename
# ./calc_avg_structure.py config_file

# ----------------------------------------
# PREAMBLE:

import sys
import numpy as np
from numpy.linalg import *
import MDAnalysis
from MDAnalysis.analysis.align import *
from distance_functions import *

zeros = np.zeros
dot_prod = np.dot
sqrt = np.sqrt
flush = sys.stdout.flush

# ----------------------------------------
# VARIABLE DECLARATION

config_file = sys.argv[1]	# Local or Global positon of the config file that holds all the values for the parameters

necessary_parameters = [pdb_file,traj_loc,start,end,system,Wrapped]  # LEAVING OUT alignment, thresh, maxIter variables...

# ----------------------------------------
# SUBROUTINES:

def ffprint(string):		# Useful function to use when on a system that has a large buffer memory
	print '%s' %(string)
	flush()

def config_parser(config_file):	# Function to take config file and create/fill the parameter dictionary 
	parameters = {}
	for i in range(len(necessary_parameters)):
		parameters[necessary_parameters[i]] = ''

	# SETTING DEFAULT PARAMETERS FOR OPTIONAL PARAMETERS:
	parameters['alignment'] = 'protein'
	parameters['important'] = 'protein'
	parameters['substrate'] = 'protein'	# GENERALLY A VERY BAD SELECTION FOR THE SUBSTRATE... 
	parameters['thresh'] = 1E-5
	parameters['maxIter'] = 100
	parameters['write_dcd'] = False
	parameters['write_summary'] = False
	parameters['write_overview'] = False
	
	# GRABBING PARAMETER VALUES FROM THE CONFIG FILE:
	execfile(config_file,parameters)
	for key, value in parameters.iteritems:
		if value == '':
			print '%s has not been assigned a value. This variable is necessary for the script to run. Please declare this variable within the config file.'
			sys.exit()

def summary():
	with open('%03d.%03d.%s.avg_structure.summary' %(parameters['start'],parameters['end'],parameters['system']),'w') as f:
		f.write('Using MDAnalysis version: %s\n' %(MDAnalysis.version.__version__))
		f.write('To recreate this analysis, run this line:\n')
		for i in range(len(sys.argv)):
			f.write('%s ' %(sys.argv[i]))
		f.write('Used the parameters:\n')
		for key,value in parameters.iteritems:
			f.write('%s = %s\n' %(key,value))
		f.write('\n\n')

# ----------------------------------------
# MAIN PROGRAM:
# CREATING PARAMETER DICTIONARY
config_parser(config_file)

# INITIATE THE PDB FILE TO BE USED TO SAVE THE AVERAGE STRUCTURE
avg_pdb = MDAnalysis.Universe(parameters['pdb_file'])
avg_important = avg_pdb.select_atoms(parameters['important'])

# LOAD IN THE UNIVERSE TO BE ANALYZED AND INITIATE THE NECESSARY ATOM SELECTIONS
u = MDAnalysis.Universe(parameters['pdb_file'])
u_align = u.select_atoms(parameters['alignment'])
u_important = u.select_atoms(parameters['important'])
if parameters['Wrapped'] != True:
	u_substrate = u.select_atoms(parameters['substrate'])

# GRABBING IMPORTANT NUMBERS FROM THE UNIVERSE
u_important_atoms = u_important.n_atoms  	#len(u_important.atoms)
u_align_atoms = u_align.n_atoms			#len(u_align.atoms)
if parameters['Wrapped'] != True:
	u_substrate_res = u_substrate.n_residues	#len(u_substrate.residues)

# DETERMINING THE NUMBER OF STEPS TO BE AVERAGED OVER
nSteps = 0
start = parameters['start']
while start <= parameters['end']:
	u.load_new('%sproduction.%s/production.%s.dcd' %(parameters['traj_loc'],start,start))
	nSteps += len(u.trajectory)
	start += 1

ffprint('Number of steps to be averaged over : %d' %(nSteps))

# ARRAY DECLARATION
all_coord = zeros((nSteps,u_important_atoms,3),dtype=np.float64)
avgCoord = zeros((u_important_atoms,3),dtype=np.float64)
all_align = zeros((nSteps,u_align_atoms,3),dtype=np.float64)
avgAlign = zeros((u_align_atoms,3),dtype=np.float64)

# Trajectory Analysis: 
ffprint('Beginning trajectory analysis')
start = parameters['start']
temp = 0 
while start <= parameters['end']:
	ffprint('Loading trajectory %s' %(start))
	u.load_new('%sproduction.%s/production.%s.dcd' %(parameters['traj_loc'],start,start))

	for ts in u.trajectory:
		u_important.translate(-u_align.center_of_mass())
		
		if parameters['Wrapped'] != True:
			# CALCULATIONS that are unnecessary if the trajectory is wrapped.
			dims = u.dimensions[:3]		
			dims2 = dims/2.0

			for i in range(u_substrate_res):
				COM = u_substrate.residues[i].center_of_mass()
				t = wrapping(COM,dims,dims2)
				u_substrate.residues[i].atoms.translate(t)

		avgCoord += u_important.positions
		avgAlign += u_align.positions
		all_coord[temp] = u_important.positions
		all_align[temp] = u_align.positions
		temp += 1
	start += 1

ffprint(nSteps)
if temp != nSteps:
	ffprint('Failed to analyze all timesteps.')
	sys.exit()

avgCoord /= float(nSteps)
avgAlign /= float(nSteps)
ffprint('Finished with the trajectory analysis')

# Calculating and Aligning to the average positions
iteration = 0
residual = parameters['thresh'] + 10.0 					# arbitrary assignment greater than thresh
ffprint('Beginning iterative process of calculating average positions and aligning to the average')
while residual > parameters['thresh'] and iteration < parameters['maxIter']:
	tempAvgCoord = zeros((u_important_atoms,3),dtype=np.float64)		# zeroing out the tempAvgCoord array every iteration
	tempAvgAlign = zeros((u_align_atoms,3),dtype=np.float64)
	for i in range(nSteps):
		R, d = rotation_matrix(all_align[i,:,:],avgAlign)
		all_align[i,:,:] = dot_prod(all_align[i,:,:],R.T)
		all_coord[i,:,:] = dot_prod(all_coord[i,:,:],R.T)
		tempAvgAlign += all_align[i,:,:]
		tempAvgCoord += all_coord[i,:,:]			# recalculate the average coordinates to optimize the average position
	tempAvgCoord /= float(nSteps)				# finishing the average
	tempAvgAlign /= float(nSteps)
	residual = RMSD(avgAlign, tempAvgAlign, u_align_atoms)			# calculating the rmsd between avg and tempAvg to quantify our iterative optimization of the average positions	
	rmsd_all = RMSD(avgCoord, tempAvgCoord, u_important_atoms)
	iteration += 1
	avgCoord = tempAvgCoord
	avgAlign = tempAvgAlign
	ffprint('Steps: %d, RMSD btw alignment atoms: %e, RMSD btw all atoms: %e' %(iteration,residual,rmsd_all))
ffprint('Average structure has converged')				# Now have the iteratively aligned avgCoord array, as well as the iteratively aligned (COG-corrected and rotated) allCoord array

# Print out pdb of average structure
ffprint('Writing a pdb of the average structure.')
avg_important.positions = avgCoord
avg_important.write('%03d.%03d.%s.avg_structure.pdb' %(parameters['start'],parameters['end'],parameters['system']))
ffprint('Finished writing pdb of the average structure')

# PRINT out dcd frame of average structure; has more precision than the pdb format
if parameters['write_dcd'] == True:
	ffprint('Writing a dcd frame of the average structure.')
	avg_important.write('%03d.%03d.%s.avg_structure.dcd' %(parameters['start'],parameters['end']),parameters['system'])
	ffprint('Finished writing dcd of the average structure')

# APPENDING INFORMATION TO THE OVERVIEW FILE
if parameters['write_overview'] == True:
	with open('%s_averaging.output' %(parameters['system']),'a') as f:
		f.write('%d   %d   %d\n' %(parameters['start'],parameters['end'], nSteps))

if parameters['write_summary'] == True:
	summary()

