#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# ----------------------------------------
# USAGE:

# ----------------------------------------
# PREAMBLE:

import sys
import numpy as np
from numpy.linalg import *
import MDAnalysis
from MDAnalysis.analysis.align import *
from distance_functions import *

# ----------------------------------------
# VARIABLE DECLARATION

input_file = sys.argv[1]		# Format: 	start_traj	end_traj	number_of_frames
pdb_file = sys.argv[2]			# Read in a pdb file to be used for atom selections and assigning positions
system = sys.argv[3]                    # system of trajectory

alignment = 'protein and name CA and (resid 15:23 or resid 50:60 or resid 68:73 or resid 79:86 or resid 90:96 or resid 117:120 or resid 135:138 or resid 150:161 or resid 170:175 or resid 182:193 or resid 197:199 or resid 211:214 or resid 222:236)'

zeros = np.zeros
dot_prod = np.dot
sqrt = np.sqrt
flush = sys.stdout.flush

thresh = 1E-5
maxIter = 100

# ----------------------------------------
# SUBROUTINES:

def ffprint(string):
	print '%s' %(string)
	flush()

# ----------------------------------------
# MAIN PROGRAM:

# LOAD IN AVERAGE SUMMARY FILE AND CALCULATE THE IMPORTANT NUMBERS
average_list = np.loadtxt('%s' %(input_file))

nAverages = len(average_list)
nSteps = np.sum(average_list[:,2])

# LOAD IN THE PDB FILE TO BE USED AS THE AVERAGE STRUCTURE
u = MDAnalysis.Universe(pdb_file)
u_all = u.select_atoms('all')
u_align = u.select_atoms(alignment)
u_important = u.select_atoms('protein or resname GTP or resname SAH or resname MG')
u_substrate = u.select_atoms('resname GTP or resname SAH or resname MG')
pos0 = u_align.positions

# GRABBING IMPORTANT NUMBERS FROM THE UNIVERSE
u_important_atoms = u_important.n_atoms		#len(u_important.atoms)
u_align_atoms = u_align.n_atoms			#len(u_align.atoms)
u_substrate_res = u_substrate.n_residues	#len(u_substrate.residues)

# ARRAY DECLARATION
all_coord = zeros((nAverages,u_important_atoms,3),dtype=np.float64)
avgCoord = zeros((u_important_atoms,3),dtype=np.float64)
all_align = zeros((nAverages,u_align_atoms,3),dtype=np.float64)
avgAlign = zeros((u_align_atoms,3),dtype=np.float64)

# AVERAGE STRUCTURE ANALYSIS
ffprint('Beginning the averaging of averages process')
for i in range(nAverages):
	ffprint('Loading in the average structure of Trajectories %d to %d' %(average_list[i][0],average_list[i][1]))
	temp = MDAnalysis.Universe('%03d.%03d.avg_structure/%03d.%03d.%s.avg_structure.pdb' %(average_list[i][0],average_list[i][1],average_list[i][0],average_list[i][1],system))
	# INITIATING ATOM SELECTIONS
	temp_all = temp.select_atoms('all')
	temp_align = temp.select_atoms(alignment)
	temp_important = temp.select_atoms('protein or resname GTP or resname SAH or resname MG')
	temp_substrate = temp.select_atoms('resname GTP or resname SAH or resname MG')

	# TRANSLATING AVERAGE STRUCTURES TO THE PDB STRUCTURE
	temp_all.translate(-temp_align.center_of_mass())

	# UNWEIGHTING THE AVERAGE STRUCTURE
	all_coord[i] = temp_important.positions
	all_align[i] = temp_align.positions
	avgCoord += temp_important.positions*average_list[i][2]
	avgAlign += temp_align.positions*average_list[i][2]

avgCoord /= float(nSteps)
avgAlign /= float(nSteps)
ffprint('Finished collecting the average coordinates and averaging the averages...')

# Calculating and Aligning to the average positions
iteration = 0
residual = thresh + 10.0 					# arbitrary assignment greater than thresh
ffprint('Beginning iterative process of calculating average positions and aligning to the average')
while residual > thresh and iteration < maxIter:		
	tempAvgCoord = zeros((u_important_atoms,3),dtype=np.float64)		# zeroing out the tempAvgCoord array every iteration
	tempAvgAlign = zeros((u_align_atoms,3),dtype=np.float64)
	for i in range(nAverages):
		R, d = rotation_matrix(all_align[i,:,:],avgAlign)
		all_align[i,:,:] = dot_prod(all_align[i,:,:],R.T)
		all_coord[i,:,:] = dot_prod(all_coord[i,:,:],R.T)
		tempAvgAlign += all_align[i,:,:]
		tempAvgCoord += all_coord[i,:,:]
	tempAvgCoord /= float(nAverages)
	tempAvgAlign /= float(nAverages)
	residual = RMSD(avgAlign,tempAvgAlign,u_align_atoms)
	rmsd_all = RMSD(avgCoord,tempAvgCoord,u_important_atoms)
	iteration += 1
	avgCoord = tempAvgCoord
	avgAlign = tempAvgAlign
	ffprint('Steps: %d, RMSD btw alignment atoms: %e, RMSD btw all atoms: %e' %(iteration,residual,rmsd_all))

ffprint('Average structure has converged')

u_important.positions = avgCoord

with MDAnalysis.Writer('%03d.%03d.average_structure.pdb' %(average_list[0][0],average_list[-1][1]),u_important_atoms) as W:
	W.write(u_important)
ffprint('Finished writing pdb of the average structure')

with MDAnalysis.Writer('%03d.%03d.average_structure.dcd' %(average_list[0][0],average_list[-1][1]),u_important_atoms) as W:
	W.write(u_important)
ffprint('Finished writing dcd of the average structure')

