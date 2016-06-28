#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
##!/mnt/lustre_fs/users/mjmcc/apps/python2.7/bin/python
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

pdb_file = sys.argv[1]
traj_loc = sys.argv[2]
#start = int(sys.argv[3])
#end = int(sys.argv[4])
out_file = sys.argv[3]

#alignment = 'protein and name CA and (resid 20:25 or resid 50:55 or resid 73:75 or resid 90:94 or resid 112:116 or resid 142:147 or resid 165:169 or resid 190:194 or resid 214:218 or resid 236:240 or resid 253:258 or resid 303:307)'
alignment = 'name C*'

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

# LOAD IN THE UNIVERSE TO BE ANALYZED AND INITIATE THE NECESSARY ATOM SELECTIONS
u = MDAnalysis.Universe(pdb_file,traj_loc)
u_all = u.select_atoms('all')
#u_backbone = u.select_atoms('backbone')
u_align = u.select_atoms(alignment)
#u_important = u.select_atoms('protein or nucleic or resname A5 or resname A3 or resname U5 or resname atp or resname adp or resname PHX or resname MG')
#u_substrate = u.select_atoms('nucleic or resname A5 or resname A3 or resname U5 or resname atp or resname adp or resname PHX or resname MG')

# INITIATE THE PDB FILE TO BE USED TO SAVE THE AVERAGE STRUCTURE
avg_pdb = MDAnalysis.Universe(pdb_file)
avg_all = avg_pdb.select_atoms('all')
#avg_important = avg_pdb.select_atoms('protein or nucleic or resname A5 or resname A3 or resname U5 or resname atp or resname adp or resname PHX or resname MG')

# GRABBING IMPORTANT NUMBERS FROM THE UNIVERSE
u_important_atoms = len(u_all.atoms)
u_align_atoms = len(u_align.atoms)
#u_substrate_res = len(u_substrate.residues)

# DETERMINING THE NUMBER OF STEPS TO BE AVERAGED OVER
#temp = start
nSteps = 2
#while temp <= end:
#	u.load_new('%sproduction.%s/production.%s.dcd' %(traj_loc,start,start))
#	nSteps += len(u.trajectory)
#	temp += 1

ffprint('Number of steps to be averaged over : %d' %(nSteps))

# ARRAY DECLARATION
all_coord = zeros((nSteps,u_important_atoms,3),dtype=np.float32)
avgCoord = zeros((u_important_atoms,3),dtype=np.float32)
all_align = zeros((nSteps,u_align_atoms,3),dtype=np.float32)
avgAlign = zeros((u_align_atoms,3),dtype=np.float32)

# Trajectory Analysis: 
ffprint('Beginning trajectory analysis')
temp = 0 
#while start <= end:
#	ffprint('Loading trajectory %s' %(start))
#	u.load_new('%sproduction.%s/production.%s.dcd' %(traj_loc,start,start))

for ts in u.trajectory:
#		dimensions = u.dimensions[:3]

#	u_all.translate(-u_all.center_of_mass())

#		for i in range(u_substrate_res):
#			COM = np.zeros(3)
#			COM = u_substrate.residues[i].center_of_mass()
#			t = wrapping(COM,dimensions)
#			u_substrate.residues[i].atoms.translate(t)

#	u_all.write('%02d.structure.pdb' %(ts.frame))
	avgCoord += u_all.positions
	avgAlign += u_align.positions
	all_coord[temp] = u_all.positions
	all_align[temp] = u_align.positions
	temp += 1
#	start += 1

ffprint(nSteps)
if temp != nSteps:
	ffprint('Failed to analyze all timesteps; fucked shit up')

avgCoord /= float(nSteps)
avgAlign /= float(nSteps)
ffprint(len(all_coord))
ffprint(len(all_coord[0]))
ffprint('Finished with the trajectory analysis')

# Calculating and Aligning to the average positions
iteration = 0
residual = thresh + 10.0 					# arbitrary assignment greater than thresh
ffprint('Beginning iterative process of calculating average positions and aligning to the average')
while residual > thresh and iteration < maxIter:		
	tempAvgCoord = zeros((u_important_atoms,3),dtype=np.float32)		# zeroing out the tempAvgCoord array every iteration
	tempAvgAlign = zeros((u_align_atoms,3),dtype=np.float32)
	for i in range(nSteps):
		R, d = rotation_matrix(all_align[i,:,:],avgAlign)
		all_align[i,:,:] = dot_prod(all_align[i,:,:],R.T)
		all_coord[i,:,:] = dot_prod(all_coord[i,:,:],R.T)
		tempAvgAlign += all_align[i,:,:]
		tempAvgCoord += all_coord[i,:,:]			# recalculate the average coordinates to optimize the average position
	tempAvgCoord /= float(nSteps)				# finishing the average
	tempAvgAlign /= float(nSteps)
	residual = RMSD(avgAlign, tempAvgAlign)			# calculating the rmsd between avg and tempAvg to quantify our iterative optimization of the average positions	
	rmsd_all = RMSD(avgCoord, tempAvgCoord)
	iteration += 1
	avgCoord = tempAvgCoord
	avgAlign = tempAvgAlign
	ffprint('Steps: %d, RMSD btw alignment atoms: %e, RMSD btw all atoms: %e' %(iteration,residual,rmsd_all))
ffprint('Average structure has converged')				# Now have the iteratively aligned avgCoord array, as well as the iteratively aligned (COG-corrected and rotated) allCoord array

# Print out pdb of average structure
ffprint('Writing a pdb of the average structure.')
#avg_important.residues.set_positions(avgCoord)
avg_all.set_positions(avgCoord)
avg_all.write('test.avg_structure.pdb')
ffprint('Finished writing pdb of the average structure')

u_all.set_positions(all_coord[0,:,:])
u_all.translate(-u_all.center_of_mass())
u_all.write('00.frame.pdb')
u_all.set_positions(all_coord[1,:,:])
u_all.write('01.frame.pdb')


# APPENDING INFORMATION TO A SUMMARY FILE
out = open('%s' %(out_file),'a')
out.write('Test: %d\n' %(nSteps))
out.close()

