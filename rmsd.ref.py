#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# ----------------------------------------
# USAGE:


# ----------------------------------------
# PREAMBLE:

import numpy as np
import sys
import os
import MDAnalysis
from MDAnalysis.analysis.align import *
import sel_list
from distance_functions import *

pdb = sys.argv[1]
traj_loc = sys.argv[2]
start = int(sys.argv[3])
end = int(sys.argv[4])
ref_pdb = sys.argv[5]
system = sys.argv[6]

alignment = 'protein and name CA and (resid 20:25 or resid 50:55 or resid 73:75 or resid 90:94 or resid 112:116 or resid 142:147 or resid 165:169 or resid 190:194 or resid 214:218 or resid 236:240 or resid 253:258 or resid 303:307)'

nSel = len(sel_list.sel)

flush = sys.stdout.flush

# ----------------------------------------
# SUBROUTINES:

def ffprint(string):
	print '%s' %(string)
        flush()

# ----------------------------------------
# MAIN PROGRAM:

ref = MDAnalysis.Universe('%s' %(ref_pdb))
ref_align = ref.select_atoms(alignment)
ref_all = ref.select_atoms('all')
ref_backbone = ref.select_atoms('backbone')
ref_all.translate(-ref_backbone.center_of_mass())
pos0 = ref_align.coordinates()

u = MDAnalysis.Universe('%s' %(pdb))
u_align = u.select_atoms(alignment)
u_all = u.select_atoms('all')
u_backbone = u.select_atoms('backbone')

if len(u_align) != len(ref_align):
	ffprint('Alignment atom selections do not have the same number of atoms.')
	sys.exit()

rest = u.select_atoms('not (resname WAT or resname Na+ or resname Cl- or protein)')
num_res = len(rest.residues)

# make selections to compute RMSD for
u_sel = ['']*nSel
ref_sel = ['']*nSel
for i in range(nSel):
	selection = sel_list.sel[i][1]

	ref_temp = ref.select_atoms(selection)
	u_sel[i] = u.select_atoms(selection)

	if len(u_sel[i]) != len(ref_temp):
		ffprint('Number of atoms do not match for selection %s' %(i))
		sys.exit()

	ref_sel[i] = ref_temp.coordinates()

# open output files
out1 = open('%s.rmsd.dat' %(system), 'w')
nSteps = 0

ffprint('Beginning trajectory analysis')
while start <= end:
	ffprint('Loading trajectory %s' %(start))
	u.load_new('%sproduction.%s/production.%s.dcd' %(traj_loc,start,start))
	nSteps += len(u.trajectory)
	
	# Loop through trajectory
	for ts in u.trajectory:
		# obtain dimension values to be used for unwrapping atoms
		dimensions = u.dimensions[:3]
	
		# Align to reference (moves COM of backbone to origin)
		u_all.translate(-u_backbone.center_of_mass())
		# Calculate the rotational matrix to align u to the ref
		R, rmsd = rotation_matrix(u_align.coordinates(), pos0)
		# Apply rotation matrix to atoms within u
		u_all.rotate(R)

		# Fix the wrapping issues
		for i in range(num_res):
			COM = np.zeros(3)
			# Calculate the COM of residues;
			COM = rest.residues[i].center_of_mass()
			
			# CALCULATING AND APPLYING THE TRANSLATIONAL MATRIX TO RESIDUE i
			t = wrapping(COM,dimensions)
			rest.residues[i].atoms.translate(t)

		# loop through selections and compute RMSD
		for i in range(nSel):
			u_coords = u_sel[i].coordinates()
			ref_coords = ref_sel[i]
			
			rmsd = RMSD(u_coords,ref_coords)
			out1.write('%10.6f   ' %(rmsd))

		out1.write('\n')
	start += 1

out1.close()
print 'Analyzed %d steps' %(nSteps)

