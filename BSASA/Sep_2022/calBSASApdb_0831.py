# This script is meant to calculate the BSASA in pymol through python code 


from pymol import cmd
import os
import csv

# ============================================
# --- LIST OF ALL FUNCTIONS NEEDED ---
# ============================================

# ============================================
# --- Function for calcualting the Kd values --- 
# ============================================
def calcKDs(pdbID):
	# this function will calcualte the Kd values based on stacked trimer topology
	# --- turn pdbID into a string --- 
	pdbName_pre = str(pdbID); pdbName = pdbName_pre[0:4]
	# --- create objects for each subunit ---
	# all structures considered are stacked trimers so 6 subunits total
	# ex: alpha1, beta1 and alpha1,beta1 pair of subunits
	# - top and bottom rings -
	cmd.create("a1", pdbName + " and chain A"); cmd.create("b1",pdbName + " and chain B"); cmd.create("c1",pdbName + " and chain C")
	cmd.create("d1",pdbName + " and chain D"); cmd.create("e1",pdbName + " and chain E"); cmd.create("f1", pdbName + " and chain F")
	# --- create the dimer objects --- 
	# sub1+sub2 complex and all other combinations needed
	# since we know this is a stacked trimer, there are 9 dimer combos
	# - top ring -
	cmd.create("ab1", pdbName + " and chain A+B"); cmd.create("ac1", pdbName + " and chain A+C"); cmd.create("bc1", pdbName + " and chain B+C"); 
	# - bottom ring -
	cmd.create("de1", pdbName + " and chain D+E"); cmd.create("df1", pdbName + " and chain D+F"); cmd.create("ef1", pdbName + " and chain F+E");
	# - between rings- 
	cmd.create("af1", pdbName + " and chain A+F"); cmd.create("be1", pdbName + " and chain E+B"); cmd.create("cd1", pdbName + " and chain C+D")

	# --- get hydrogens onto everything ---
	# (NOTE: must have valid valences on e.g. small organic molecules)
	cmd.h_add()
	
	# # --- make sure all atoms including HETATM within an object occlude one another ---
	#  but ignore solvent
	cmd.flag("ignore", "none")
	cmd.flag("ignore", "solvent")
	# --- use solvent-accessible surface with high sampling density ---
	cmd.set("dot_solvent", 1); cmd.set("dot_density", 3)

	# --- measure the components individually storing the results for later ---
	# get area for each subunit 
	a1_area=cmd.get_area("a1"); b1_area=cmd.get_area("b1"); c1_area=cmd.get_area("c1")
	d1_area=cmd.get_area("d1"); e1_area=cmd.get_area("e1"); f1_area=cmd.get_area("f1")
	
	# --- measure the dimer areas ---
	# - top ring - 
	ab1_area=cmd.get_area("ab1"); bc1_area=cmd.get_area("bc1"); ac1_area=cmd.get_area("ac1")
	# - bottom ring - 
	df1_area=cmd.get_area("df1"); de1_area=cmd.get_area("de1"); ef1_area=cmd.get_area("ef1")
	# - between rings -
	af1_area=cmd.get_area("af1"); be1_area=cmd.get_area("be1"); cd1_area=cmd.get_area("cd1")

	# --- calculate all buried surface areas --- 
	# - top ring - 
	bsasa_ab = (a1_area + b1_area) - ab1_area; bsasa_bc = (b1_area + c1_area) - bc1_area; bsasa_ac = (a1_area + c1_area) - ac1_area
	# - bottom ring - 
	bsasa_de = (d1_area + e1_area) - de1_area; bsasa_df = (d1_area + f1_area) - df1_area; bsasa_ef = (e1_area + f1_area) - ef1_area
	# - between ring - 
	bsasa_af = (a1_area + f1_area) - af1_area; bsasa_be = (b1_area + e1_area) - be1_area; bsasa_cd = (c1_area + d1_area) - cd1_area

	# --- ouput the Kd1 and Kd2 ---
	# Kd1 = mean within all rings 
	Kd1 = (bsasa_ab+bsasa_bc + bsasa_ac+bsasa_de+bsasa_df+bsasa_ef)/6
	# Kd2 = mean between rings 
	Kd2 = (bsasa_af+bsasa_be+bsasa_cd )/3
	
	return Kd1, Kd2

# ============================================
# --- Make array of data --- 
# ============================================
def get_Kds(pdbs):
	# initializations 
	arrayKDs = [] # array that will hold [PDB_ID, Kd1, Kd2]
	# --- load the pdb file --- 
	for fil in pdbs: # for each file, calc the Kds
		cmd.load(fil)
		# --- calcualte Kd1 & Kd2 --- 
		kds = calcKDs(fil)
		# --- store values in array for output ---
		pdbID = fil; pdbName_pre = str(pdbID); pdbName = pdbName_pre[0:4]
		addArrayVal = [pdbName, kds[0], kds[1]]
		arrayKDs.append(addArrayVal)
		#print(arrayKDs)
		#print("Array to add: ", addArrayVal)

	# --- export the array into a .txt file --- 
	with open('PDB_KDcalcs_082022.csv', 'w') as f:
    		csv.writer(f, delimiter=',').writerows(arrayKDs)   


# =============================================
# ---- Lines to run functions ---
# ============================================
# list of all PDB ID's
pdbs = ["4lkb.pdb","2wym.pdb","7fjl.pdb","1li1.pdb","5lq3.pdb"]
get_Kds(pdbs)
# get_Kds(["4lkb.pdb","2wym.pdb"])

