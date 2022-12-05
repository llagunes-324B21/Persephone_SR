# This code is my attempt at downloading and working with XML files/data

# I downloaded a .dat file from PDBePISA that contains ALL the hexamers available on the database. 
# We need to parse out that data and make sure we are only considering stacked rings

# Name: Leo Lagunes and Paige M.
# Date: 09/07/22 

# Last update: 09/07/22

# ===============================================
# --- PACKAGES --- 
# ===============================================
import sys
import numpy as np

# ===============================================
# --- FUNCTIONS --- 
# ===============================================

Filename = 'dbsres.dat' # raw data 

# --- read in data --- 
EntryData = []
EntryWhiteSpaceCount = 0
CurrentEntryNum = 0

with open(Filename, 'r') as File:
		InData = False
		for Line in File:
			print(Line)
			if ("pdb_entry" in Line) and not ("/pdb_entry" in Line):
				print("Enter if statement!")















