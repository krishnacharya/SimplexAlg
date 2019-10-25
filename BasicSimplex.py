import numpy as np
from fractions import Fraction
## Basic Tableau generator
# Basic Simplex
# Get objective function in OG Variables(x1,x2,...xn)
# Create Slack Variables for each constraint
# Find improver in obj
# Ratio test, Pivot
# End when no more -ve variables in tableau
# Or if there is a -ve variable and the ratio test gives no bound, unbounded.

def get_extra_vars(vec, ineq, b):
	# input is numpy vector (contraints) the ineq and the constraint rhs
	slack = surplus = artficial = 0
	if b >= 0:
		if ineq == '<':
			slack += 1
		elif ineq == '>':
			surplus += 1
			artficial += 1
		else:
			artficial += 1
		return vec, [slack,surplus,artficial]
	else:
		if ineq == '<':
			surplus += 1
			artficial += 1
		elif ineq == '>':
			slack += 1
		else:
			artficial += 1
		return -1*vec, [slack,surplus,artficial]
	# the output will always be a vec and [1,0,0] or [0,1,1] or [0,0,1]

def parse_input(input_file):
	# reading complete file is not expensive, can read ~100 MB within a second
	with open(input_file) as f:
		content = f.read()
	lines = content.split('\n')
	n = int(lines[0])
	m = int(lines[1])
	obj = np.array([Fraction(f) for f in lines[2].split(' ')])
	cst = np.array([Fraction(f) for f in lines[3].split(' ')])
	updated_cst = np.zeros((m, n), dtype=Fraction) # stores the constraint vector LHS after RHS has been made positive
	ssa_cst = np.zeros((m, 3), dtype=int) #stores number of slack, surplus,artificial for each constraint 
	# tsl,tsu,ta = 0,0,0 # total slack, surplus and artificial variables
	# print(obj, cst)
	j = 0 # variable to loop over constraints
	for i in range(4, len(lines)):
		temp = lines[i].split(' ')
		c = np.array([Fraction(f) for f in temp[0:-1]])
		#print(c, temp[-1]) 
		v, triple = get_extra_vars(c, temp[-1], cst[j])
		updated_cst[j] = v
		ssa_cst[j] = triple
		# tsl += triple[0]; tsu +=triple[1]; ta += triple[2]
		j += 1
	# print("Printing updated constraints\n",updated_cst)
	# print("Printing number of extra\n", ssa_cst)
	
	#note now cst must be positive
	return n, m, obj, abs(cst), updated_cst, ssa_cst

def phase1():
	pass

def phase2():
	pass

def pivot():
	pass
	# should also swap basis and non basis variable

def max_coe():
	pass
def bland_choice():
	pass
def random_choice():
	pass

def solve(input_file):
	n, m, obj, abs(cst), updated_cst, ssa_cst = parse_input(input_file)
	if np.sum(ssa_cst[:,1:]) == 0:
		#Go directly to Phase2, as there are no surplus or artificial variables
	else:
		#Go to Phase1
		#Process
		#Go to Phase2


solve('ip2.txt')
