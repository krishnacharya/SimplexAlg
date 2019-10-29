import numpy as np
from fractions import Fraction
import math
from bidict import bidict
## Basic Tableau generator
# Basic Simplex
# Get objective function in OG Variables(x1,x2,...xn)
# Create Slack Variables for each constraint
# Find improver in obj
# Ratio test, Pivot
# End when no more -ve variables in tableau
# Or if there is a -ve variable and the ratio test gives no bound, unbounded.

###TO IMPLEMENT - dictionary based basic_in_rows, will be faster for the verbose case

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


def find_increaser(vector, rule = "max"):
	column = 0
	if rule == "max":
		column = max_coe(vector) # gives which x_c to choose.
	elif rule == "random": # randomly pick among the -ve elements
		column = random_sel(vector)
	elif rule == "bland":
		column = bland_choice(vector)
	else:
		print("NO OTHER RULES")# possibly throw an error
	if column == -1:
		return -1 # NO MORE INCREASE
	# if vector[column] >= 0:
	#   return -1 #no -ve coeff found, therefore no more increase possible  
	return column

def max_coe(vector):
	# return np.argmin(vector[1:-1]) + 1 #gives which x_c to pivot on
	col = np.argmin(vector[:-1])
	if vector[col] >= 0:
		return -1
	return col

def bland_choice(vector):
	negval_indices = np.where(vector[:-1] < 0)[0]
	if np.array_equal(negval_indices, []):
		return -1
	return negval_indices[0] # among all candidates for entering variable in BI, choose the one with smallest index

def random_sel(vector):
	negval_indices = np.where(vector[:-1] < 0)[0]# a[0] is always 0, dummy val
	if np.array_equal(negval_indices, []):
		return -1
	index = np.random.choice(len(negval_indices), 1) # pick one index
	return negval_indices[index]

def ratio_test_phase2(v1, v2):
	# YET TO HANDLE, pick row witht the lowest numbered variable in the basis
	# ratios = v2[1:] / v1[1:]
	v1 = v1[1:]; v2 = v2[1:]
	# can be optimized
	min_posval = math.inf
	flag = 0 # flag for existence of upper bound till where the variable can be increased
	for i in range(len(v1)):
		lhs = v1[i]; rhs = v2[i]
		if lhs <= 0:
			continue # that variable doesnt appear in this constraint pass, or gives a useless constraint

		ratio = rhs / lhs
		if ratio < 0:
			print("BIZARRE!!")
			continue

		if ratio < min_posval: # if we reach here it implies the ratio is >=0, yet to implement blands condition
			min_posval = ratio
			min_index = i
			flag = 1
		# if rhs/lhs >= 0:
		#   if ratios[i] < min_posval:
		#       min_posval = ratios[i]
		#       min_index = i
		#       flag = 1

	if flag == 0: # no x_c < bound found
		return -1
	return min_index + 1


def phase2(tableau, basic_in_row, rule = "max"):
	while(find_increaser(tableau[0], rule) != -1): # while there exists a non basic variable to increase along
		col = find_increaser(tableau[0], rule)
		pivot_row = ratio_test_phase2(tableau[:,col], tableau[:,-1])
		if pivot_row == -1:
			return "UB"
		tableau, basic_in_row =  pivot(tableau, pivot_row, col, basic_in_row)        
	return tableau[0,-1] #maximum value

def pivot(tableau, pivot_row, col, basic_in_row):
	# should also swap basis and non basis variable
	basic_in_row[pivot_row] = col
	tableau[pivot_row] = tableau[pivot_row] / tableau[pivot_row, col]# generates a 1 at (pr,col)
	# can be optimized
	for i in range(len(tableau)):
		if i == pivot_row:
			continue
		tableau[i] = tableau[i] - (tableau[i,col]) * tableau[pivot_row]
	return tableau, basic_in_row

def phase1(tableau, basic_in_row, rule="max"):
	while(find_increaser(tableau[0], rule) != -1):
		col = find_increaser(tableau[0], rule)
		pivot_row = ratio_test_phase2(tableau[:,col], tableau[:,-1])
		tableau, basic_in_row = pivot(tableau, pivot_row, col, basic_in_row)
	if tableau[0,-1] != 0:
		return tableau, basic_in_row, "infeasible"
	return tableau, basic_in_row, "feasible"

def normalize(n, tableau, basic_in_row):
	#normalization after phase1 gives feasible
	# this could case the issue
	# for i in range(1,n+1):
	# 	row_pos = np.where(basic_in_row == i)[0]
	# 	if not np.array_equal(row_pos,[]):
	# 		tableau[0] = tableau[0] - tableau[0, i] * tableau[row_pos] ##
	# return tableau 
	for i in range(1, n+1):
		if i in basic_in_row.inverse:
			tableau[0] = tableau[0] - tableau[0, i] * tableau[basic_in_row.inverse[i]]
	return tableau

def interpret_solution(tableau, basic_in_row):
	pass


# TO IMPLEMENT: also keep a variable to row mapping like x_1: rowc etc
def solve(input_file, rule="max"):
	n, m, obj, cst, updated_cst, ssa_cst = parse_input(input_file)
	# print(n,m)
	# print(obj, type(obj))
	# print(cst, type(cst))
	# print(updated_cst, type(updated_cst))
	# print(ssa_cst, type(ssa_cst))

	if np.sum(ssa_cst[:,1:]) == 0:
	  # print("Yeah boii")

	  tableau_direct = np.zeros((m+1, n+m+2), dtype=Fraction)
	  # nb_vars = #all the OG variables
	  # b_vars = # all the slack
	  tableau_direct[0,1:n+1] = -1 * obj#objective function
	  for j in range(1,m+1):
		  tableau_direct[j,1:n+1] = updated_cst[j-1]
		  tableau_direct[j,n+j] = 1
		  tableau_direct[j,-1] = cst[j-1]
	  basic_in_row = {j: n+j for j in range(1,m+1)}; basic_in_row[0] = -2
	  basic_in_row = bidict(basic_in_row)
	  #basic_in_row  = np.array(np.append([-2],[n+j for j in range(1,m+1)]),dtype=int) # basic variable associated with particular row
	  # print(tableau_direct)
	  # print(basic_in_row)
	  Opval = phase2(tableau_direct, basic_in_row, rule)
	  if Opval == "UB":
	  	print("UNBOUNDED")
	  else:
	  	#interpret_solution(tableau_direct, basic_in_row)
	  	print(Opval)
	
	else:
		tot_extra = np.sum(ssa_cst)
		tot_slsu = np.sum(ssa_cst[:,0:2]) # total number of slack and surplus vars  
		tableau = np.zeros((m+1, n + tot_extra + 2), dtype = Fraction)
		basic_in_row = bidict({0:-2})
		#basic_in_row = np.zeros(m+1,dtype=int); basic_in_row[0] = -2 # some dummy value
		tableau[0, n+tot_slsu+1 : -1] = np.ones(tot_extra - tot_slsu, dtype=Fraction) # artifical vars
		k = 0; l = 0 # shifters for the slack cols and artficial cols
		for j in range(1, m+1):
			tableau[j,1:n+1] = updated_cst[j-1]
			tableau[j,-1] = cst[j-1]
			if np.array_equal(ssa_cst[j-1], [1,0,0]):
				tableau[j, n+k+1] =  1
				basic_in_row[j] = n+k+1
				k+=1
			elif np.array_equal(ssa_cst[j-1], [0,1,1]):
				tableau[j, n+k+1] =  -1
				tableau[j, n+tot_slsu+l+1] = 1
				basic_in_row[j] = n+tot_slsu+l+1
				tableau[0] = tableau[0] - tableau[j]  # removing the artifical var from obj fun (only non basic in obj func)
				k+=1
				l+=1
			else: # the case [0,0,1]
				tableau[j, n+tot_slsu+l+1] = 1
				basic_in_row[j] = n+tot_slsu+l+1
				tableau[0] = tableau[0] - tableau[j] # removing the artifical var from obj fun (only non basic in obj func)
				l+=1
		tableau, basic_in_row, status = phase1(tableau, basic_in_row, rule)
		if status == "infeasible": #also all the artificial vars must be basic
			print("INFEASIBLE")
		else:
			tableau[0] = 0 * tableau[0] 
			tableau[0, 1:n+1] = -1 * obj
			tableau = np.delete(tableau,range(n+tot_slsu+1, n+tot_extra+1), axis=1) # must drop all the artificial variables
			tableau = normalize(n, tableau, basic_in_row)
			Opval = phase2(tableau, basic_in_row, rule)
			if Opval == "UB":
				print("UNBOUNDED")
			else:
				print(Opval)

solve('knapsack.txt', rule="bland")
