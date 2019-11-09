import numpy as np
from fractions import Fraction
import math
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

def parse_input(input_file):
    # reading complete file is not expensive, can read ~100 MB within a second
    with open(input_file) as f:
        content = f.read()
    lines = content.strip().split('\n')
    for l in range(len(lines)):
        lines[l] = lines[l].strip() # removing trailing and leading on each line
    n = int(lines[0])
    m = int(lines[1])
    obj = np.array([Fraction(f) for f in lines[2].split(' ')])
    cst = np.array([Fraction(f) for f in lines[3].split(' ')])
    updated_cst = np.zeros((m, n), dtype=Fraction) # stores the constraint vector LHS after RHS has been made positive
    ssa_cst = np.zeros((m, 3), dtype=int) #stores number of slack, surplus,artificial for each constraint 
    j = 0 # variable to loop over constraints
    for i in range(4, len(lines)):
        temp = lines[i].split(' ')
        c = np.array([Fraction(f) for f in temp[0:-1]])
        v, triple = get_extra_vars(c, temp[-1], cst[j])
        updated_cst[j] = v
        ssa_cst[j] = triple
        j += 1
    return n, m, obj, abs(cst), updated_cst, ssa_cst

def find_increaser(tableau, rule = "max"):
    column = 0
    if rule == "max":
        column = max_coe(tableau[0]) # gives which x_c to choose.
    elif rule == "random": # randomly pick among the -ve elements
        column = random_sel(tableau[0])
    elif rule == "bland":
        column = bland_choice(tableau[0])
    elif rule == "step_max":
        column = step_max(tableau)
    else:
        print("NO OTHER RULES")# possibly throw an error
    if column == -1:
        return -1 # NO MORE INCREASE 
    return column

def max_coe(vector):
    # return np.argmin(vector[1:-1]) + 1 #gives which x_c to pivot on
    col = np.argmin(vector[:-1])
    if vector[col] >= 0:
        return -1
    return col

def step_max(tableau):
    vector = tableau[0]
    negval_indices = np.where(vector[:-1] < 0)[0]
    if np.array_equal(negval_indices, []):
        return -1
    neg_minprod = math.inf
    col_index = np.argmin(vector[:-1])
    for ind in negval_indices:
        min_index, min_posval = ratio_test(tableau[:,ind], tableau[:,-1], "step_max",[])
        if min_posval == -1: # no x_c < bound found
            return ind
        if(tableau[0,ind] * min_posval < neg_minprod):
            neg_minprod = tableau[0,ind] * min_posval
            col_index = ind
    return col_index

def bland_choice(vector):
    negval_indices = np.where(vector[:-1] < 0)[0]
    if np.array_equal(negval_indices, []):
        return -1
    return negval_indices[0] # among all candidates for entering variable in BI, choose the one with smallest index

def random_sel(vector):
    negval_indices = np.where(vector[:-1] < 0)[0]# a[0] is always 0, dummy val
    if np.array_equal(negval_indices, []):
        return -1
    index = np.random.choice(len(negval_indices), 1)[0] # pick one index
    return negval_indices[index]

def ratio_test(v1, v2, rule, basic_in_row):
    min_posval = math.inf
    flag = 0 # flag for existence of upper bound till where the variable can be increased
    for i in range(1, len(v1)):
        lhs = v1[i]; rhs = v2[i]
        if lhs <= 0:
            continue # that variable doesnt appear in this constraint pass, or gives a useless constraint

        ratio = rhs / lhs
        if ratio < 0:
            print("BIZARRE!!")
            continue

        if rule == "bland" and ratio == min_posval and basic_in_row[i] < basic_in_row[min_index]:
            min_index = i
            flag = 1

        if ratio < min_posval: # if its not blands rule just pick the first guy thats less
            min_posval = ratio
            min_index = i
            flag = 1
    if flag == 0: # no x_c < bound found
        return -1,-1
    return min_index, min_posval