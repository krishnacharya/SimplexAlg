''' LP Solver '''
__author__ = "Krishna ACHARYA"
__license__ = "GPL"
__version__ = "1.0.1"
__email__ = "krishnacharya97@gmail.com"

import numpy as np
from fractions import Fraction
import math
from bidict import bidict
from printing import *
from helper_funcs import *
import pathlib

def phase2(tableau, basic_in_row, rule = "max", verbose=False):
    status = "bounded"; num_pivots = 0# number of pivots used in phase2
    while(find_increaser(tableau, rule) != -1): # while there exists a non basic variable to increase along
        col = find_increaser(tableau, rule)
        pivot_row, ___ = ratio_test(tableau[:,col], tableau[:,-1], rule, basic_in_row)
        if pivot_row == -1:
            status = "unbounded"
            break
        tableau, basic_in_row =  pivot(tableau, pivot_row, col, basic_in_row, verbose)
        num_pivots += 1
    return tableau, basic_in_row, status, num_pivots

def phase1(tableau, basic_in_row, rule="max", verbose=False):
    num_pivots = 0
    while(find_increaser(tableau, rule) != -1):
        col = find_increaser(tableau, rule)
        pivot_row, ___ = ratio_test(tableau[:,col], tableau[:,-1], rule, basic_in_row)
        tableau, basic_in_row = pivot(tableau, pivot_row, col, basic_in_row, verbose)
        num_pivots += 1
    if tableau[0,-1] != 0:
        return tableau, basic_in_row, "infeasible", num_pivots
    return tableau, basic_in_row, "feasible", num_pivots

def pivot(tableau, pivot_row, col, basic_in_row, verbose = False):
    # should also swap basis and non basis variable
    basis_leaving_variable = basic_in_row[pivot_row]
    basic_in_row[pivot_row] = col # col is the basis entering variable
    tableau[pivot_row] = tableau[pivot_row] / tableau[pivot_row, col]# generates a 1 at (pr,col)
    # can be vectorized?
    for i in range(len(tableau)):
        if i == pivot_row:
            continue
        tableau[i] = tableau[i] - (tableau[i,col]) * tableau[pivot_row]
    if verbose:
        pivot_printing(tableau, basis_leaving_variable, col, basic_in_row)
    return tableau, basic_in_row

def interpret_solution(tableau, basic_in_row, n):
    # get the values for each of the OG variables
    # outputs cordinates and solution value
    soln = np.zeros(n+1, dtype=Fraction);soln[0] = 'dum'# dummy
    for i in range(1,n+1):
        if i in basic_in_row.inverse:
            soln[i] = tableau[basic_in_row.inverse[i],-1]
    return soln, tableau[0,-1]

def phase1_preproc(n, m, cst, updated_cst, ssa_cst):
    tot_extra = np.sum(ssa_cst)
    tot_slsu = np.sum(ssa_cst[:,0:2]) # total number of slack and surplus vars  
    tableau = np.zeros((m+1, n + tot_extra + 2), dtype = Fraction)
    basic_in_row = bidict()
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
    return tableau, basic_in_row

def normalize(n, tableau, basic_in_row):
    for i in range(1, n+1):
        if i in basic_in_row.inverse:
            tableau[0] = tableau[0] - tableau[0, i] * tableau[basic_in_row.inverse[i]]
    return tableau

def remove_degen_afterphase1(tableau, basic_in_row, ssa_cst, n):
    tot_extra = np.sum(ssa_cst)
    tot_slsu = np.sum(ssa_cst[:, 0:2])
    for i in range(n+tot_slsu+1, n+tot_extra+1): # looping over artifical variables
        if i in basic_in_row.inverse:
            pivot_row = basic_in_row.inverse[i]
            col = 0
            for j in range(1, n+tot_slsu+1):
                if tableau[pivot_row, j] != 0:
                    col = j
            if col == 0:
                print("No variable to replace artificial?? strange")
            tableau, basic_in_row = pivot(tableau, pivot_row, col, basic_in_row)
    return tableau, basic_in_row

def phase2_preproc(tableau, basic_in_row, ssa_cst, obj, n):
    tot_extra = np.sum(ssa_cst)
    tot_slsu = np.sum(ssa_cst[:, 0:2])
    tableau[0] = 0 * tableau[0]
    tableau[0, 1:n + 1] = -1 * obj
    tableau = np.delete(tableau, range(n + tot_slsu + 1, n + tot_extra + 1), axis=1)  # must drop all the artificial variables, nonbasic
    tableau = normalize(n, tableau, basic_in_row)
    return tableau
    
def direct_phase2_preproc(n, m, obj, cst, updated_cst, ssa_cst):
    tableau = np.zeros((m+1, n+m+2), dtype=Fraction)
    tableau[0, 1:n+1] = -1 * obj #objective function
    for j in range(1,m+1):
        tableau[j,1:n+1] = updated_cst[j-1]
        tableau[j,n+j] = 1
        tableau[j,-1] = cst[j-1]
    basic_in_row = {j: n+j for j in range(1,m+1)}
    basic_in_row = bidict(basic_in_row)
    return tableau, basic_in_row

def solve(input_file, rule="max", verbose=False):
    '''
        Algorithm: If there are are only slack variable go directly to phase2, else do phase1 then phase2
                   The rest of the code in between is to make the tableau into the correct format and for
                   printing
    '''
    n, m, obj, cst, updated_cst, ssa_cst = parse_input(input_file)
    soln = "notYetDefined"
    opval = "notYetDefined"
    if np.sum(ssa_cst[:,1:]) == 0:
        tableau, basic_in_row = direct_phase2_preproc(n, m, obj, cst, updated_cst, ssa_cst)
        initial_printing(obj, cst, updated_cst, ssa_cst, tableau, basic_in_row, verbose) 
        print("-------------------NO PHASE I--------------------\n\n")
        print("-------------------PHASE II--------------------\n\n")
        tableau, basic_in_row, status, num_pivots = phase2(tableau, basic_in_row, rule, verbose)
        if status == "unbounded":
            print("This linear program is UNBOUNDED")
        else:
            soln, opval = interpret_solution(tableau, basic_in_row, n)
            final_printing(soln, opval, num_pivots, n, rule)    
    else:
        tableau, basic_in_row = phase1_preproc(n, m, cst, updated_cst, ssa_cst)
        initial_printing(obj, cst, updated_cst, ssa_cst, tableau, basic_in_row, verbose)
        print("-------------------PHASE I--------------------\n\n")
        tableau, basic_in_row, status, num_pivots_p1 = phase1(tableau, basic_in_row, rule, verbose)
        if status == "infeasible": #also all the artificial vars must be basic
            print("This linear program is INFEASIBLE")
        else:
            print("This linear program is FEASIBLE\n\n")
            tableau, basic_in_row = remove_degen_afterphase1(tableau, basic_in_row, ssa_cst, n)
            tableau = phase2_preproc(tableau, basic_in_row, ssa_cst, obj, n)
            print("-------------------PHASE II--------------------\n\n")
            tableau, basic_in_row, status, num_pivots_p2 = phase2(tableau, basic_in_row, rule, verbose)
            if status == "unbounded":
                print("This linear program is UNBOUNDED")
            else:   
                soln, opval = interpret_solution(tableau, basic_in_row, n)
                final_printing(soln, opval, num_pivots_p1 + num_pivots_p2, n, rule)
    return status, soln, opval
    
#blockPrint()
# status, soln, opval = solve("degen-ex1.txt", verbose=False)
#enablePrint()
# print(status)
# print(soln)
# print(opval)
# solve("knapsack.txt", rule="max", verbose=True)