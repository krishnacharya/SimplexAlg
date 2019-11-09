import numpy as np
import sys, os
def pivot_printing(tableau, basis_leaving_variable, col, basic_in_row):
    print("The entering variable is x_"+str(col))
    print("The leaving variable is x_"+str(basis_leaving_variable))
    print("\n")
    br_inv = basic_in_row.inverse
    for basic_var in sorted(br_inv):
        row = tableau[br_inv[basic_var]]
        print("x_" + str(basic_var)+" = "+str(row[-1]), end = " ")     
        for i in range(1, len(tableau[0])-1):
            if i == basic_var:
                continue
            if row[i] > 0:
                print(str(-1*row[i])+"x_"+str(i),end=" ")
            elif row[i] < 0:
                print("+" + str(-1*row[i])+"x_"+str(i),end=" ")
        print()
    print("-------------------------------------------------")
    print("z = "+ str(tableau[0,-1]),end = " ")
    for i in range(1,len(tableau[0])-1):
        if tableau[0,i] > 0:
            print(str(-1*tableau[0, i])+"x_"+str(i),end=" ")
        elif tableau[0,i] < 0:
            print("+" + str(-1*tableau[0, i])+"x_"+str(i),end=" ")
    print("\n\n")    

def initial_printing(obj, cst, updated_cst, ssa_cst, tableau, basic_in_row, verbose=False):
    #Lot of boiler plate printing, no serious algorithms here.
    print("The input linear program is:\n")
    print("Maximize", end="  ")
    for i in range(len(obj)):
        if obj[i] < 0:
            print(str(obj[i])+"x_"+str(i+1), end=" ")
        elif obj[i] > 0:
            print("+"+str(obj[i])+"x_"+str(i+1), end=" ")
    print()
    print("Such that")
    for i in range(len(cst)):
        allzeros_flag = 0
        for j in range(len(obj)):
            if updated_cst[i,j] < 0:
                print(str(updated_cst[i,j])+"x_"+str(j+1), end=" ")
                allzeros_flag = 1
            elif updated_cst[i,j] > 0:
                print("+"+str(updated_cst[i,j])+"x_"+str(j+1), end=" ")
                allzeros_flag = 1
        if allzeros_flag != 0:
            if np.array_equal(ssa_cst[i], [1,0,0]):
                print("<=",end = " ")
            elif np.array_equal(ssa_cst[i], [0,1,1]):
                print(">=",end = " ")
            else:
                print("=", end=" ")
            print(cst[i])
    for i in range(len(obj)):
        print("x_"+str(i+1),end=" ; ")
    print("are non negative\n\n")
    if verbose:
        print("The initial dictionary is:\n")
        br_inv = basic_in_row.inverse
        for basic_var in sorted(br_inv):
            row = tableau[br_inv[basic_var]]
            print("x_" + str(basic_var)+" = "+str(row[-1]), end = " ")     
            for i in range(1, len(tableau[0])-1):
                if i == basic_var:
                    continue
                if row[i] > 0:
                    print(str(-1*row[i])+"x_"+str(i),end=" ")
                elif row[i] < 0:
                    print("+" + str(-1*row[i])+"x_"+str(i),end=" ")
            print()
        print("-------------------------------------------------")
        print("z = "+ str(tableau[0,-1]),end = " ")
        for i in range(1,len(tableau[0])-1):
            if tableau[0,i] > 0:
                print(str(-1*tableau[0, i])+"x_"+str(i),end=" ")
            elif tableau[0,i] < 0:
                print("+" + str(-1*tableau[0, i])+"x_"+str(i),end=" ")
        print("\n\n") 

def final_printing(soln, opval, num_pivots, n, rule):
    print("An optimal solution is:", end=" ")
    for i in range(1,n+1):
        print("x_"+str(i),"=",soln[i],end="; ")
    print()
    print("The value of the objective for this solution is: ", opval)
    print("The number of pivoting operations is: ", num_pivots)
    print("The pivot rule used is: ", rule)

def blockPrint():
    sys.stdout = open(os.devnull, 'w')

def enablePrint():
    sys.stdout = sys.__stdout__