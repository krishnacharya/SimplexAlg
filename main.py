#! /usr/bin/env python3
from argparse import ArgumentParser
from simp import *
import time

if __name__ == "__main__":
    rules = {"max", "bland", "random", "step_max"}
    parser = ArgumentParser(description="A simplex solver")
    parser.add_argument("input", help="the simplex input file")
    parser.add_argument("-v", "--verbose", action="store_true",help ="increase output verbosity", default = False)
    parser.add_argument("-r", "--rule", choices=rules, default="max",
                        help="specify the pivot rule to be used")
    arguments = parser.parse_args()
    
    start_time = time.time()
    solve(arguments.input, rule=arguments.rule, verbose=arguments.verbose)
    print("--- %s seconds ---" % (time.time() - start_time))