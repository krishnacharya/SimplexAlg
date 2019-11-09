

# Simplex algorithm - Krishna Acharya

## Running the code

To run the default algorithm, do `python ./main.py <input file path>`. 
For e.g run `python ./main.py -v -r random <input file path>` for running the simplex algorithm with random pivot rule in verbose mode. 
Run `python ./main.py -h` for additional explanation. 
Note: verbose mode is off by default. Also the output is printed directly on the terminal, if this is undesired just redirect the ouput with the `>`  operator in the shell i.e do `python ./main.py -v -r random <input file path>  > <output file path>`

## Code structure
`main.py` is the main file accepts command line arguments and then calls `simp.py`

`simp.py` has the simplex algorithm, phase I, II, pivoting etc.

`helper_funcs.py` has all the helper functions used in `simp.py` e.g  ratio test, parsing the input, pivoting rules -  max, random, bland and step_max. See these to really understand my logic for pivoting.

`printing.py` has all the functions for printing the output properly. As such it has no "deep logic" in it

## Implementation, Design choice

I have coded in Python3. The `fractions` module is used to represent values in fractional format, [Numpy](http://www.numpy.org/) to maintain the tableau and [bidict](https://bidict.readthedocs.io/en/master/) to represent the mapping from constraint row to a basic variable.

bidict is a bidirectional dictionary data structure, so it is randomized O(1) for `key -> value` and `value -> key` This is useful as sometimes I need the basic variable corresponding to a particular constraint, and sometimes the constraint corresponding to a basic variable. 

For some LP instance a particular pivot rule may be non-terminating. I have left it up to the user to terminate a run directly from the console.

## Pivoting rules

There are four pivoting rules

 * max (default) - maximum coeffiecent rule
 * bland
 * random
 * step_max

## Drawbacks/Limitation

The full power of numpy isn't being harnessed. Since numpy treats `fractions` as objects, it cannot do the floating point optimizations it is well known for. 

Also since python is interpreted it is already at a drawback compared to compiled languages. I have tested on upto 120 variables and 190 constraints. 


## Possible improvement

If I had more time I would like to see if floating point results are far from the exact fractional values. 

Also what the tradeoff for computation time vs accuracy(bits of precision, exact fractional values) is.

## Note on input and output fileformat
Other than for the TD problems I havent prefixed my surname(acharya) before every file name.

The tags random, bland, stepmax  at the end of the output file specify which rule was used while generating it

No tag means the default rule i.e max  

## Acknowledgement for ideas and test cases

 * Ralph for the Bipartite matching examples
 * Gabriel for pointing out how to handle artificial variables in basis after Phase 1
 * Gaetan(aka Gaga) for the random lp examples(orf, osterone, ament) and for our discord discussion on the idea behind stepmax
 * Emile for sorci-comp
