# SNAC_PSNAC
The Sample Non-Anticipativity Constraint (SNAC) algorithm and Parallelization version (PSNAC) for minimum cardinality NAC generation

## Summary
This repository contains the code for The Sample Non-Anticipativity Constraint (SNAC) algorithm and Parallelization version (PSNAC) for minimum cardinality NAC generation and the data to replicate the results of manuscript
  >Zuo Zeng, Brianna Christian, Alexander Vinel, Selen Cremaschi. <br />
  >A Graph Theoretic Approach to Non-Anticipativity Constraint Generation in Multistage Stochastic Programs with Incomplete Scenario Sets. <br />
  >Computers & Operation Research. 

Please note that Core is a work-in-progress. The current version as of December 21, 2020 has been copied into this repository for result reproducibility. An up-to-date version can be found at the public repo https://github.com/CremaschiLab/SNAC_PSNAC

## Quickstart
###1. Clone file to obtain a copy of the source code

###2. How to run algorithm <br />
The command line: python X.py <br />
Example: python SNAC.py 
Please note that you will need to have Pyomo 5.6.6 (http://www.pyomo.org/), Python 3.5.0 and CPLEX 12.6.3. (as the solver) installed before you run it. 

X: speficies the approach that is used to generate NACs. <br />
* Option 1: Full<br />
* Option 2: SNAC
* Option 3: SNAC_alter

## Data
The test data is in Problem Files.
New created data file should be located in Problem Files.

