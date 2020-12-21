### This is the PH function called to solve instances of the pharmaceutical R&D pipeline problem
### The file requires info from the bounding approach

import os
import sys
from pyomo.environ import *
from pyomo.pysp.scenariotree.manager import \
    (ScenarioTreeManagerClientPyro,
     InvocationType)
from pyomo.pysp.ef import create_ef_instance
from pyomo.opt import SolverFactory
import pyro4

def PH_Bound_Func(scenario_iterable):
	# declare the number of scenarios over which to construct a simple
	# two-stage scenario tree
	num_scenarios = len(scenario_iterable)


	
#
# Define the scenario tree structure as well as stage
# costs and variables
#
def pysp_scenario_tree_model_callback():
    from pyomo.pysp.scenariotree.tree_structure_model import \
        CreateConcreteTwoStageScenarioTreeModel

    st_model = CreateConcreteTwoStageScenarioTreeModel(num_scenarios)

    first_stage = st_model.Stages.first()
    second_stage = st_model.Stages.last()

    # First Stage
    st_model.StageCost[first_stage] = 'FirstStageCost'
    st_model.StageVariables[first_stage].add('x_first_stage')

    # Second Stage
    st_model.StageCost[second_stage] = 'SecondStageCost'
    st_model.StageVariables[second_stage].add('y_second_stage')

    return st_model
