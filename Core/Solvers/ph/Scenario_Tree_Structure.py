import sys
import os
from pyomo.environ import *
from pyomo.pysp.scenariotree.tree_structure_model import ScenarioTreeModelFromNetworkX
from pyomo.pysp.phinit import exec_runph as exec_runph
from pyomo.pysp.phinit import main as ph_instance
from Core.Solvers.PH.Decision_Tree_Object import Multistage_PYSP_Object as DT
import networkx
import pdb
from pyomo.pysp.scenariotree.manager_solver import \
    ScenarioTreeManagerSolverClientSerial
import builtins
import Core.Solvers.MTSSP.M2S_item as M2S_item
import copy

#
# Define the scenario tree structure as well as stage
# costs and variables
#
def pysp_scenario_tree_model_callback():
	from pyomo.pysp.scenariotree.tree_structure_model import ScenarioTreeModelFromNetworkX
	
    ## build a decision tree object
	Tree_Graph = DT(dec_GLOBAL,MD_GLOBAL)
	builtins.Tree_Graph_GLOBAL = Tree_Graph
	
	### Initialize scenario tree model
	try:
		stm = ScenarioTreeModelFromNetworkX(Tree_Graph.G,edge_probability_attribute='probability', stage_names = Tree_Graph.stage_names)
	except:
		pdb.set_trace()
		
    ### Define Variables for each Node
	for N in Tree_Graph.G.node:
		if N.startswith('s') or N.startswith('pn'):
			stm.NodeVariables[N].add(str(N)+'VarXX')
			stm.NodeDerivedVariables[N].add(str(N)+'VarY')
			stm.NodeDerivedVariables[N].add(str(N)+'VarZ')
		else:
			###NAC Constrained Decision Variables
			stm.NodeVariables[N].add(str(N)+'VarX')
			stm.NodeVariables[N].add(str(N)+'VarXX')
			stm.NodeDerivedVariables[N].add(str(N)+'VarY')
			stm.NodeDerivedVariables[N].add(str(N)+'VarZ')
		
		
	for stage_name in Tree_Graph.stage_names:	
		stm.StageCost[stage_name] = str(stage_name) + "StageCost"

	return stm	
