from .NAC_Graph import NAC_Graph
from .Kruskal_Class import Kruskal_MST as MST
import pdb

def SubsetNAC(row_vec, dim, Scenarios,Scenario_Realizations,mst):
	### Does it need to be cut?
	if len(row_vec) == 1:
		### Is it a MSP tree?
		subset = set()
		for s in Scenarios:
			if Scenario_Realizations[s][dim] == row_vec[0]:
				subset.add(s)
				
		exedge_sub = []
		for (s,sp) in mst:
			if s in subset and sp in subset:
				exedge_sub.append((s,sp))
				
		if len(exedge_sub) < len(subset) - 1:
			graph_sub = NAC_Graph({})
			
			for s in subset:
				graph_sub.add_vertex(s)
				graph_sub.all_edges()
				
			### Find MST
			MST_sub = MST(graph_sub, Scenario_Realizations, exedge_sub)
				
			### Update Overall Graph
			mst = mst.union(MST_sub.mst)
			
				
		
	else:
		### Define Cut
		for i in row_vec[:-1]:
			### Define Scenario Subsets
			subset1 = set()
			subset2 = set()
			
			### Given Cut Calculate Subsets
			for s in Scenarios:
				if Scenario_Realizations[s][dim] < i + 1:
					subset1.add(s)
				else:
					subset2.add(s)
		
			### Determine existing edges for each subset
			exedge_sub1 = []
			exedge_sub2 = []
	
			for (s,sp) in mst:
				if s in subset1 and sp in subset1:
					exedge_sub1.append((s,sp))
				elif s in subset2 and sp in subset2:
					exedge_sub2.append((s,sp))

			if len(exedge_sub1) < len(subset1) - 1:
				### Create Graph for subset1
				graph_sub1 = NAC_Graph({})
				
				for s in subset1:
					graph_sub1.add_vertex(s)
				graph_sub1.all_edges()
				
				### Find MST
				MST_sub1 = MST(graph_sub1, Scenario_Realizations, exedge_sub1)
				
				### Update Overall Graph
				mst = mst.union(MST_sub1.mst)
				
			if len(exedge_sub2) < len(subset2) - 1:
				### Create Graph for subset2
				graph_sub2 = NAC_Graph({})
		
				for s in subset2:
					graph_sub2.add_vertex(s)
				graph_sub2.all_edges()
				
				### Find MST		
				MST_sub2 = MST(graph_sub2, Scenario_Realizations, exedge_sub2)
				
				### Update Overall Graph
				mst = mst.union(MST_sub2.mst)
			
	return mst
		
	

def NAC_Generator(Scenarios, Scenario_Realizations, Outcomes):
	"""
		Scenarios is a list of scenario numbers i.e. [1,2,3,...S]
		
		Scenario_Realizations is a dictionary mapping S to its outcome 
			i.e. { 1:(1,0,1)}, notice the realization of uncertain 
			parameters is a vector. The vector length should be equal to the 
			number of uncertain parameters
			
		Outcomes is a list corresponding to the number of realizations for
			each uncertain parameter, in the same order as the Scenario_Realizations
			vectors!!! It should have a length equal to the number of uncertain
			parameters
	
	"""
	### Create undirected graph with scenarios connecting all verticies
	graph = NAC_Graph()

	for s in Scenarios:
		graph.add_vertex(s)
	
	graph.all_edges() 
	MST_Graph = MST(graph,Scenario_Realizations)
	mst = MST_Graph.mst
	
	### For all dimensions
	for i in range(len(Outcomes)):
		rows = range(Outcomes[i] + 1)
		while len(rows) > 0:
			mst = SubsetNAC(rows,i,Scenarios,Scenario_Realizations,mst)
			try:
				rows = rows[1:-1]
			except:
				rows = []
	
	
	return mst


				
