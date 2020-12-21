import sys
import os
import numpy as np
import pdb

class Kruskal_MST:
	def __init__(self,graph,vert_inf,exist_edge=[]):
		self.graph = graph
		self.edges = edge_weight_add(vert_inf)
		self._parent = {}
		self._rank = {}
		self.mst = MST_Kruskal(exist_edge)
	
	def __repr__(self):
		return repr(self.mst)
		
		
	def edge_weight_add(self, vert_inf):
		### Find Edge Weights
		edge_weights = []
		
		edges = graph.edges()
		for i in edges:
			k = tuple(i)
			w = np.linalg.norm(np.array(vert_inf[k[0]])-np.array(vert_inf[k[1]]))
			edge_weights.append((w,) + k)
		return edge_weights
	
	def MST_Kruskal(self,exist_edge):
		""" Implements Kruskal's algorithm for MST """
		
		### Make sets for each vertex
		for verts in self.graph.verticies:
			self.make_set(vert)
			
		### Define Empty MST
		mst = set()
		
		self.edges.sort()
		for e in self.edges:
			weight, vert1, vert2 = e
			if find(vert1) != find(vert2):
				union(vert1,vert2)
				self.mst.add((vert1,vert2))
		return 
		

	def make_set(vert):
		par[vert] = vert
		rank[vert] = 0
	
	def find(vert):
		if self.par[vert] != vert:
			self.par[vert] = find(self.par[vert])
		return self.par[vert]
	
	def uni(vert1,vert2):
		root1 = find(vert1)
		root2 = find(vert2)
		if root1 != root2:
			if self.rank[root1] > self.rank[root2]:
				self.par[root2] = root1
			else:
				self.par[root1] = root2
				if self.rank[root1] == self.rank[root2]:self.rank[root2] += 1
	
