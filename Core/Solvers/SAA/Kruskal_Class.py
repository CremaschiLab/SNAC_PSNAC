import sys
import os
import numpy as np
import pdb

class Kruskal_MST:
	def __init__(self,graph,vert_inf,exist_edge=[]):
		self.graph = graph
		self.edges = self.edge_weight_add(vert_inf)
		self._parent = {}
		self._rank = {}
		self.mst = set() 
		self.MST_Kruskal(exist_edge)
	
	def __repr__(self):
		return repr(self.mst)
		
		
	def edge_weight_add(self, vert_inf):
		### Find Edge Weights
		edge_weights = []
		
		edges = self.graph.edges()
		for i in edges:
			k = tuple(i)
			w = np.linalg.norm(np.array(vert_inf[k[0]])-np.array(vert_inf[k[1]]))
			edge_weights.append((w,) + k)
		return edge_weights
	
	def MST_Kruskal(self,exist_edge):
		""" Implements Kruskal's algorithm for MST """
		
		### Make sets for each vertex
		for vert in self.graph.verticies():
			self.make_set(vert)
			
		### Define Empty MST
		mst = set()
		
		### Add existing edges
		for e in exist_edge:
			vert1,vert2 = e
			### Add edge
			self._union(vert1,vert2)
			
		
		
		### Add edges to form MST
		self.edges.sort()
		for e in self.edges:
			weight, vert1, vert2 = e
			if self.find(vert1) != self.find(vert2):
				self._union(vert1,vert2)
				self.mst.add((vert1,vert2))
		return 
		

	def make_set(self,vert):
		self._parent[vert] = vert
		self._rank[vert] = 0
	
	def find(self,vert):
		if self._parent[vert] != vert:
			self._parent[vert] = self.find(self._parent[vert])
		return self._parent[vert]
	
	def _union(self,vert1,vert2):
		root1 = self.find(vert1)
		root2 = self.find(vert2)
		if root1 != root2:
			if self._rank[root1] > self._rank[root2]:
				self._parent[root2] = root1
			else:
				self._parent[root1] = root2
				if self._rank[root1] == self._rank[root2]:self._rank[root2] += 1
	
