class NAC_Graph:
	def __init__(self, graph_dict = {}):
		### Initializes Graph Object
		self._graph_data = graph_dict
		
	def verticies(self):
		return list(self._graph_data.keys())
		
	def edges(self):
		# Returns the edges of the graph
		return self._get_edges()
		
	def add_vertex(self, vertex):
		### Add Vertex to graph
		if vertex not in self._graph_data:
			self._graph_data[vertex] = []	
		
	def add_edge(self, edge):
		### Add Edge to graph
		edge = set(edge)
		(vertex1, vertex2) = tuple(edge)
		if vertex1 in self._graph_data:
			self._graph_data[vertex1].append(vertex2)
		else:
			self._graph_data[vertex1] = [vertex2]		
		
	def _get_edges(self):
		### Find Edges of Graph
		edges = []
		for vertex in self._graph_data:
			for neighbour in self._graph_data[vertex]:
				if {neighbour, vertex} not in edges:
					edges.append({vertex, neighbour})
		return edges
        
        
	def __str__(sellf):
		res = "vertices: "
		for k in self._graph_data:
			res += str(k) + " "
			res += "\nedges: "
		for edge in self.__generate_edges():
			res += str(edge) + " "
		return res
		
	def plot(self):
		return print("Currently Not Working")
	
	def find_path(self,start_vertex,end_vertex, path=[]):
		### Find a path from the start vertex to the end vertex
		graph = self._graph_data
		path = path + [start_vertex]
		if start_vertex == end_vertex:
			return path
		if start_vertex not in graph:
			return None
		for vertex in graph[start_vertex]:
			if vertex not in path:
				extended_path = self.find_path(vertex,end_vertex,path)
				if extended_path:
					return extended_path
		return None
		
	def all_edges(self):
		""" Adds all the edges between all samples"""
		for i in self._graph_data:
			for j in self._graph_data:
				if i != j and j > i:
					self.add_edge((i,j))
					
		

		
		
		
		
		
        
        
  
