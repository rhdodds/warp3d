from collections import OrderedDict

import numpy as np

class FEModel(object):
	def __init__(self):
		self.title = ''

class SimpleFEModel(FEModel):
	"""
		A class representing a finite element model with results data.
		Probably want this to be as flexible as possible...
		Want: variable nodes, elements, results at nodes, results at elements,
		integration results at elements, node properties, and element properties.
	"""
	def __init__(self):
		super(SimpleFEModel,self).__init__()
		self.timesteps = 0
		self.times = []
		self.nodes = []
		self.elements = []
		self.nsets = OrderedDict()
		self.eblks = OrderedDict()
		self.efields = OrderedDict()
		self.nfields = OrderedDict()
		self.ifields = OrderedDict()

	def read(self, reader):
		"""
			Read in from a file via the reader (which must quack like a duck by
			having the relevant methods).
		"""
		self.title = reader.title
		self.timesteps = reader.num_steps
		self.dim = reader.dim
		self.times = reader.times
		
		print(" ...reading nodes")
		for coords in reader.node_iterator():
			self.nodes.append(SimpleNode(coords))

		print(" ...reading elements")
		for conn, etype in reader.elem_iterator():
			self.elements.append(SimpleElement([ self.nodes[node] for node in conn ],
				etype))
	
		print(" ...reading node coordinates")
		for name, nodes in reader.nset_iterator():
			self.nsets[name] = [self.nodes[node] for node in nodes]

		print(" ...reading element data" )
		for name, elems in reader.eblk_iterator():
			self.eblks[name] = [self.elements[elem] for elem in elems]
		print(" ...done")	

		print("Reading nodal results...")
		for name, shape, fieldit in reader.node_field_iterator():
			self.nfields[name] = shape
			for i,values in enumerate(fieldit):
				self.nodes[i].fields[name] = values
		print(" ...done")	

		print("Reading element results...")
		for name, shape, fieldit in reader.element_field_iterator():
			self.efields[name] = shape
			for i,values in enumerate(fieldit):
				self.elements[i].fields[name] = values
		print(" ...done")	

		print("Reading integration point results...")
		for name, shape, fieldit in reader.integration_field_iterator():
			self.ifields[name] = shape
			for i,values in enumerate(fieldit):
				raise NotImplementedError("Err, how to do this?")
 
	@property
	def num_nodes(self):
		return len(self.nodes)

	@property
	def num_elems(self):
		return len(self.elements)

	@property
	def num_eblks(self):
		return len(self.eblks)

	@property
	def num_nsets(self):
		return len(self.nsets)
	
	@property
	def num_elem_fields(self):
		return len(self.efields)

	@property
	def num_node_fields(self):
		return len(self.nfields)

	@property
	def num_integration_fields(self):
		return len(self.ifields)

class Node(object):
	def __init__(self):
		pass

class SimpleNode(Node):
	"""
		A node class, containing:
			coordinates
			properties
			fields
	"""
	def __init__(self, coords):
		self.coords = coords
		self.properties = dict()
		self.fields = dict()

class Element(object):
	def __init__(self):
		pass

class SimpleElement(Element):
	"""
		An element class, containing:
			connectivity
			element type (from standard list, eventually)
			properties
			fields
			*integration point fields (do later)
	"""
	def __init__(self, connectivity, etype):
		self.conn = connectivity
		self.etype = etype
		self.properties = dict()
		self.fields = dict()
