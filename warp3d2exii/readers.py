import collections
import itertools
import os.path

import numpy as np
import scipy as sp
import scipy.io as sio

from patran import *

# Print progress every wo packets of patran neutral file
wo = 1000

# Primary reader methods are:
# node_iterator: returns (in file order) the coordinates of each node
# elem_iterator: returns (in file order) the connectivities and type of each
# 	element

class Reader(object):
	def __init__(self):
		self.opened = False

	def __del__(self):
		self.close()

	def open(self):
		self.opened = True

	def close(self):
		self.opened = False

class PatranResultsReader(Reader):
	"""
		Superclass for everything that reads patran-type results files.
	"""
	def __init__(self, noderesults=[], elemresults=[], timestepfile=None):
		"""
			noderesults -- list of nodal results structures
			elemresults -- list of element results structures
			timestepfile -- file giving the times corresponding to each step
		"""
		self.noderesults = noderesults
		self.elemresults = elemresults
		total_results = len(noderesults) + len(elemresults)
		if total_results > 0:
			self.with_results = True
		else:
			self.with_results = False

		if timestepfile:
			self.with_times = True
			self.tf = timestepfile
		else:
			self.with_times = False

	def node_field_iterator(self):
		"""
			For each field (result column) , assemble a big numpy array.  
			We're going to slurp the
			whole thing in anyway, so we may as well make it a (nfields,nnodes,nsteps)
			massive array.
		"""
		class FNodeTypeIt(object):
			def __init__(self, patobj, nfieldobj):
				self.patobj = patobj
				self.nfieldobj = nfieldobj
				self.data = np.zeros((len(nfieldobj), patobj.num_nodes, 
					patobj.num_steps))
				self.n = 0
				self.load_data()

			def load_data(self):
				"""
					Actually read data into the array
				"""
				count = 0; message_at = 50
				for i,j in enumerate(self.patobj.stepnums):
					if j in self.nfieldobj.steps:
						k = self.nfieldobj.steps.index(j)
						fname = self.nfieldobj.fpaths[k]
						ftype = self.nfieldobj.ftypes[k]
						count += 1
						if ftype == "text":
							if( count % message_at == 0 ):
							   print((" %s  files read: %i" % (os.path.basename(fname), count) ))
							self.data[:,:,i] = np.loadtxt(fname).T
						elif ftype == "stream":
							if( count % message_at == 0 ): 
							   print((" %s  files read: %i" % (fname, count) ))
							self.data[:,:,i] = np.fromfile(fname).reshape((
								self.patobj.num_nodes, len(self.nfieldobj))).T
						else:
							raise ValueError("Cannot read results file of type %s." % ftype)

			def __iter__(self):
				return self

			def __next__(self):
				if self.n == len(self.nfieldobj):
					raise StopIteration()

				self.n += 1

				return self.nfieldobj.labels[self.n-1], (1,), self.data[self.n-1,:,:]

		return itertools.chain.from_iterable(
				FNodeTypeIt(self, f) for f in self.noderesults)

	def element_field_iterator(self):
		"""
			For each field, assemble a big numpy array.  We're going to slurp the
			whole thing in anyway, so we may as well make it a (nfields,nelems,nsteps)
			massive array.
		"""
		class FElemTypeIt(object):
			def __init__(self, patobj, nfieldobj):
				self.patobj = patobj
				self.nfieldobj = nfieldobj
				self.data = np.zeros((len(nfieldobj), patobj.num_elems, 
					patobj.num_steps))
				self.n = 0
				self.load_data()

			def load_data(self):
				"""
					Actually read data into the array
				"""
				count = 0; message_at = 50
				for i,j in enumerate(self.patobj.stepnums):
					if j in self.nfieldobj.steps:
						k = self.nfieldobj.steps.index(j)
						fname = self.nfieldobj.fpaths[k]
						ftype = self.nfieldobj.ftypes[k]
						count += 1
						if ftype == "text":
							if( count % message_at == 0 ): 
							   print((" %s  files read: %i" % (os.path.basename(fname), count) ))
							ldata = np.loadtxt(fname).T
							if len(ldata.shape) == 1:
								ldata = ldata.reshape(len(ldata),1)
							self.data[:,:,i] = ldata
						elif ftype == "stream":
							if( count % message_at == 0 ): 
							   print((" %s  files read: %i" % (fname, count) ))
							self.data[:,:,i] = np.fromfile(fname).reshape((
								self.patobj.num_elems, len(self.nfieldobj))).T
						else:
							raise ValueError("Cannot read results file of type %s." % ftype)

			def __iter__(self):
				return self

			def __next__(self):
				if self.n == len(self.nfieldobj):
					raise StopIteration()

				self.n += 1

				return self.nfieldobj.labels[self.n-1], (1,), self.data[self.n-1,:,:]

		return itertools.chain.from_iterable(
				FElemTypeIt(self, f) for f in self.elemresults)

	def integration_field_iterator(self):
		# Not supported
		return []

	def setup_results(self):
#		Parse result file info
		if self.with_results:
			stepnums = set()
			self.num_nvars = 0
			self.num_evars = 0
			for r in self.noderesults:
				stepnums |= set(r.steps)
				self.num_nvars += len(r)
			for r in self.elemresults:
				stepnums |= set(r.steps)
				self.num_evars += len(r)
			self.stepnums = sorted(list(stepnums))
			self.num_steps = len(self.stepnums)
		else:
			self.num_nvars = 0
			self.num_evars = 0
			self.num_steps = 0

#		Read in timestep file, if provided
		if self.with_times:
			self.times = np.loadtxt(self.tf)
			if len(self.times) != self.num_steps:
				raise ValueError("Length of time array from %s does not match model"
						" number of steps!" % self.tf)
		elif self.with_results:
#			Else, default to [0,1,2,...]
			self.times = list(range(self.num_steps))
		else:
			self.times = []

class SimpleReader(PatranResultsReader):
	"""
		Read from new, simple model file and patran-type results files.
	"""
	def __init__(self, filename, *args, **kwargs):
		super(SimpleReader, self).__init__(*args, **kwargs)
		self.filename = filename

		self.open()
		self.read_data()
	
	def read_data(self):
#		Dimensions, nsets, and ivars are fixed
		self.dim = 3
		self.num_nsets = 0
		self.num_ivars = 0

#		Need to read: title, num_nodes, num_elems, num_eblocks
		self.title = 'Generated from simple model'
		psplit = self.filename.split('.')
		if psplit[-1] == 'text':
			self.read_text_data()
		elif psplit[-1] == 'str':
			self.read_stream_data()
		else:
			raise ValueError("Unknown simple model file with suffix %s." % psplit[-1])
		
#		Setup the results reader
		self.setup_results()

	def node_iterator(self):
		"""
			Return the coordinates of each node, in order of ID.
		"""
		return self.nodes

	def elem_iterator(self):
		"""
			Return connectivity, etype tuples for each element, in ID order.
		"""
		return zip(self.conn, self.etype)

	def eblk_iterator(self):
		"""
			Return tuples consisting of the name of the block and a list of
			associated elements.
		"""
		sgroups = np.array(self.group)
		inorder = np.argsort(sgroups)
		sgroups = sgroups[inorder]

		groups = []
		names = []
		
		g = -1
		ngroup = []
		for i,e in enumerate(inorder):
			if sgroups[i] != g:
				if len(ngroup) > 0:
					groups.append(ngroup)
				g = sgroups[i]
				ngroup = [e]
				names.append(str(g))
			else:
				ngroup.append(e)

		groups.append(ngroup)

		return zip(names, groups)

	def nset_iterator(self):
		# Simple model doesn't support node sets.  However, hardcode in a node set called
		# "All Nodes" with all the nodes.
		return (("All Nodes", list(range(self.num_nodes))),)

	def read_text_data(self):
		"""
			Read data in from the text-type files.
		"""
		self.elements = []
		with open(self.filename, 'r') as f:
			reading = 1
			for line in f:
				line = line.strip()
				if line[0] == '#':
					continue
				sline = line.split()
				if len(sline) == 2 and reading == 1:
					self.num_nodes = int(sline[0])
					self.num_elems = int(sline[1])
					self.nodes = np.zeros((self.num_nodes,3))
					reading = 2
					ncount = 0
				elif len(sline) == 3 and reading == 2:
					self.nodes[ncount] = [ float(i) for i in sline ]
					ncount += 1
					if ncount == self.num_nodes:
						reading = 3
						ecount = 0
						self.conn = []
						self.group = []
						self.etype = []
				elif reading == 3:
					self.etype.append(patran_types[int(sline[0])])
					self.group.append(int(sline[1]))
#					Figure out where the connectivity line actually ends
					intline = [ int(i) for i in sline[:] ]
					for i, e in reversed(list(enumerate(intline))):
						if e != 0:
							n = i+1
							break
					self.conn.append([ int(i)-1 for i in sline[2:n]])
#					Special rules for some weird elements types
					if self.etype[-1] == "WEDGE" and n == 17:
						# Actually a trint12, collapse nodes
						self.conn[-1][9] = self.conn[-1][0]
						self.conn[-1][10] = self.conn[-1][1]
						self.conn[-1][11] = self.conn[-1][2]
				else:
					raise ValueError("Found unknown state while reading file")

class PatranReader(PatranResultsReader):
	"""
		Reads a patran neutral file + patran-type results files

		IMPORTANT NOTE: This assumes that the node and element packets of the
		neutral file are in order.  We do raise an exception if this isn't
		the case.
	"""
	def __init__(self, neutralfile, *args, **kwargs):
		"""
			neutralfile -- filename of patran neutral file
		"""
		super(PatranReader, self).__init__(*args, **kwargs)
		self.neut_file = neutralfile

		self.open()
		self.read_opening()

	def read_opening(self):
		"""
			Read the number of nodes, elements, and timesteps.
		"""
		print("Reading basic data...")
#		Dimensions, nsets, and ivars are fixed
		self.dim = 3
		self.num_nsets = 0
		self.num_ivars = 0

#		Read remaining header data
		total = 0
		for packet,data in PatranNeutralIt(self.neut_file, writeout=wo):
			if packet[0] == 25:
				self.title = ''.join(data[0]).strip()
				total += 1
			elif packet[0] == 26:
				self.num_nodes = packet[4]
				self.num_elems = packet[5]
				total += 1
			if total == 2:
				break

		if total != 2:
			raise ValueError("Could not read header packets 25 and 26!")

#		Unfortunately, we don't write consistently the number of element configs
#		Read through the file, looking at element packets, data card 1, slot 2
#		to total up the configurations
		configs = set()
		total = 0
		for packet,data in PatranNeutralIt(self.neut_file, writeout=wo):
			if packet[0] == 2:
				configs.add(data[0][1])
				total += 1
			if total == self.num_elems:
				break
		self.num_eblocks = len(configs)

#		Setup for the patran results
		self.setup_results()

	def node_iterator(self):
		"""
			Return the coordinates of each node, in order of ID.

			The iterator assumes that the neutral file has sorted the node packets
			into order by ID.  If not, it will raise an exception.
		"""
		class PatranNodeIt(PatranNeutralIt):
			def __init__(self, pfile, nnodes, *args, **kwargs):
				super(PatranNodeIt,self).__init__(pfile, *args, **kwargs)
				self.nnodes = nnodes
				self.nn = 0

			def __next__(self):
				if self.nn == self.nnodes:
					raise StopIteration()

				packet, data = next(super(PatranNodeIt, self))
				while packet[0] != 1:
					packet, data = next(super(PatranNodeIt, self))
				
				self.nn += 1
				if packet[1] != self.nn:
					raise ValueError("Node %i in neutral file is not sorted by ID!"
							% packet[1])

				return np.array(data[0][:3])

		return PatranNodeIt(self.neut_file, self.num_nodes, writeout=wo)
	
	def elem_iterator(self):
		"""
			Return connectivity, etype tuples for each element, in ID order.

			If the elements in the file aren't in ID order, this will raise 
			an exception.
		"""
		class PatranElemIt(PatranNeutralIt):
			def __init__(self, pfile, nelem, *args, **kwargs):
				super(PatranElemIt,self).__init__(pfile, *args, **kwargs)
				self.nelem = nelem
				self.ne = 0

			def __next__(self):
				if self.ne == self.nelem:
					raise StopIteration()

				packet, data = next(super(PatranElemIt, self))
				while packet[0] != 2:
					packet, data = next(super(PatranElemIt, self))
				
				self.ne += 1
				if packet[1] != self.ne:
					raise ValueError("Element %i in neutral file is not sorted by ID!"
							% packet[1])
				
				cnode = data[0][0]
				joined = np.hstack(data[1:])
				conn = joined[:cnode]

				raw_etype = packet[2]
				etype = patran_types[raw_etype]

				if etype == 'UNKNOWN':
					raise ValueError("Element %i has unknown type %i." % (packet[1],
						raw_etype))

				nconn = np.array(conn, dtype=int)-1

#				Special rules for some weird elements types
				if etype == "WEDGE" and cnode == 15:
					# Actually a trint12, collapse nodes
					nconn[9] = nconn[0]
					nconn[10] = nconn[1]
					nconn[11] = nconn[2]

				return (nconn,etype)

		return PatranElemIt(self.neut_file, self.num_elems, writeout=wo)	

	def nset_iterator(self):
		# Patran doesn't support node sets.  However, hardcode in a node set called
		# "All Nodes" with all the nodes.
		return (("All Nodes", list(range(self.num_nodes))),)

	def eblk_iterator(self):
		"""
			Assume elements are blocked by the **config** ID (not property ID)

			Unfortunately, all we can do is run through and add elements to
			the appropriate list.
		"""
		blocks = collections.defaultdict(list)
		nelem = 0
		for packet,data in PatranNeutralIt(self.neut_file, writeout=wo):
			if packet[0] == 2:
				nelem += 1
				blocks[str(data[0][1])].append(packet[1]-1)

			if nelem == self.num_elems:
				break

		return iter(blocks.items())

class ExodusIIReader(Reader):
	"""
		Reads from an ExodusII file.
	"""
	def __init__(self, filename):
		super(ExodusIIReader, self).__init__()
		self.filename = filename
		self.open()

	def open(self):
		self.ncdf = sio.netcdf_file(self.filename, mode='r')
		self.opened = True

	def close(self):
		self.ncdf.close()
		self.opened = False

	def node_iterator(self):
		"""
			Return the coordinates of each node, in the file order.
		"""
		class ExodusNodeIt(object):
			def __init__(self, exo_obj):
				self.exo = exo_obj
				self.n = 0

			def __iter__(self):
				return self

			def __next__(self):
				if self.n == self.exo.num_nodes:
					raise StopIteration
				else:
					ccoord = np.zeros((self.exo.dim,))
					cnames = ['coordx', 'coordy', 'coordz']
					for i,name in zip(list(range(self.exo.dim)), cnames):
						ccoord[i] = self.exo.ncdf.variables[name][self.n]
					self.n += 1
					return ccoord

		return ExodusNodeIt(self)

	def elem_iterator(self):
		"""
			Return an iterator which will pass through each element in order 
			(regardless of block) and return the connectivity 
			(in terms of zero-start node numbers) and the element type (for now, the raw exodus type)

			Unfortunately, the Exodus files don't store the elements sequentially
			anywhere.  Instead, it essentially stores each element block in order,
			and then the true element number must be dereferenced from a map.

			This means the iterator is a bit more complicated.  It's also cheaper
			to sort the map once and get a reverse map.
		"""
		class ExodusElemIt(object):
			def __init__(self, exo_obj):
				self.exo = exo_obj
				self.inv_map = list(range(self.exo.num_elems))
				self.inv_map.sort(key = lambda i: 
						self.exo.ncdf.variables['elem_num_map'][i])
				self.n = 0

			def __iter__(self):
				return self

			def __next__(self):
				if self.n == self.exo.num_elems:
					raise StopIteration
				else:
					stored_enum = self.inv_map[self.n]
					blk, offset = self.exo.enum_to_blk(stored_enum)
					raw_connect = self.exo.ncdf.variables['connect'+str(blk)]
					connect = raw_connect[offset]-1
					conn_size = self.exo.ncdf.dimensions['num_nod_per_el'+str(blk)]
					self.n += 1
					return connect[0:conn_size], raw_connect.elem_type

		return ExodusElemIt(self)

	def nset_iterator(self):
		"""
			Return an iterator which returns a tuple of (name, list of nodes) for
			each node set
		"""
		class ExodusNSetIt(object):
			def __init__(self, exo_obj):
				self.exo = exo_obj
				self.n = 0

			def __iter__(self):
				return self

			def __next__(self):
				if self.n == self.exo.num_nsets:
					raise StopIteration
				else:
					name = ''.join(self.exo.ncdf.variables['ns_names'][0])
					nodes = self.exo.ncdf.variables['node_ns'+str(self.n+1)][:]-1
					self.n += 1
					return name, nodes

		return ExodusNSetIt(self)

	def eblk_iterator(self):
		"""
			Return an iterator which returns a tuple of (name, list of elements) for
			each element block
		"""
		class ExodusEBlkIt(object):
			def __init__(self, exo_obj):
				self.exo = exo_obj
				self.n = 0

			def __iter__(self):
				return self

			def __next__(self):
				if self.n == self.exo.num_eblocks:
					raise StopIteration
				else:
					name = ''.join(self.exo.ncdf.variables['eb_names'][0])
					blkszs = self.exo.eblk_sizes
					blkszs.insert(0,0)
					blkszs = np.cumsum(blkszs)
					elems = self.exo.ncdf.variables['elem_num_map'][blkszs[self.n]:blkszs[self.n+1]]-1
					self.n += 1
					return name, elems

		return ExodusEBlkIt(self)

	def node_field_iterator(self):
		"""
			This one is annoying.  First we need to iterate over all fields, then
			over all nodes/elements (in stored order), then over time steps.  Or
			some other combination of the three.
		"""
		class ExodusNodeFieldsIt(object):
			"""
				Returns field name, field shape, and field iterator
			"""
			def __init__(self, exo_obj):
				self.exo = exo_obj
				self.n = 0

			def __iter__(self):
				return self

			def __next__(self):
				if self.n == self.exo.num_nvars:
					raise StopIteration
				else:
					cfieldname = ''.join(self.exo.ncdf.variables['name_nod_var'][self.n])
					self.n += 1
					return (cfieldname, (1,), ExodusNodeFieldIt(self.exo, self.n-1))

		class ExodusNodeFieldIt(object):
			def __init__(self, exo_obj, field_num):
				self.exo = exo_obj
				self.n = 0
				self.fn = field_num

			def __iter__(self):
				return self

			def __next__(self):
				if self.n == self.exo.num_nodes:
					raise StopIteration()
				else:
					ff_name = 'vals_nod_var'+str(self.fn+1)
					raw_data = self.exo.ncdf.variables[ff_name][:, self.n]
					self.n += 1
					return np.reshape(raw_data, (self.exo.num_steps,1))
		
		return ExodusNodeFieldsIt(self)

	def element_field_iterator(self):
		"""
			This one is annoying.  First we need to iterate over all fields, then
			over all nodes/elements (in stored order), then over time steps.  Or
			some other combination of the three.

			Additionally, we have the problem of the data possibly being "masked" --
			not present for a particular element block.
		"""
		class ExodusElementFieldsIt(object):
			"""
				Returns field name, field shape, and field iterator
			"""
			def __init__(self, exo_obj):
				self.exo = exo_obj
				self.n = 0

			def __iter__(self):
				return self

			def __next__(self):
				if self.n == self.exo.num_evars:
					raise StopIteration
				else:
					cfieldname = ''.join(self.exo.ncdf.variables['name_elem_var'][self.n])
					self.n += 1
					return (cfieldname, (1,), ExodusElementFieldIt(self.exo, self.n-1))

		class ExodusElementFieldIt(object):
			def __init__(self, exo_obj, field_num):
				self.exo = exo_obj
				self.inv_map = list(range(self.exo.num_elems))
				self.inv_map.sort(key = lambda i: 
						self.exo.ncdf.variables['elem_num_map'][i])
				self.n = 0
				self.fn = field_num

			def __iter__(self):
				return self

			def __next__(self):
				if self.n == self.exo.num_elems:
					raise StopIteration()
				else:
					stored_enum = self.inv_map[self.n]
					blk, offset = self.exo.enum_to_blk(stored_enum)
					if self.exo.ncdf.variables["elem_var_tab"][blk-1][self.fn] != 1:
						return np.zeros((self.exo.num_steps,1))
					else:
						ff_name = 'vals_elem_var'+str(self.fn+1)+'eb'+str(blk)
						raw_data = self.exo.ncdf.variables[ff_name][:, offset]
						self.n += 1
						return np.reshape(raw_data, (self.exo.num_steps,1))
		
		return ExodusElementFieldsIt(self)
	
	def integration_field_iterator(self):
		return []

	@property
	def title(self):
		return self.ncdf.title.strip()

	@property
	def dim(self):
		return self.ncdf.dimensions['num_dim']

	@property
	def num_steps(self):
		return len(self.ncdf.variables['time_whole'][:])

	@property
	def num_nodes(self):
		return self.ncdf.dimensions['num_nodes']

	@property
	def num_elems(self):
		return self.ncdf.dimensions['num_elem']

	@property
	def num_eblocks(self):
		return self.ncdf.dimensions['num_el_blk']

	@property
	def num_nsets(self):
		if 'num_node_sets' not in self.ncdf.dimensions:
			return 0
		else:
			return self.ncdf.dimensions['num_node_sets']

	@property
	def num_nvars(self):
		if 'num_nod_var' not in self.ncdf.dimensions:
			return 0
		else:
			return self.ncdf.dimensions['num_nod_var']
	
	@property
	def num_evars(self):
		if 'num_elem_var' not in self.ncdf.dimensions:
			return 0
		else:
			return self.ncdf.dimensions['num_elem_var']
	
	@property
	def num_ivars(self):
		"""
			Exodus cannot store integration point variables.
		"""
		return 0

	@property
	def times(self):
		"""
			Return the actual time step times.
		"""
		return self.ncdf.variables['time_whole'][:]

	def enum_to_blk(self, el):
		"""
			Take a zero-indexed *stored* element number and return its block number
			(one-indexed) and offset into the block (zero-indexed)
		"""
		blkszs = self.eblk_sizes
		blkszs.insert(0,0)
		blkszs = np.cumsum(blkszs)
		for i in range(len(blkszs)-1):
			s = blkszs[i]
			e = blkszs[i+1]
			if s <= el < e:
				return (i+1,el-blkszs[i])

		raise ValueError("Element #%i does not seem to be in a block." % e)

	@property
	def eblk_sizes(self):
		szs = []
		for i in range(self.num_eblocks):
			szs.append(self.ncdf.dimensions['num_el_in_blk'+str(i+1)])

		return szs


