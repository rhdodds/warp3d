import numpy as np
import scipy as sp
import scipy.io as sio

# Basic data about ExodusII files
char_dtype = '|S1'

api_version = 5.14
fp_w_sz = 8 # Does this need to change?
m_name_l = 32
version = api_version

four = 4
len_line = 81
len_name = 33
len_string = 33

coor_names = ['X','Y','Z']
coord_labels = ['coordx', 'coordy', 'coordz']

def transform_string(raw, l):
  """
    Transform a string to an np array of dtype '|S1' of length l.
    Also append a null-terminal to it just in case
  """
  if len(raw) > l-1:
    raw = raw[:l-1]
  
  raw += '\0'
  
  ns = np.empty((l,), dtype=char_dtype)
  ns[:] = ''
  for i in range(len(raw)):
    ns[i] = raw[i]
  
  return ns

class Writer(object):
  def __init__(self, filename):
    self.filename = filename
    self.opened = False

class ExodusIIWriter(Writer):
  def __init__(self, filename, large_file = True):
    super(ExodusIIWriter, self).__init__(filename)
    self.large_file = large_file
    self.open()

  def write(self, model):
    """
      Write everything to file, in stages.
    """
    self.setup(model)
    self.write_basic(model)
    self.write_size_dimensions(model)
    self.write_nodes(model)
    print " ...coordinates written"
    self.write_elements_by_block(model)
    print " ...element data written"
    self.write_node_blocks(model)
    self.write_timestep_basic(model)
    self.write_nodal_variables(model)
    print " ...nodal results written"
    self.write_element_variables(model)
    print " ...element results written"
  
  def setup(self, model):
    """
      Do some optimization on accessing nodes by name.
    """
#     Hash the nodes to map indexing into their new names easy
    self.nmap = { node: i for (i, node) in enumerate(model.nodes) }
    self.emap = { elem: i for (i, elem) in enumerate(model.elements) }

  def write_basic(self, model):
    """
      Write basic data to the exodus file.  Only thing actually coming
      from model is the title.
    """
    self.ncdf.title = model.title.ljust(len_line-1)+"\0"
    self.ncdf.api_version = api_version
    if self.large_file:
      self.ncdf.file_size = 1
    else:
      self.ncdf.file_size = 0
    self.ncdf.floating_point_word_size = fp_w_sz
    self.ncdf.maximum_name_length = m_name_l
    self.ncdf.version = version
    self.ncdf.sync()
  
  def write_size_dimensions(self, model):
    """
      Write out all the string/line/name size dimensions
    """
    self.ncdf.createDimension('four', four)
    self.ncdf.createDimension('len_line', len_line)
    self.ncdf.createDimension('len_name', len_name)
    self.ncdf.createDimension('len_string', len_string)

    self.ncdf.sync()

  def write_nodes(self, model):
    """
      Write out all the basic nodal data (coordinates).
    """
    # Dimensions and dimension names
    self.ncdf.createDimension('num_dim', model.dim)
    if model.dim > 3:
      raise ValueError("No default names for dimensions > 3!")
    cn = self.ncdf.createVariable('coor_names', char_dtype,
        ('num_dim','len_name'))
    names = np.vstack((transform_string(coor_names[i], len_name) for i 
      in range(model.dim)))
    cn[:] = names

    # Nodes
    self.ncdf.createDimension('num_nodes', model.num_nodes)
    
    if self.large_file:
      for i in range(model.dim):
        coords = np.array([node.coords[i] for node in model.nodes])
        cfile = self.ncdf.createVariable(coord_labels[i], 'd', ('num_nodes',))
        cfile[:] = coords
    else:
      coords = np.array([node.coords for node in model.nodes])
      cfile = self.ncdf.createVariable('coord', 'd', ('num_dim', 'num_nodes'))
      cfile[:] = coords.T
    self.ncdf.sync()

  def write_elements_by_block(self, model):
    """
      Exodus files use blocks of elements.  Now we have to assume that
      each element is present in one and only one block.  May need to
      enforce this in the model.
    """
    self.ncdf.createDimension('num_el_blk', model.num_eblks)
    self.ncdf.createDimension('num_elem', model.num_elems)
    
    # Names of all blocks
    enames = self.ncdf.createVariable('eb_names', char_dtype, 
        ('num_el_blk', 'len_name'))
    names = np.vstack((transform_string(name, len_name) for name
      in model.eblks.keys() ))
    enames[:] = names

    # Status of all blocks (may be wrong...)
    estat = self.ncdf.createVariable('eb_status', 'i',
        ('num_el_blk',))
    estat[:] = np.array([1], dtype='i')

#   Currently no real properties -- just block IDs
    eprop = self.ncdf.createVariable('eb_prop1', 'i',
        ('num_el_blk',))
    eprop[:] = np.array(range(1,model.num_eblks+1), dtype='i')
    eprop.name = "ID"+"\0"

    e_map = np.zeros((model.num_elems,), dtype='i')
    j = 0
    for i, (name, elems) in enumerate(model.eblks.iteritems()):
#     Basics
      self.ncdf.createDimension('num_el_in_blk'+str(i+1), 
          len(elems))
      nconn = len(elems[0].conn)
      self.ncdf.createDimension('num_nod_per_el'+str(i+1), nconn)

#     Connectivity
      blk_conn = np.zeros((len(elems), nconn), dtype='i')

      for k,e in enumerate(elems):
        e_map[j] = self.emap[e]
#       e_map[j] = model.elements.index(e)
        j += 1
#       blk_conn[k,:] = [ model.nodes.index(n) for n in e.conn ]
        blk_conn[k,:] = [ self.nmap[n] for n in e.conn ]

      file_conn = self.ncdf.createVariable('connect'+str(i+1),
          'i', ('num_el_in_blk'+str(i+1), 'num_nod_per_el'+str(i+1)))
      file_conn[:] = blk_conn+1
#     Assign element type.  Eventually translate back and forth b/w
#     my element type and the Exodus types.  For testing, just assume
#     the Exodus type.
      file_conn.elem_type = elems[0].etype+"\0"

#   Map
    fmap = self.ncdf.createVariable('elem_num_map', 'i', ('num_elem',))
    fmap[:] = e_map + 1
    
    self.ncdf.sync()

  def write_node_blocks(self, model):
    """
      Write out node block data
    """
    self.ncdf.createDimension('num_node_sets', len(model.nsets))

    if len(model.nsets) == 0:
      return

#   Names
    nnames = self.ncdf.createVariable('ns_names', char_dtype, 
        ('num_node_sets', 'len_name'))
    names = np.vstack((transform_string(name, len_name) for name
      in model.nsets.keys() ))
    nnames[:] = names

#   Pointless status
    ns = self.ncdf.createVariable('ns_status', 'i', ('num_node_sets',))
    ns[:] = np.array([1], dtype='i')

#   Pointless properties
    nop = self.ncdf.createVariable('ns_prop1', 'i',
        ('num_node_sets',))
    nop[:] = np.array([1], dtype='i')
    nop.name = "ID"+"\0"

    for i,(name,nodes) in enumerate(model.nsets.iteritems()):
#     Actual nodes
      self.ncdf.createDimension('num_nod_ns'+str(i+1), len(nodes))
#     nn = np.array([model.nodes.index(node) for node in nodes])+1
      nn = np.array([self.nmap[node] for node in nodes]) + 1
      nfile = self.ncdf.createVariable('node_ns'+str(i+1), 'i',
          ('num_nod_ns'+str(i+1),))
      nfile[:] = nn

    self.ncdf.sync()

  def write_timestep_basic(self, model):
    """
      Write basic timestep data.
    """
    self.ncdf.createDimension('time_step', model.timesteps)
    tw = self.ncdf.createVariable('time_whole', 'd', ('time_step',))
    if model.timesteps != 0:
      tw[:] = model.times
    self.ncdf.sync()

  def write_nodal_variables(self, model):
    """
      Write out all the nodal variables.
    """
    if len(model.nfields) == 0:
      return

    self.ncdf.createDimension('num_nod_var', len(model.nfields))

    # Names of all the variables
    nvnames = self.ncdf.createVariable('name_nod_var', char_dtype, 
        ('num_nod_var', 'len_name'))
    names = np.vstack((transform_string(name, len_name) for name
      in model.nfields.keys() ))
    nvnames[:] = names

    # Dump each variable
    if self.large_file:
      for i,field in enumerate(model.nfields.keys()):
        nv = self.ncdf.createVariable('vals_nod_var'+str(i+1), 'd',
            ('time_step','num_nodes'))
        var = np.zeros((model.timesteps,model.num_nodes))
        for j,node in enumerate(model.nodes):
          var[:,j] = np.reshape(node.fields[field], (model.timesteps,))
        nv[:] = var
    else:
      fields = np.zeros((model.timesteps, len(model.nfields), model.num_nodes))
      for i, field in enumerate(model.nfields.keys()):
        fi = np.array([node.fields[field] for node in model.nodes])
        fields[:,i,:] = fi.T
      nv = self.ncdf.createVariable('vals_nod_var', 'd', ('time_step',
        'num_nod_var', 'num_nodes'))
      nv[:] = fields

    self.ncdf.sync()

  def write_element_variables(self, model):
    """
      Write out all the element variables.  This is complicated by
      the fact that we need to do it by block.

      Note the Exodus format lets you "mask out" element blocks if you
      don't have results data for that block.  I'm ignoring that and have
      already written zeros if necessary in the intermediate representation.
    """
    if len(model.efields) == 0:
      return
    self.ncdf.createDimension('num_elem_var', len(model.efields))

    # Names of all the variables
    evnames = self.ncdf.createVariable('name_elem_var', char_dtype, 
        ('num_elem_var', 'len_name'))
    names = np.vstack((transform_string(name, len_name) for name
      in model.efields.keys() ))
    evnames[:] = names

    # Write a fake truth table (all blocks have all variables)
    etab = self.ncdf.createVariable('elem_var_tab', 'i',
        ('num_el_blk', 'num_elem_var'))
    etab[:] = np.ones((len(model.eblks),len(model.efields)), dtype='i')

    # For each block...
    for b, block in enumerate(model.eblks.values()):
      for v, fname in enumerate(model.efields.keys()):
        ev = self.ncdf.createVariable(
            'vals_elem_var'+str(v+1)+'eb'+str(b+1),
            'd', ('time_step','num_el_in_blk'+str(b+1)))
        mev = np.zeros((model.timesteps,len(block)))
        for k,e in enumerate(block):
          mev[:,k] = np.reshape(e.fields[fname], (model.timesteps,))
        ev[:] = mev

    self.ncdf.sync()

  def open(self):
    if self.large_file:
      self.ncdf = sio.netcdf_file(self.filename, mode='w',version=2)
    else:
      self.ncdf = sio.netcdf_file(self.filename, mode='w', version=2) # Sigh
    self.opened = True

  def close(self):
    self.ncdf.flush()
    self.ncdf.close()
    self.opened = False

  def __del__(self):
    if self.opened:
      self.close()
  
