# Various patran helper objects
import collections
import re

# Valid warp3d result file rules:
initial = ['w']
location = ['n', 'e']
types = ['r','e','v','t','a','d','s','m']
suffixes = ['_text', '_stream']

mats = []

# Big lookup table: for each description (i.e. e = strains) given a huge
# list of the field names
results_lookup = {
		'r': ['Reactionx', 'Reactiony', 'Reactionz'],
		'e': ['Strain_xx', 'Strain_yy', 'Strain_zz', 'Strain_xy', 'Strain_yz',
			'Strain_xz', 'Effective_mises_strain', 'Strain_invariant_1',
			'Strain_invariant_2', 'Strain_invariant_3', 'Strain_principle_min',
			'Strain_principle_inter', 'Strain_principle_max',
			'Strain_dir_cosine_l1', 'Strain_dir_cosine_m1', 'Strain_dir_cosine_n1',
			'Strain_dir_cosine_l2', 'Strain_dir_cosine_m2', 'Strain_dir_cosine_n2',
			'Strain_dir_cosine_l3', 'Strain_dir_cosine_m3', 'Strain_dir_cosine_n3'],
		'v': ['Velocityx', 'Velocityy', 'Velocityz'],
		't': ['Temperature', 'Reference temperature'],
		'a': ['Accelerationx', 'Accelerationy', 'Accelerationz'],
		'd': ['Displacementx', 'Displacementy', 'Displacementz'],
		's': ['Stress_xx', 'Stress_yy', 'Stress_zz', 'Stress_xy', 'Stress_yz',
			'Stress_xz', 'Stress_work_density', 'mises_stress', 'C1_Mat_val',
			'C2_Mat_val', 'C3_Mat_val', 'Stress_invariant_1', 'Stress_invariant_2',
			'Stress_invariant_3', 'Stress_principle_min', 'Stress_principle_inter',
			'Stress_principle_max',
			'Stress_dir_cosine_l1', 'Stress_dir_cosine_m1', 'Stress_dir_cosine_n1',
			'Stress_dir_cosine_l2', 'Stress_dir_cosine_m2', 'Stress_dir_cosine_n2',
			'Stress_dir_cosine_l3', 'Stress_dir_cosine_m3', 'Stress_dir_cosine_n3']
		}

# Don't change, just compiles the regular expression
reor = lambda s: '('+'|'.join(s)+')'
restr = ''.join([reor(initial), reor(location), reor(types), r'(\d*)',
	reor(suffixes),reor(mats)+'?',r'(?:.gz)?','$'])
reg = re.compile(restr)

def recompile():
	"""
		Blargh: after changing the list of valid mats you need to recompile the
		stupid regex...
	"""
	reor = lambda s: '('+'|'.join(s)+')'
	restr = ''.join([reor(initial), reor(location), reor(types), r'(\d*)',
		reor(suffixes),reor(mats)+'?',r'(?:.gz)?','$'])
	reg = re.compile(restr)
	return reg

def valid_result_file(fname):
	"""
		Determine if a file fname is a valid Warp3D results file, according to 
		the above rules.
	"""
	if reg.search(fname):
		return True
	else:
		return False

def parse_result_fname(fname):
	"""
		Separate a result file name into the location, type, step #, and file type
	"""
	mtchs = reg.match(fname)
	if mtchs.group(6):
		typ = mtchs.group(3)+mtchs.group(6)
	else:
		typ = mtchs.group(3)

	if not mtchs:
		raise ValueError("Filename %s is not a valid results file!" % fname)

	return (mtchs.group(2), typ, int(mtchs.group(4)), 
			mtchs.group(5).lstrip('_'))

# Element type definitions
pt = {
		2: 'TRUSS',
		3: 'TRIANGLE',
		4: 'QUAD',
		5: 'TETRA',
		7: 'WEDGE',
		8: 'HEX'}

patran_types = collections.defaultdict(lambda: 'UNKNOWN', pt)

# Fortran packet definitions.  Dictionary linking the packet number to 
# the relevant data formats.  If KC > # of formats just repeat the last one.
packets = { 	25: ('(20A4)',),
							26: ('(3A4,2A4,3A4)',),
							1:  ('(3E16.9)','(I1,1A1,I8,I8,I8,2X,6I1)'),
							2:  ('(I8,I8,I8,I8,3E16.9)', '(10I8)'),
							3:  ('(5E16.9)',),
							4:	('raw',),
							5:	('(5E16.9)',),
							6:  ('(I1,I1,I1,6I1,8I1,I2)', '(5E16.9)'),
							7:  ('(I8,6I1)', '(5E16.9)'),
							8:  ('(I8,6I1)', '(5E16.9)'),
							10: ('(E16.9',),
							11: ('(E16.9',),
							14: ('(3A12)','(2I8,E16.9)','2(2I8,E16.9)'),
							15: ('(E16.9)',),
							16: ('(I1,1X,8I1)','(5E16.9)'),
							17: ('(I1,1X,8I1)','(5E16.9)'),
							18: ('(I1,1X,8I1)','(5E16.9)'),
							19: ('(6I8,2X,8I1)',),
							21: ('(A12)','(10I8)'),
							31: ('(3E16.9)',),
							32: 3*('(5E16.9,5E16.9,2E16.9,2I8)',),
							33: 10*('9(5E16.9/5), 3E16.9/2E16.9,4I8',),
							42: ('(I8,3I8,I8,5I8)','(I8,I8,7X,1A1,7I8)','(3E16.9)','(3E16.9)','(8I8, I8, I8)','(3E16.9, I8, I8)'),
							43: ('(I8,3I8,I8,5I8)','(I8,I8,7X,1A1,7I8)','(3E16.9)','(3E16.9)','(8I8, I8, I8)','(3E16.9, I8, I8)'),
							44: ('(I8,3I8,I8,5I8)','(I8,I8,7X,1A1,7I8)','(3E16.9)','(3E16.9)','(8I8, I8, I8)','(3E16.9, I8, I8)'),
							45: ('(6I8)','(10I8)'),
							99: tuple()
							}

class PatranResultsDesc(object):
	"""
		Descriptor class for a set of results files at a common location and type.
		Contains step numbers, full paths of datafiles, and file types.
	"""
	def __init__(self, loc, rtype):
		self.loc = loc
		self.rtype = rtype
		self.steps = []
		self.fpaths = []
		self.ftypes = []
		self.labels = []
	
	def __len__(self):
		return len(self.labels)

	@property
	def description(self):
		"""
			Return the 2 character string describing the type.
		"""
		return self.loc+self.rtype

	def add_file_at_step(self, step, ftype, path):
		"""
			Add datafile with path at step #.
		"""
		if step in self.steps:
			raise ValueError("Results file for type %s at step %i already exists."
					% (self.description,step))

		# Find the correct place to insert, this gives the index beyond the spot
		# to insert
		i = next((i for i,v in enumerate(self.steps) if v > step), len(self.steps))
		self.steps.insert(i,step)
		self.ftypes.insert(i,ftype)
		self.fpaths.insert(i,path)

