from array import array
import math as magic
import numpy
import warnings

# DEPRECATED as it has been replaced by a 3d array in class Model
class Macrostate(object):
	"""
	Deprecated as it has been replaced by a 3d array in class Model
	"""
	thisState = -1;								# enum value of this macrostate
	ensembleSize = 1;							# int
	nPositions = 1;								# int, number of positions on sequence examined
	positionOffset = 0;							# offset for indexing
	backrubTemp = 1.0;							# float
	boltzmannTemp = 1.0;						# float
	steepness = 1.0;							# float

	# TODO: convert this to a true 2D array with numpy
	residueEnergies;							# double[position][residue]

	#
	def __init__(self, state, ensembleSize, nPositions, backrubTemp, boltzmannTemp, positionOffset):
		warnings.warn("This class has been superceded by a 3D numeric array", DeprecationWarning);
		self.thisState = state;
		self.ensembleSize = ensembleSize;
		self.nPositions = nPositions;
		self.backrubTemp = backrubTemp;
		self.boltzmannTemp = boltzmannTemp;
		self.positionOffset = positionOffset;

		# allocate arrays
		self.residueEnergies = numpy.ones([self.nPositions, 20]);

	def setPositionEnergies(self, position, energies):
		"""
		Sets the energies at a certain position

		@param position		position of this set of energies
		@param energies		length-20 array of energy values
		@throws Exception	when an energies list already exists for this location
		"""
		if self.residueEnergies[position - self.positionOffset][0] != 1:
			raise Exception("Something something already copied here");

		for i in range(20):
			self.residueEnergies[position - self.positionOffset][i] = energies[i];

	# TODO: where doe this belong? prolly in Model. 
	def evalFitness(self, offset):
		"""
		Evaluates the macrostate fitness

		@param offset		offset to be used
		@return float[] corresponding to fitness for each of the residues.
		"""
		# TODO: implement this
		fitness = numpy.array([20]);