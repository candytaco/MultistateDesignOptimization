import math as magic
import ast
import numpy
from io import *
from enumeration import enum

# should the macrostates be hard-coded? probably not if this ends up being actually used for tuning other models...
MACROSTATES_T = enum("E-DHF-NADPH", "E-NADPH", "E-OPEN", "E-THF", "E-THF-NADPX", "TS");
RESIDUES = enum('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');

class Model:
	"""
	A multistate design model
	"""
	
	# I need these statements for intellisense.

	# vars used by all instances
	MACROSTATES = MACROSTATES_T;					# enum of the macrostates TODO: replace evertually with assignment in constructor
	nMacrostates = 6;								# int
	ensembleSize = 0;								# int
	nPositions = 0;									# int, number of positions on sequence examined
	positionOffset = 0;								# offset for positions
	backrubTemp = 0.0;								# float
	boltzmannTemp = 0.0;							# float
	weights = numpy.array(0);						# float[] relative weights of the macrostates
	steepness = 0.0;								# float steepness of siggmoid (s)
	fitnesses = numpy.array(0);						# double[position][residue fitness] calculated fitnesses of residues at each location
	frequencies = numpy.array(0);					# double[position][residue frequency] calculated frequencies of residues at each location
	macrostateResidueEnergies = numpy.array(0);		# double[position][residue energy][macrostate]
	recovery = -1.0;								# float assigned by the outside similiartyMeasure to how well this mode recovers the sequence
	# TODO: check if macrostatesUsed is actually useful
	macrostatesUsed = numpy.array(False);			# bool[], macrostates examined during optimization
	
	# vars used when microstate data is involved
	useAltAveragingMethod = False;					# use accepted avg method or (the alt) xinjie's expression?
	isFrequenciesCalculated = False;				# prevent unnecessary calculations
	useMicrostateData = False;						# do we have data from individual microstates?
	microstateResidueEnergies = numpy.array(0);		# double[position][residue energy][macrostate][microstate]
	microstateCounts = numpy.array(0);				# double[position][macrostate] number of microstates
	microstatesUsed = numpy.array(0);				# int[macrostate][microstate index], the microstates used to calculate the macrostates

	def __init__(self, ensembleSize:int, backrubTemp:float, boltzmannTemp:float, weights:numpy.array, steepness:float, positions:int, positionOffset:int, useMicrostateData:bool = False, useAltAverageMethod:bool = False):
		"""
		Default constructor

		"""
		self.ensembleSize = ensembleSize;
		self.backrubTemp = backrubTemp;
		self.boltzmannTemp = boltzmannTemp;
		self.weights = weights;
		self.steepness = steepness;
		self.nPositions = positions;
		self.positionOffset = positionOffset;
		self.isFrequenciesCalculated = False;
		self.useMicrostateData = useMicrostateData;
		self.useAltAveragingMethod = useAltAverageMethod;
		self.macrostatesUsed = numpy.array([True] * self.nMacrostates);

		# allocate arrays
		self.macrostateResidueEnergies = numpy.zeros([self.nPositions, 20, self.nMacrostates], dtype = float);
		self.fitnesses = numpy.ones([self.nPositions, 20], dtype = float);
		self.frequencies = numpy.zeros([self.nPositions, 20], dtype = float);

		# a little checker to prevent two identical entries, used in function addMacrostateData()
		for i in range(self.nMacrostates):
			for j in range(self.nPositions):
				self.macrostateResidueEnergies[j][0][i] = 65536.0; 

		if self.useMicrostateData:
			self.microstatesUsed = numpy.zeros([0]);
			self.microstateCounts = numpy.zeros([self.nPositions, 20], dtype = int);
			self.microstateResidueEnergies = numpy.zeros([self.nPositions, 20, self.nMacrostates, 512], dtype = float); # magic number 512 - max expected microstates
		

	def constructFromExisting(existing, ensembleSize:int, backrubTemp:float, boltzmannTemp:float, weights:numpy.array, steepness:float):
		"""
		"Overloaded" "constructor" that uses a pre-existing Model as a template

		@param existing			pre-exisiting Model object
		@param ensembleSize		int, new ensemble size
		@param backrubTemp		float, new backrub temperate
		@param boltzmannTemp	float, new Boltzmann temerature
		@param weights			float[], new weights
		@return Model
		"""
		# TODO: implement microstate construction
		new = Model(ensembleSize, backrubTemp, boltzmannTemp, weights, steepness, existing.nPositions, existing.positionOffset, existing.useMicrostateData);
		new.MACROSTATES = existing.MACROSTATES;
		new.nMacrostates = existing.nMacrostates;
		new.macrostatesUsed = existing.macrostatesUsed;
		if not existing.useMicrostateData:
			new.macrostateResidueEnergies = numpy.array(existing.macrostateResidueEnergies);
		else:
			new.microstateResidueEnergies = numpy.array(existing.microstateResidueEnergies);
		
		return new;

	def addMacrostateData(self, macrostate:int, position:int, energies:"float[]") -> None:
		"""
		Inserts a macrostate_position set of fitness values to this model

		@param macrostate		int, the macrostate this corresponds to
		@param position			int, the position the energies corresponds to
		@param energies			float[] of length-20 of the energies
		@return void
		"""
		if self.macrostateResidueEnergies[position - self.positionOffset][0][macrostate] != 65536.0:
			raise Exception("Something something this entry already full");

		for i in range(20):
			self.macrostateResidueEnergies[position - self.positionOffset][i][macrostate] = energies[i];

	def addMicrostateData(self, macrostate:int, position:int, energies:"float") -> None:
		"""
		Inserts a microstate_position fitness into this model

		@param macrostate		int, the macrostate this microstate belongs to
		@param position			int, the position the energy corresponds to
		@param energy			float, the energy of this mutation
		@retun void
		"""
		position -= self.positionOffset;

		# TODO: do I need a overwrite check as in adding macrostate data?
		for i in range(20):
			self.microstateResidueEnergies[position][i][macrostate][self.microstateCounts[position]];

		self.microstateCounts[position] += 1;

		return None;

	# change weights sets
	# TODO ascartain that this is actually necessary
	def setWeights(self, newWeights):
		self.weights = newWeights;

	# PRIVATE
	def averageMicrostates(self) -> None:
		"""
		Boltzmann averages the microstates to calculate the energy for the macrostate
		"""
		# pick backbones to use for the ensemble.
		self.microstatesUsed = numpy.zeros([self.nMacrostates, self.ensembleSize]);
		for i in range(self.nMacrostates):
			self.microstatesUsed[i] = numpy.random.randint(0, self.microstateCounts[i], [self.ensembleSize]);

		# cherry-pick out the selected microstates
		energies = numpy.zeros([self.nPositions, 20, self.nMacrostates, self.ensembleSize]);
		for i in range(self.nMacrostates):
			for j in range(self.ensembleSize):
				for k in range(self.nPositions):
					for l in range(20):
						energies[k][l][i][j] = self.microstateResidueEnergies[k][l][i][self.microstatesUsed[i][j]];

		# TODO: these fxns, copied from Smauel's email, calc only for a single residue
		
		if not self.useAltAveragingMethod:

			if (temp == 0.0):
				self.macrostateResidueEnergies = numpy.amin(energies, axis = 3);
			elif (temp == -1.0):
				self.macrostateResidueEnergies = numpy.mean(energies, axis = 3); 
			else:
				self.macrostateResidueEnergies = numpy.sum(energies * numpy.exp(energies / -temp), axis = 3) / numpy.sum(numpy.exp(energies / -temp), axis = 3);

		else:
			if (temp == 0.0):
				self.macrostateResidueEnergies = numpy.amin(energies, axis = 3)
			elif (temp == -1.0):
				self.macrostateResidueEnergies = numpy.mean(energies, axis = 3)
			else:
				self.macrostateResidueEnergies = -numpy.log(sum(numpy.exp(energies / -temp), axis = 3));

		# After averaging, delete the 4D array to save space and flip the microstate flag
		self.microstateResidueEnergies = numpy.array(0);
		self.useMicrostateData = False;
		return None;

	# PRIVATE
	def calcFitness(self):
		"""
		Calculates the fitnesses of the each residue at each position.
		There is no need for this function to be externally called

		@param void
		@return void
		"""

		# collapse microstates into macrostates
		if self.useMicrostateData:
			self.averageMicrostates();

		# TODO: verify implementation
		minEnergies = numpy.amin(self.macrostateResidueEnergies, axis = 1);	# for each position and macrostate, which residue had min energy?
		offsets = minEnergies + numpy.divide(numpy.log(99), self.steepness);		# calculate offset is double[position][macrostate]

		for i in range(self.nPositions):
			for j in range(20):
				for k in range(self.nMacrostates):
					f = 1.0 / (1.0 + numpy.exp(self.steepness * (self.macrostateResidueEnergies[i][j][k] - offsets[i][k])));
					self.fitnesses[i][j] *= (1 - self.weights[k] + self.weights[k] * f);

		# PRIVATE
	def calcFrequencies(self) -> None:
		"""
		Calculates the fitnesses and the frequencies of residues at each location

		@param void
		@return void
		"""
		if not self.isFrequenciesCalculated: # only do it once!
			self.isFrequenciesCalculated = True;
			self.calcFitness();
			# TODO: Verify implementation
			self.frequencies = numpy.divide(self.fitnesses, numpy.subtract(1.0, self.fitnesses));	# non-normalized frequencies
			sums = numpy.sum(self.frequencies, axis = 1);
			for i in range(self.nPositions):	# normalize
				self.frequencies[i] = numpy.divide(self.frequencies[i], sums[i]);

	# get functions
	# member fields should not be directly accessed; use these get funtions instead
	def getEnsembleSize(self):
		"""
		Self-explanatory name

		@return int
		"""
		return self.ensembleSize;

	def getBackrubTemp(self):
		"""
		Self-explanatory name

		@return double
		"""
		return self.backrubTemp;

	def getBoltzmannTemp(self):
		"""
		Self-explanatory name

		@return double
		"""
		return self.boltzmannTemp;

	def getWeights(self):
		"""
		Self-explanatory name

		@return float[]
		"""
		return self.weights;

	def getSteepness(self):
		"""
		Self-explanatory name

		@return float
		"""
		return self.steepness;

	def getFrequencies(self):
		"""
		Self-explanatory name

		@return float[][]
		"""
		
		# this is a special instance for storing data
		if self.ensembleSize == 0 and self.useMicrostateData:
			raise PermissionError("This object is a raw data storage instance and this call should not have been made");
			return None;

		if not self.isFrequenciesCalculated:
			self.calcFrequencies();
		return numpy.array(self.frequencies);

	# comparison operators based on similarity measure
	def __eq__(self, other):
		assert(isinstance(other, Model));
		return self.recovery == other.recovery;		# see, this is where static typing comes in useful. it lets autocomplete see the fields of other
	
	def __le__(self, other):
		assert(isinstance(other, Model));
		return self.recovery <= other.recovery;
	
	def __lt__(self, other):
		assert(isinstance(other, Model));
		return self.recovery < other.recovery;

	def __ge__(self, other):
		assert(isinstance(other, Model));
		return self.recovery >= other.recovery;

	def __gt__(self, other):
		assert(isinstance(other, Model));
		return self.recovery < other.recovery;

	def __ne__(self, other):
		assert(isinstance(other, Model));
		return self.recovery != other.recovery;