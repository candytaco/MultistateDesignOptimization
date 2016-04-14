from model import Model
from SearchAlgorithm import SearchAlgorithm
from SimilarityMeasure import SimilarityMeasure
from io import *
from enumeration import enum
import numpy
import ast
import datetime

# THESE TWO FUNCTIONS HAVE BEEN ADDED AS STATIC FXNS IN CLASS OPTIMIZER
#def positionReindexer(indices:str) -> {}:
#	"""
#	Used to convert arbitrary positions to 0-based indices since the positions
#	we look at may not start at zero
#	"""
#	infile = open(indices, 'r');
#	indices = {};
#	# the two indices at the two terminii are thrown out for alignment
#	# but we need these to read all data
#	n = 0;
#	for line in infile:
#		entries = line.split(' ');
#		i = entries[2].strip('\n');
#		indices[i] = n;
#		n += 1;
#	indices['nPos'] = n;
#	infile.close();
#	return indices;

#def convertMicrostateDataIndexCorrection(source:str, out:str, indexer:str, entries:int, minPosition:int, nPositions:int, stateIndex:dict = None) -> None:
#	"""
#	Serializes a large tsv to a formatted dat and re-indixes positions
#	"""
#	indexToRes = {0:'A', 1:'C', 2:'D', 3:'E', 4:'F', 5:'G', 6:'H', 7:'I', 8:'K', 9:'L', 10:'M', 11:'N', 12:'P', 13:'Q', 14:'R', 15:'S', 16:'T', 17:'V', 18:'W', 19:'Y'};
#	if stateIndex is None:
#		macStateToIndex = {"Stability":0, "Affinity":1};
#	else:
#		macStateToIndex = stateIndex;

#	indexCorrection = positionReindexer(indexer);
#	infile = open(source, 'r');
#	outfile = open(out, 'w');

#	outfile.write("{:d} ".format(minPosition));
#	outfile.write("{:d} ".format(indexCorrection['nPos']));
#	outfile.write("{:d}\n".format(entries));

#	n = 0;
#	#line = infile.readline();
#	#while true:
#	for line in infile:
#		if numpy.mod(n, 1000) == 0:
#			print(n);
#		n += 1;
#		entries = line.split('\t');
#		macrostate = entries[0];
#		backrubT = entries[1];
#		try:
#			position = indexCorrection[entries[2]];
#		except KeyError:	# if the key wasn't found in the dict, it's one of the discarded positions
#			continue;		# skip it
#		#backbone = entries[3];
#		energies = ast.literal_eval(entries[4]);

#		# convert from string to useful data types
#		backrubT = float(backrubT);
#		macrostate = macStateToIndex[macrostate];
#		position = int(position);
#		temp = numpy.zeros([20]);
#		for i in range(20):
#			temp[i] = energies[indexToRes[i]];
#		energies = temp;

#		outfile.write("{:d} ".format(macrostate));
#		outfile.write("{:.1f} ".format(backrubT));
#		outfile.write("{:d} ".format(position));
#		for i in range(20):
#			outfile.write("{:.6f} ".format(energies[i]));
#		outfile.write('\n');

#	infile.close();
#	outfile.close();
#	return None;