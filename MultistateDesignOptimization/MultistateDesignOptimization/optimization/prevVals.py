from Optimizer import Optimizer;
from SimilarityMeasure import SimilarityMeasure;
from SearchAlgorithm import SearchAlgorithm;
from CuckooSearch import CuckooSearch;
from KLDivergence import KLDivergence;
from JensenShannonDistance import JensenShannonDistance;
from CosineSimilarity import CosineSimilarity;
from EntropyWeightsMixedSimilarity import EntropyWeightsMixedSimilarity;
from enumeration import enum;
from datetime import *
import numpy;
from io import *

def readParams(source:str) -> {}:
	infile = open(source, 'r');

	params = {};

	line = infile.readline();
	entries = line.split(' ');
	params['ensembleSize'] = numpy.array([int(entries[0])]);
	params['backrubT'] = numpy.array([float(entries[1])]);
	params['boltzmannT'] = numpy.array([float(entries[2])]);
	params['steepness'] = numpy.array([float(entries[3]), float(entries[3])]);

	line = infile.readline();
	entries = line.split(' ');
	arr = numpy.zeros([6]);
	for i in range(len(entries)):
		arr[i] = float(entries[i]);
	params['weights'] = arr;

	return params;

def prevOptValsTest():
	MACROSTATES = enum("E-DHF-NADPH", "E-NADPH", "E-OPEN", "E-THF", "E-THF-NADPX", "TS");
	RESIDUES = enum('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');
	
	optimParams5 = readParams('C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\params5.txt');
	optimParams10 = readParams('C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\params10.txt');
	optimParams20 = readParams('C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\params20.txt');
	optimParams50 = readParams('C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\params50.txt');

	targetFreqs = "C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\ecDHFR_openseq_bacterial_representative_final_align_trim.fasta";

	for i in range(1, 11):
		data = "C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\DHFR_MSD_M20loop\\DHFR_MSD_M20loop_repeat" + str(i) + ".tsv";
		optimizer = Optimizer(MACROSTATES);
		optimizer.readTargetFrequencies(targetFreqs);	
		optimizer.readData(data);

		ensembleSizes = optimParams20['ensembleSize'];
		backrubTemps = optimParams20['backrubT'];
		boltzmannTemps = optimParams20['boltzmannT'];
		steepnessRange = optimParams20['steepness'];
		minWeights = optimParams20['weights'];
		maxWeights = optimParams20['weights'];

		search = CuckooSearch(optimizer.models, JensenShannonDistance(optimizer.targetFrequencies), False, 1, 1, 0.25);
		search.setMaxIterations(1);
		search.suppressOutputs = True;
		search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
		search.setSearchParameters(False, False, False, False, numpy.array([False, False, False, False, False, False]));
		optimizer.useAlgorithm(search);

		optimizer.optimize();
		optimizer.writeFrequenciesToFASTA(optimizer.getBestFrequencies(), "Prev optim vals ensemble size 5 on macrostate " + str(i) + ".fasta");
		optimizer.writeBestParamsToText("Prev optim vals ensemble size 5 on macrostate " + str(i));
		
		ensembleSizes = optimParams50['ensembleSize'];
		backrubTemps = optimParams50['backrubT'];
		boltzmannTemps = optimParams50['boltzmannT'];
		steepnessRange = optimParams50['steepness'];
		minWeights = optimParams50['weights'];
		maxWeights = optimParams50['weights'];
		search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
		optimizer.optimize();
		optimizer.writeFrequenciesToFASTA(optimizer.getBestFrequencies(), "Prev optim vals ensemble size 50 on macrostate " + str(i) + ".fasta");
		optimizer.writeBestParamsToText("Prev optim vals ensemble size 50 on macrostate " + str(i));

def prevOptValsMicrostates():
	MACROSTATES = enum("E-DHF-NADPH", "E-NADPH", "E-THF", "E-THF-NADPX", "TS");
	RESIDUES = enum('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');
	
	optimParams20 = readParams('C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\params20.txt');
	optimParams50 = readParams('C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\params50.txt');

	data = "C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\20160120_M20_enumeration_scores\\microstates.dat";
	targetFreqs = "C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\ecDHFR_openseq_bacterial_representative_final_align_trim.fasta";

	optimizer = Optimizer(MACROSTATES, True);
	optimizer.readTargetFrequencies(targetFreqs);	
	optimizer.readFormattedMicrostateData(data);

	ensembleSizes = optimParams20['ensembleSize'];
	backrubTemps = optimParams20['backrubT'];
	boltzmannTemps = optimParams20['boltzmannT'];
	steepnessRange = optimParams20['steepness'];
	minWeights = optimParams20['weights'];
	maxWeights = optimParams20['weights'];

	
	search = CuckooSearch(optimizer.models, JensenShannonDistance(optimizer.targetFrequencies), False, 1, 1, 0.25);
	search.setMaxIterations(1);
	search.suppressOutputs = True;
	search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
	search.setSearchParameters(False, False, False, False, numpy.array([False, False, False, False, False, False]));
	optimizer.useAlgorithm(search);

	for i in range(16):
		optimizer.optimize();
		optimizer.writeFrequenciesToFASTA(optimizer.getBestFrequencies(), "Prev optim vals ensemble size 20 on microstates rep " + str(i) + ".fasta");
		optimizer.writeBestParamsToText("Prev optim vals ensemble size 20 on microstates rep " + str(i));

	ensembleSizes = optimParams50['ensembleSize'];
	backrubTemps = optimParams50['backrubT'];
	boltzmannTemps = optimParams50['boltzmannT'];
	steepnessRange = optimParams50['steepness'];
	minWeights = optimParams50['weights'];
	maxWeights = optimParams50['weights'];
	search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);

	for i in range(16):
		optimizer.optimize();
		optimizer.writeFrequenciesToFASTA(optimizer.getBestFrequencies(), "Prev optim vals ensemble size 50 on microstates rep " + str(i) + ".fasta");
		optimizer.writeBestParamsToText("Prev optim vals ensemble size 50 on microstates rep " + str(i));