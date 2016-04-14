from Optimizer import Optimizer;
from SimilarityMeasure import SimilarityMeasure;
from SearchAlgorithm import SearchAlgorithm;
from CuckooSearch import CuckooSearch;
from KLDivergence import KLDivergence;
from JensenShannonDistance import JensenShannonDistance;
from CosineSimilarity import CosineSimilarity;
from CosineSimilarityByPosition import CosineSimilarityByPosition;
from EntropyWeightsMixedSimilarity import EntropyWeightsMixedSimilarity;
from model import Model;
from enumeration import enum;
from datetime import *
import numpy;
import threading
from io import *

def evalCosMatchWithJSD():
	MACROSTATES = enum("E-DHF-NADPH", "E-NADPH", "E-OPEN", "E-THF", "E-THF-NADPX", "TS");
	RESIDUES = enum('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');

	# only looking at MACROSTATE.TS
	# only optimizing backrub temperature and steepness
	ensembleSizes = numpy.array([50]);
	backrubTemps = numpy.array([1.8]);
	boltzmannTemps = numpy.array([-1.0]);
	steepnessRange = numpy.array([2.0541]);
	minWeights = numpy.array([0, 0, 0, 0, 0, 1]);
	maxWeights = numpy.array([0, 0, 0, 0, 0, 1]);


	targetFreqs = "C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\ecDHFR_openseq_bacterial_representative_final_align_trim.fasta";
	data = "C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\DHFR_MSD_M20loop\\DHFR_MSD_M20loop_repeat7.tsv";
	dataMicro = "C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\20160120_M20_enumeration_scores\\20160120_M20_enumeration_scores.tsv";
	
	optimizer = Optimizer(MACROSTATES);

	# slightly different paths on my two computers
	try:
		optimizer.readTargetFrequencies(targetFreqs);	
		optimizer.readData(data);
	except:
		optimizer.readTargetFrequencies(targetFreqsAlt);	
		optimizer.readData(dataAlt);

	print("Files read in");
	
	model = Model.constructFromExisting(optimizer.getModelByParams(backrubTemps[0], ensembleSizes[0], boltzmannTemps[0]), ensembleSizes[0], backrubTemps[0], boltzmannTemps[0], numpy.array([0, 0, 0, 0, 0, 1]), steepnessRange[0]);

	JSD = JensenShannonDistance(optimizer.targetFrequencies);
	cos = CosineSimilarity(optimizer.targetFrequencies);

	now = datetime.now();
	optimizer.writeFrequenciesToFASTA(model.getFrequencies(), " 7jsd vals direct freqs" + now.strftime('%Y%m%d%H%M') + ".fasta");
	print(JSD.getSimilarityMeasure(model.getFrequencies()));
	print(cos.getSimilarityMeasure(model.getFrequencies()));

	s = CuckooSearch(optimizer.models, CosineSimilarity(optimizer.targetFrequencies), False, 1, 1, 0.25);
	#s.setAllSearchToFalse();
	s.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
	s.setSearchParameters(False, False, False, False, numpy.array([False, False, False, False, False, False]));
	s.setMaxIterations(1);
	s.suppressOutputs = True;
	optimizer.useAlgorithm(s);
	optimizer.optimize();
	#optimizer.writeFrequenciesToFASTA(optimizer.getBestFrequencies(), " 1cos vals optim freqs" + now.strftime('%Y%m%d%H%M') + ".fasta");
	print(optimizer.getBestParameters()['match']);
	
	return None;

def microstates():
	MACROSTATES = enum("E-DHF-NADPH", "E-NADPH", "E-THF", "E-THF-NADPX", "TS");
	RESIDUES = enum('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');

	ensembleSizes = numpy.array([32, 64, 128, 256]);
	backrubTemps = numpy.array([0.3, 0.6, 0.9, 1.2, 1.5, 1.8]);
	boltzmannTemps = numpy.array([-1, 5]);
	steepnessRange = numpy.array([1, 7]);
	minWeights = numpy.array([0, 0, 0, 0, 0]);
	maxWeights = numpy.array([1, 1, 1, 1, 1]);

	data = "/netapp/home/tianjiao.zhang/data/microstates.dat";
	targetFreqs = "/netapp/home/tianjiao.zhang/data/ecDHFR_openseq_bacterial_representative_final_align_trim.fasta";
	
	optimizer = Optimizer(MACROSTATES, True);

	optimizer.readTargetFrequencies(targetFreqs);	
	optimizer.readFormattedMicrostateData(data);
		
	search = CuckooSearch(optimizer.models, JensenShannonDistance(optimizer.targetFrequencies), True, 128, 1, 0.25);
	search.setMaxIterations(4096);
	search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
	search.setAllSearchToTrue();
	optimizer.useAlgorithm(search);
	optimizer.optimize();
	now = datetime.now();
	optimizer.writeFrequenciesToFASTA(optimizer.getBestFrequencies(), "var microstates all states" + now.strftime('%Y%m%d%H%M') + ".fasta");
	optimizer.writeBestParamsToText("var Microstates all states" + now.strftime('%Y%m%d%H%M'));

def TSonlyMacroAltCos(i:int):
	MACROSTATES = enum("E-DHF-NADPH", "E-NADPH", "E-OPEN", "E-THF", "E-THF-NADPX", "TS");
	RESIDUES = enum('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');

	ensembleSizes = numpy.array([20, 50]);
	backrubTemps = numpy.array([0.3, 0.6, 0.9, 1.2, 1.5, 1.8]);
	boltzmannTemps = numpy.array([0, -1, 1, 5.0]);
	steepnessRange = numpy.array([0.5, 5]);
	minWeights = numpy.array([0, 0, 0, 0, 0, 1]);
	maxWeights = numpy.array([0, 0, 0, 0, 0, 1]);

	data = "/netapp/home/tianjiao.zhang/data/DHFR_MSD_M20loop_repeat" + str(i + 1) + ".tsv";
	targetFreqs = "/netapp/home/tianjiao.zhang/data/ecDHFR_openseq_bacterial_representative_final_align_trim.fasta";
	optimizer = Optimizer(MACROSTATES);
	optimizer.readTargetFrequencies(targetFreqs);
	optimizer.readData(data);

	search = CuckooSearch(optimizer.models, CosineSimilarityByPosition(optimizer.targetFrequencies), False, 64, 1, 0.25);
	measure = "cos alt";
	search.setMaxIterations(2048);
	search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
	search.setSearchParameters(True, True, True, True, numpy.array([False, False, False, False, False, True]));
	optimizer.useAlgorithm(search);
	optimizer.optimize();

	name = "Macrostates " + str(i + 1) + measure + datetime.now().strftime('%Y%m%d%H%M');
	optimizer.writeFrequenciesToFASTA(optimizer.getBestFrequencies(), name + ".fasta", 3);
	optimizer.writeBestParamsToText(name + ".txt");
	
	return None;