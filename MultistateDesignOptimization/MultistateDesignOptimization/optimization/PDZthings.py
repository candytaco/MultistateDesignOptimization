from Optimizer import Optimizer
from model import Model
from SimilarityMeasure import SimilarityMeasure
from SearchAlgorithm import SearchAlgorithm
from CuckooSearch import CuckooSearch
from KLDivergence import KLDivergence
from JensenShannonDistance import JensenShannonDistance
from CosineSimilarity import CosineSimilarity
from EntropyWeightsMixedSimilarity import EntropyWeightsMixedSimilarity
from EntropyWeightedSimilarity import EntropyWeightedSimilarity
from ChunkedSimilarity import ChunkedSimilarity
from JSDistByPos import JSDistByPos
from enumeration import enum
from datetime import *
import numpy
from io import *

# running things on PDZ data

def PDZsmallTest():

	data = "C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\PDZ\\PDZ1all.dat";
	targetFreqs = "C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\PDZ\\PF00595_80.fasta";
	positionInd = "C:\\Users\\Candy\\SkyDrive\\Documents\\rotation 2\\PDZ\\1BE9.INDICES_MSA";

	dataAlt = "C:\\Users\\Candy_000\\SkyDrive\\Documents\\rotation 2\\PDZ\\PDZ1all.dat";
	targetFreqsAlt = "C:\\Users\\Candy_000\\SkyDrive\\Documents\\rotation 2\\PDZ\\PF00595_80.fasta";
	positionIndAlt = "C:\\Users\\Candy_000\\SkyDrive\\Documents\\rotation 2\\PDZ\\1BE9.INDICES_MSA";

	macrostates = enum("Stability", "Affinity");
	RESIDUES = enum('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');

	optimizer = Optimizer(macrostates, True, False);

	try:
		optimizer.readTargetFrequencies(targetFreqs, positionInd);
		optimizer.readFormattedMicrostateData(data);
	except FileNotFoundError:
		optimizer.readTargetFrequencies(targetFreqsAlt, positionIndAlt);
		optimizer.readFormattedMicrostateData(dataAlt);

	ensembleSizes = numpy.array([20]);
	backrubTemps = numpy.array([0.3, 0.6, 0.9, 1.2, 1.5, 1.8]);
	boltzmannTemps = numpy.array([-1, 5.0]);
	steepnessRange = numpy.array([0.5, 5]);
	minWeights = numpy.array([0, 0]);
	maxWeights = numpy.array([1, 1]);

	search = CuckooSearch(optimizer.models, JensenShannonDistance(optimizer.targetFrequencies), True, 16, 1, 0.25);
	search.setMaxIterations(128);
	search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
	search.setSearchParameters(False, True, True, True, numpy.array([True, True]));
	optimizer.useAlgorithm(search);

	print("\nJS Dist");
	#search.setSimilarityMeasure(JensenShannonDistance(optimizer.targetFrequencies));
	optimizer.optimize();
	now = datetime.now();
	optimizer.writeFrequenciesToFASTA(optimizer.getBestFrequencies(), "PDZ test " + now.strftime('%Y%m%d%H%M%S') + ".fasta");
	optimizer.writeBestParamsToText("PDZ test " + now.strftime('%Y%m%d%H%M%S'));
	print(optimizer.getBestParameters()['match']);

	print("Bloop");
	print(optimizer.nPositions);
	return None;

def PDZVarGens(generations:int = 128):

	data = "/netapp/home/tianjiao.zhang/data/PDZselected.dat";
	targetFreqs = "/netapp/home/tianjiao.zhang/data/PDZselected.fasta";
	positionInd = "/netapp/home/tianjiao.zhang/data/1BE9.INDICES_MSA";

	macrostates = enum("Stability", "Affinity");
	RESIDUES = enum('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');

	optimizer = Optimizer(macrostates, True);

	optimizer.readTargetFrequencies(targetFreqs);
	optimizer.readFormattedMicrostateData(data);

	ensembleSizes = numpy.array([25]);
	backrubTemps = numpy.array([0.3, 0.6, 0.9, 1.2, 1.5, 1.8]);
	boltzmannTemps = numpy.array([-1, 5.0]);
	steepnessRange = numpy.array([0.5, 5]);
	minWeights = numpy.array([0, 0]);
	maxWeights = numpy.array([1, 1]);

	search = CuckooSearch(optimizer.models, JSDistByPos(optimizer.targetFrequencies), True, 32, 1, 0.25);
	search.setMaxIterations(generations);
	search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
	search.setSearchParameters(False, True, True, True, numpy.array([True, True]));
	search.suppressOutputs = True;
	optimizer.useAlgorithm(search);

	optimizer.optimize();
	now = datetime.now();
	optimizer.writeFrequenciesToFASTA(optimizer.getBestFrequencies(), "PDZ varGens " + str(generations) + " " + now.strftime('%Y%m%d%H%M%S') + ".fasta");
	optimizer.writeBestParamsToText("PDZ varGens "  + str(generations) + " " + now.strftime('%Y%m%d%H%M%S'));

	return None;

def PDZVarGensCos(generations:int = 128):

	data = "/netapp/home/tianjiao.zhang/data/PDZselected.dat";
	targetFreqs = "/netapp/home/tianjiao.zhang/data/PDZselected.fasta";
	positionInd = "/netapp/home/tianjiao.zhang/data/1BE9.INDICES_MSA";

	macrostates = enum("Stability", "Affinity");
	RESIDUES = enum('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');

	optimizer = Optimizer(macrostates, True);

	optimizer.readTargetFrequencies(targetFreqs);
	optimizer.readFormattedMicrostateData(data);

	ensembleSizes = numpy.array([25]);
	backrubTemps = numpy.array([0.3, 0.6, 0.9, 1.2, 1.5, 1.8]);
	boltzmannTemps = numpy.array([-1, 5.0]);
	steepnessRange = numpy.array([0.5, 5]);
	minWeights = numpy.array([0, 0]);
	maxWeights = numpy.array([1, 1]);

	search = CuckooSearch(optimizer.models, CosineSimilarity(optimizer.targetFrequencies), True, 32, 1, 0.25);
	search.setMaxIterations(generations);
	search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
	search.setSearchParameters(False, True, True, True, numpy.array([True, True]));
	search.suppressOutputs = True;
	optimizer.useAlgorithm(search);

	optimizer.optimize();
	now = datetime.now();
	optimizer.writeFrequenciesToFASTA(optimizer.getBestFrequencies(), "PDZ cosine " + str(generations) + " " + now.strftime('%Y%m%d%H%M%S') + ".fasta");
	optimizer.writeBestParamsToText("PDZ cosine "  + str(generations) + " " + now.strftime('%Y%m%d%H%M%S'));

	return None;

def PDZVarGensWeightedJSD(generations:int = 128):

	data = "/netapp/home/tianjiao.zhang/data/PDZselected.dat";
	targetFreqs = "/netapp/home/tianjiao.zhang/data/PDZselected.fasta";
	positionInd = "/netapp/home/tianjiao.zhang/data/1BE9.INDICES_MSA";

	macrostates = enum("Stability", "Affinity");
	RESIDUES = enum('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');

	optimizer = Optimizer(macrostates, True);

	optimizer.readTargetFrequencies(targetFreqs);
	optimizer.readFormattedMicrostateData(data);

	ensembleSizes = numpy.array([25, 50, 100]);
	backrubTemps = numpy.array([0.3, 0.6, 0.9, 1.2, 1.5, 1.8]);
	boltzmannTemps = numpy.array([-1, 5.0]);
	steepnessRange = numpy.array([0.5, 5]);
	minWeights = numpy.array([0, 0]);
	maxWeights = numpy.array([1, 1]);

	search = CuckooSearch(optimizer.models, EntropyWeightedSimilarity(JensenShannonDistance(), optimizer.targetFrequencies), True, 32, 1, 0.25);
	search.setMaxIterations(generations);
	search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
	search.setAllSearchToTrue();
	optimizer.useAlgorithm(search);

	optimizer.optimize();
	now = datetime.now();
	optimizer.writeFrequenciesToFASTA(optimizer.getBestFrequencies(), "PDZ weighted JSD " + str(generations) + " " + now.strftime('%Y%m%d%H%M%S') + ".fasta");
	optimizer.writeBestParamsToText("PDZ weighted JSD "  + str(generations) + " " + now.strftime('%Y%m%d%H%M%S'));

	return None;

def PDZVarGensChunkedJSD(generations:int = 128):

	data = "/netapp/home/tianjiao.zhang/data/PDZselected.dat";
	targetFreqs = "/netapp/home/tianjiao.zhang/data/PDZselected.fasta";
	positionInd = "/netapp/home/tianjiao.zhang/data/1BE9.INDICES_MSA";

	macrostates = enum("Stability", "Affinity");
	RESIDUES = enum('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');

	optimizer = Optimizer(macrostates, True);

	optimizer.readTargetFrequencies(targetFreqs);
	optimizer.readFormattedMicrostateData(data);

	ensembleSizes = numpy.array([25]);
	backrubTemps = numpy.array([0.3, 0.6, 0.9, 1.2, 1.5, 1.8]);
	boltzmannTemps = numpy.array([-1, 5.0]);
	steepnessRange = numpy.array([0.5, 5]);
	minWeights = numpy.array([0, 0]);
	maxWeights = numpy.array([1, 1]);

	search = CuckooSearch(optimizer.models, ChunkedSimilarity(JensenShannonDistance(), 8, optimizer.targetFrequencies), True, 32, 1, 0.25);
	search.setMaxIterations(generations);
	search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
	search.setSearchParameters(False, True, True, True, numpy.array([True, True]));
	search.suppressOutputs = True;
	optimizer.useAlgorithm(search);

	optimizer.optimize();
	now = datetime.now();
	optimizer.writeFrequenciesToFASTA(optimizer.getBestFrequencies(), "PDZ chunked JSD " + str(generations) + " " + now.strftime('%Y%m%d%H%M%S') + ".fasta");
	optimizer.writeBestParamsToText("PDZ chunked JSD "  + str(generations) + " " + now.strftime('%Y%m%d%H%M%S'));

	return None;

def PDZVarGensWeightedJSDAll(generations:int = 128):

	data = "/netapp/home/tianjiao.zhang/data/PDZ1all.dat";
	targetFreqs = "/netapp/home/tianjiao.zhang/data/PDZ.fasta";
	positionInd = "/netapp/home/tianjiao.zhang/data/1BE9.INDICES_MSA";

	macrostates = enum("Stability", "Affinity");
	RESIDUES = enum('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');

	optimizer = Optimizer(macrostates, True, False);

	optimizer.readTargetFrequencies(targetFreqs, positionInd);
	optimizer.readFormattedMicrostateData(data);

	ensembleSizes = numpy.array([25]);
	backrubTemps = numpy.array([0.3, 0.6, 0.9, 1.2, 1.5, 1.8]);
	boltzmannTemps = numpy.array([-1, 5.0]);
	steepnessRange = numpy.array([0.5, 5]);
	minWeights = numpy.array([0, 0]);
	maxWeights = numpy.array([1, 1]);

	search = CuckooSearch(optimizer.models, EntropyWeightedSimilarity(JensenShannonDistance(), optimizer.targetFrequencies), True, 32, 1, 0.25);
	search.setMaxIterations(generations);
	search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
	search.setSearchParameters(False, True, True, True, numpy.array([True, True]));
	search.suppressOutputs = True;
	optimizer.useAlgorithm(search);

	optimizer.optimize();
	now = datetime.now();
	optimizer.writeFrequenciesToFASTA(optimizer.getBestFrequencies(), "PDZ weighted JSD " + str(generations) + " " + now.strftime('%Y%m%d%H%M%S') + ".fasta");
	optimizer.writeBestParamsToText("PDZ weighted JSD "  + str(generations) + " " + now.strftime('%Y%m%d%H%M%S'));

	return None;

def PDZVarGensWeightedJSDGDeweighted(generations:int = 128):

	data = "/netapp/home/tianjiao.zhang/data/PDZ1all G deweighted.dat";
	targetFreqs = "/netapp/home/tianjiao.zhang/data/PDZ.fasta";
	positionInd = "/netapp/home/tianjiao.zhang/data/1BE9.INDICES_MSA";

	macrostates = enum("Stability", "Affinity");
	RESIDUES = enum('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');

	optimizer = Optimizer(macrostates, True, False);

	optimizer.readTargetFrequencies(targetFreqs, positionInd);
	optimizer.readFormattedMicrostateData(data);

	ensembleSizes = numpy.array([25, 50, 100]);
	backrubTemps = numpy.array([0.3, 0.6, 0.9, 1.2, 1.5, 1.8]);
	boltzmannTemps = numpy.array([-1, 5.0]);
	steepnessRange = numpy.array([0.5, 5]);
	minWeights = numpy.array([0, 0]);
	maxWeights = numpy.array([1, 1]);

	search = CuckooSearch(optimizer.models, EntropyWeightedSimilarity(JensenShannonDistance(), optimizer.targetFrequencies), True, 32, 1, 0.25);
	search.setMaxIterations(generations);
	search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
	search.setAllSearchToTrue();
	optimizer.useAlgorithm(search);

	optimizer.optimize();
	now = datetime.now();
	optimizer.writeFrequenciesToFASTA(optimizer.getBestFrequencies(), "PDZ G deweighted JSD " + str(generations) + " " + now.strftime('%Y%m%d%H%M%S') + ".fasta");
	optimizer.writeBestParamsToText("PDZ G deweighted JSD "  + str(generations) + " " + now.strftime('%Y%m%d%H%M%S'));

	return None;