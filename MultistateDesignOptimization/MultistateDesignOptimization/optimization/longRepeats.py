from Optimizer import Optimizer;
from SimilarityMeasure import SimilarityMeasure;
from SearchAlgorithm import SearchAlgorithm;
from CuckooSearch import CuckooSearch;
from KLDivergence import KLDivergence;
from JensenShannonDistance import JensenShannonDistance;
from CosineSimilarity import CosineSimilarity;
from EntropyWeightsMixedSimilarity import EntropyWeightsMixedSimilarity;
from model import Model;
from enumeration import enum;
from datetime import *
import numpy;
from io import *

# run repeats of the same fitting algorithm over long generations
data = "/netapp/home/tianjiao.zhang/data/microstates.dat";
targetFreqs = "/netapp/home/tianjiao.zhang/data/ecDHFR_openseq_bacterial_representative_final_align_trim.fasta";

MACROSTATES = enum("E-DHF-NADPH", "E-NADPH", "E-THF", "E-THF-NADPX", "TS");
RESIDUES = enum('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y');

ensembleSizes = numpy.array([16, 32, 64, 128, 256]);
backrubTemps = numpy.array([0.3, 0.6, 0.9, 1.2, 1.5, 1.8]);
boltzmannTemps = numpy.array([-1, 5]);
steepnessRange = numpy.array([1, 5]);
minWeights = numpy.array([0, 0, 0, 0, 0]);
maxWeights = numpy.array([1, 1, 1, 1, 1]);
	
optimizer = Optimizer(MACROSTATES, True);
optimizer.readTargetFrequencies(targetFreqs);	
optimizer.readFormattedMicrostateData(data);
		
search = CuckooSearch(optimizer.models, JensenShannonDistance(optimizer.targetFrequencies), True, 128, 1.25, 0.25);
search.setMaxIterations(6000);
search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
search.setAllSearchToTrue();
search.suppressOutputs = True;
optimizer.useAlgorithm(search);
optimizer.optimize();
now = datetime.now();
optimizer.writeFrequenciesToFASTA(optimizer.getBestFrequencies(), "microstateRepeats" + now.strftime('%Y%m%d%H%M') + ".fasta");
optimizer.writeBestParamsToText("microstateRepeats" + now.strftime('%Y%m%d%H%M'));