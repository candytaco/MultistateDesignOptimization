from SimilarityMeasure import SimilarityMeasure
from enumeration import enum
from JensenShannonDistance import JensenShannonDistance
from CosineSimilarity import CosineSimilarity
from KLDivergence import KLDivergence
from EntropyWeightedSimilarity import EntropyWeightedSimilarity
from EntropyWeightsMixedSimilarity import EntropyWeightsMixedSimilarity
from Chi2Kernel import Chi2Kernel
from Optimizer import Optimizer
from CuckooSearch import CuckooSearch
from datetime import *
import sys
import numpy

data = "/netapp/home/tianjiao.zhang/data/PDZselected.dat";
targetFreqs = "/netapp/home/tianjiao.zhang/data/PDZselected.fasta";

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

search = CuckooSearch(optimizer.models, JensenShannonDistance(optimizer.targetFrequencies), True, 32, 1, 0.25);
search.setMaxIterations(4096);
search.setParamBounds(ensembleSizes, backrubTemps, boltzmannTemps, steepnessRange, minWeights, maxWeights);
search.setAllSearchToTrue();
optimizer.useAlgorithm(search);

optimizer.optimize();
now = datetime.now();
optimizer.writeFrequenciesToFASTA(optimizer.getBestFrequencies(), "PDZ cosine " + now.strftime('%Y%m%d%H%M%S') + ".fasta");
optimizer.writeBestParamsToText("PDZ cosine " + now.strftime('%Y%m%d%H%M%S'));
