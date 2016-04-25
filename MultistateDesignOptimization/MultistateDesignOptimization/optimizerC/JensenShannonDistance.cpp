#include "JensenShannonDistance.h"

namespace OPTIMIZER
{
	JensenShannonDistance::JensenShannonDistance() : SimilarityMeasure() 
	{
		targetFrequencies->operator/=(sum(*targetFrequencies));
	}

	JensenShannonDistance::JensenShannonDistance(const SimilarityMeasure &other) : SimilarityMeasure(other) 
	{
		targetFrequencies->operator/=(sum(*targetFrequencies));
	}

	JensenShannonDistance::JensenShannonDistance(mat *targetFreqs) : SimilarityMeasure(targetFreqs) 
	{
		targetFrequencies->operator/=(sum(*targetFrequencies));
	}

	double JensenShannonDistance::getSimilarity(mat *expFrequencies)
	{
		// TODO: implement this
		return -1;
	}

	string *JensenShannonDistance::toString()
	{
		return new string("Jensen-Shannon Distance");
	}
}