#include "stdafx.h"
#include "SimilarityMeasure.h"
#include "JensenShannonDistance.h"

namespace OPTIMIZER
{
	JensenShannonDistance::JensenShannonDistance() : SimilarityMeasure() 
	{
		targetFrequencies->operator/=(sum(*targetFrequencies));
		hTargetFreqs = h(targetFrequencies);
	}

	JensenShannonDistance::JensenShannonDistance(const SimilarityMeasure &other) : SimilarityMeasure(other) 
	{
		targetFrequencies->operator/=(sum(*targetFrequencies));
		hTargetFreqs = h(targetFrequencies);
	}

	JensenShannonDistance::JensenShannonDistance(mat *targetFreqs) : SimilarityMeasure(targetFreqs) 
	{
		targetFrequencies->operator/=(sum(*targetFrequencies));
		hTargetFreqs = h(targetFrequencies);
	}

	double JensenShannonDistance::getSimilarity(mat *expFrequencies)
	{
		// TODO: check for correctness and memory leaks
		mat calcs = *h(expFrequencies);
		mat calcs2 = *targetFrequencies + *expFrequencies;
		calcs2 = *h(&calcs2);
		//		h(expFreq) - h(exp + tgt) + h(tgt)
		calcs = calcs - calcs2 + *hTargetFreqs;
		calcs.elem(find_nonfinite(calcs)).zeros();	// zero the NaNs
		double JSdiv = accu(calcs);

		if (JSdiv < 0 || JSdiv > 1)
			throw exception::exception("calculated value out of bounds");

		delete expFrequencies;

		return JSdiv;
	}

	mat* JensenShannonDistance::h(mat *freqs)
	{
		// TODO: check for memory leaks
		mat *out = new mat(*freqs % arma::log2(*freqs));
		out->operator*=(-1);
		return out;
	}

	string *JensenShannonDistance::toString()
	{
		return new string("Jensen-Shannon Distance");
	}
}