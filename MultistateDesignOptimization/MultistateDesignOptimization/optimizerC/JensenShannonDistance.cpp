#include "stdafx.h"
#include "SimilarityMeasure.h"
#include "JensenShannonDistance.h"

namespace OPTIMIZER
{
	JensenShannonDistance::JensenShannonDistance() : SimilarityMeasure() 
	{
		targetFrequencies->operator/=(accu(*targetFrequencies));
		hTargetFreqs = h(targetFrequencies);
	}

	JensenShannonDistance::JensenShannonDistance(const SimilarityMeasure &other) : SimilarityMeasure(other) 
	{
		targetFrequencies->operator/=(accu(*targetFrequencies));
		hTargetFreqs = h(targetFrequencies);
	}

	JensenShannonDistance::JensenShannonDistance(mat *targetFreqs) : SimilarityMeasure(targetFreqs) 
	{
		targetFrequencies->operator/=(accu(*targetFrequencies));
		hTargetFreqs = h(targetFrequencies);
	}

	double JensenShannonDistance::getSimilarity(mat *expFrequencies)
	{
		// TODO: check for correctness and memory leaks
		if (expFrequencies == NULL)
			throw std::runtime_error("Null matrix");

		mat efreqs = (*expFrequencies) / accu(*expFrequencies);
		mat calcs = *h(&efreqs);
		mat calcs2 = *targetFrequencies + efreqs;
		calcs2 = *h(&calcs2);
		//		h(expFreq) - h(exp + tgt) + h(tgt)
		calcs = calcs - calcs2 + *hTargetFreqs;
		calcs.elem(find_nonfinite(calcs)).zeros();	// zero the NaNs
		double JSdiv = sqrt(accu(calcs) * 0.5);

		if (JSdiv < 0 || JSdiv > 1)
			throw std::runtime_error("calculated value out of bounds");

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
