#include "stdafx.h"
#include "SimilarityMeasure.h"

namespace OPTIMIZER
{
	SimilarityMeasure::SimilarityMeasure()
		: targetFrequencies(NULL)
	{
	}

	SimilarityMeasure::SimilarityMeasure(const SimilarityMeasure& other)
	{
		this->targetFrequencies = other.targetFrequencies; // shallow copy ok?
	}

	SimilarityMeasure::SimilarityMeasure(mat *targetFreqs)
		: targetFrequencies(targetFreqs)
	{
	}

	void SimilarityMeasure::setTargetFrequencies(mat *newTargetFreqs)
	{
		this->targetFrequencies = newTargetFreqs;
	}

	SimilarityMeasure::~SimilarityMeasure()
	{
		if (this->targetFrequencies != NULL)
			this->targetFrequencies->~Mat();
	}
}