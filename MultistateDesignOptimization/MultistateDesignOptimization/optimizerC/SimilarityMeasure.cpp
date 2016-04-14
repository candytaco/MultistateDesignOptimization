#include "SimilarityMeasure.h"

namespace OPTIMIZER
{
	SimilarityMeasure::SimilarityMeasure()
	{
		targetFrequencies = NULL;
	}

	SimilarityMeasure::SimilarityMeasure(SimilarityMeasure *existing)
	{
		this->targetFrequencies = mxDuplicateArray(existing->targetFrequencies);
	}

	SimilarityMeasure::SimilarityMeasure(mxArray *targetFreqs)
	{
		this->targetFrequencies = targetFreqs;
	}

	void SimilarityMeasure::setTargetFrequencies(mxArray *newTargetFreqs)
	{
		this->targetFrequencies = newTargetFreqs;
	}

	SimilarityMeasure::~SimilarityMeasure()
	{
		if (this->targetFrequencies != NULL)
			mxDestroyArray(this->targetFrequencies);
	}
}