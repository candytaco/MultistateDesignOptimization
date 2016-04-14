#include "stdafx.h"

namespace OPTIMIZER
{
	class SimilarityMeasure
	{
	public:
		SimilarityMeasure();
		SimilarityMeasure(SimilarityMeasure*);
		SimilarityMeasure(mxArray *targetFreqs);

		void setTargetFrequencies(mxArray *targetFreqs);
		virtual SimilarityMeasure *clone();
		virtual double getSimilarity(mxArray *expFrequencies);

		virtual string *toString();

		~SimilarityMeasure();

	protected:
		mxArray *targetFrequencies;
	};
}