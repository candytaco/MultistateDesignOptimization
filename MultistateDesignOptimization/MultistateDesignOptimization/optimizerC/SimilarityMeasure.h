#include "stdafx.h"

namespace OPTIMIZER
{
	class SimilarityMeasure
	{
	public:
		SimilarityMeasure();
		SimilarityMeasure(SimilarityMeasure*);
		SimilarityMeasure(mat *targetFreqs);

		void setTargetFrequencies(mat *targetFreqs);
		virtual SimilarityMeasure *clone();	// TODO: not needed because copy constructor
		virtual double getSimilarity(mat *expFrequencies);

		virtual string *toString();

		~SimilarityMeasure();

	protected:
		mat *targetFrequencies;
	};
}