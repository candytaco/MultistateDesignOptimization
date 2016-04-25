#include "stdafx.h"

namespace OPTIMIZER
{
	class JensenShannonDistance : public SimilarityMeasure
	{
	public:
		JensenShannonDistance();
		JensenShannonDistance(const SimilarityMeasure &other);
		JensenShannonDistance(mat *targetFreqs);

		double getSimilarity(mat *expFrequencies) override;

		~JensenShannonDistance();

		string *toString() override;

	};
}