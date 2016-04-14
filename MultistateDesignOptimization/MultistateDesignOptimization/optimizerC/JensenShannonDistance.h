#include "stdafx.h"

namespace OPTIMIZER
{
	class JensenShannonDistance : public SimilarityMeasure
	{
	public:
		JensenShannonDistance();
		JensenShannonDistance(mxArray *targetFreqs);

		JensenShannonDistance *clone() override;

		double getSimilarity(mxArray *expFrequencies) override;

		~JensenShannonDistance();

		string *toString() override;

	};
}