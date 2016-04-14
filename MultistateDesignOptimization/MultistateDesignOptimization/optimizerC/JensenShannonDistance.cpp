#include "JensenShannonDistance.h"

namespace OPTIMIZER
{
	JensenShannonDistance::JensenShannonDistance() : SimilarityMeasure() {}

	JensenShannonDistance::JensenShannonDistance(mxArray *targetFreqs) : SimilarityMeasure(targetFreqs) {}


}