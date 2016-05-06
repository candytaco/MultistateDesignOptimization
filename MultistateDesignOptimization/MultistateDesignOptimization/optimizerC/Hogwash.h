#include "stdafx.h"
#include "SearchAlgorithm.h"
#include "CuckooSearch.h"

namespace OPTIMIZER
{	
	/// <summary>
	/// Hog̶w̶a̶s̶hwild SGD adapted to search
	/// </summary>
	/// <seealso cref="CuckooSearch" />
	class Hogwash : public CuckooSearch
	{
		Hogwash();
		Hogwash(int nMacrostates, map<int, Model> *models, SimilarityMeasure *similarityMeasure, int populationSize, double scaleParam, double elimination);
		Hogwash(int nMacrostates, map<int, Model> *models, SimilarityMeasure *similarityMeasure, int populationSize, double scaleParam, double elimination, bool continuousBoltzmann);

		void iterate() override;

		string *toString() override;

	};
}