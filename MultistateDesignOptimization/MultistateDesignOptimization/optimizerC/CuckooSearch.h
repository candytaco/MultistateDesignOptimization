#include "stdafx.h"

namespace OPTIMIZER
{
	class CuckooSearch : public SearchAlgorithm
	{
	public:
		CuckooSearch();
		CuckooSearch(int nMacrostates, map<int, Model> *models, SimilarityMeasure *similarityMeasure, int populationSize, double scaleParam, double elimination);
		CuckooSearch(int nMacrostates, map<int, Model> *models, SimilarityMeasure *similarityMeasure, int populationSize, double scaleParam, double elimination, bool continuousBoltzmann);

		void iterate() override;

		string *toString();

	private:
		void initPopulation();
		double nextLevyStep();
		double *nextLevySteps(int n);
		void recordBestParams() override;

		double scaleParam;
		double elimination;
		int populationSize;
		list<Model> *population;
	};
}