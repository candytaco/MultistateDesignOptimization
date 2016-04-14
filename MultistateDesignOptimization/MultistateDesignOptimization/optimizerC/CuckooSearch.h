#include "stdafx.h"

namespace OPTIMIZER
{
	class CuokooSearch : SearchAlgorithm
	{
	public:
		CuokooSearch();
		CuokooSearch(map<int, Model> models, SimilarityMeasure similarityMeasure, int populationSize, double scaleParam, double elimination);
		CuokooSearch(map<int, Model> models, SimilarityMeasure similarityMeasure, int populationSize, double scaleParam, double elimination, bool continuousBoltzmann);

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
		list<Model> population;
	};
}