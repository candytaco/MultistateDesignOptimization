#pragma once
#include "stdafx.h"
#include "Model.h"
#include "SimilarityMeasure.h"
#include "SearchAlgorithm.h"


namespace OPTIMIZER
{
	class CuckooSearch : public SearchAlgorithm
	{
	public:
		CuckooSearch();
		CuckooSearch(int nMacrostates, map<int, Model> *models, SimilarityMeasure *similarityMeasure, int populationSize, double scaleParam, double elimination);
		CuckooSearch(int nMacrostates, map<int, Model> *models, SimilarityMeasure *similarityMeasure, int populationSize, double scaleParam, double elimination, bool continuousBoltzmann);

		inline void startSearch() override { this->iterate(); };
		void iterate() override;

		string *toString() override;

	private:
		void initilizeMembers();
		void initPopulation();
		double nextLevyStep();
		double *nextLevySteps(int n);
		void recordBestParams() override;
		static bool sortCompModels(const Model &lhs, const Model &rhs);

		double scaleParam;
		double elimination;
		int populationSize;
		vector<Model> *population;
		mt19937 *e;
		boost::math::normal *normal_dist;
		boost::random::uniform_real_distribution<double> *uniform_dist05;
	};
}