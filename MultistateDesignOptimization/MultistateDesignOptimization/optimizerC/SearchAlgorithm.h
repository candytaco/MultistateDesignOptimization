#include "stdafx.h"

namespace OPTIMIZER
{
	class SearchAlgorithm
	{
	public:
		SearchAlgorithm();
		SearchAlgorithm(map<int, Model> models, SimilarityMeasure similarityMeasure);
		SearchAlgorithm(map<int, Model> models, SimilarityMeasure similarityMeasure, bool continuousBoltzmann);

		void setSimilarityMeasure(SimilarityMeasure similarityMeasure);
		void setMaxIterations(int iter);
		void setContinuousBoltzmann(bool yes);
		bool getContinuousBoltzmann();
		void setParameterBounds(int *ensembleSizes, double *backrubTemps, double *boltzmannTemps, double *steepnessRange, double *weightMins, double *weightMaxs);
		void setSearchparameters(bool ensemble, bool backrub, bool boltzmann, bool steepness, bool *weights);

		virtual void startSearch();
		virtual void iterate();

		Model *getModelByParams(double, double, double);

		double *getBestParams();
		float **getBestFrequencies();

		bool suppressOutputs;

		~SearchAlgorithm();

	protected:

		void boundCheckBoltzmann(double *newBoltzmann);
		void boundCheckSteepness(double *newSteep);
		void boundCheckWeights(double *weights);
		virtual void recordBestParams();

		map<int, Model> models;
		SimilarityMeasure *similarityMeasure;
		int maxIterations;

		int *ensembleSizes;
		double *backrubTemps;
		double *boltzmannTemps;
		bool continuousBoltzmann;
		double *steepnessRange;
		double *weightMins;
		double *weightMaxs;

		bool searchEnsemble;
		bool searchBackrub;
		bool searchBoltzmann;
		bool searchSteepness;
		bool *searchWeights;

		int bestEnsembleSize;
		double bestBackrubTemp;
		double bestBoltzmannTemp;
		double bestSteepness;
		double *bestWeights;
		double bestMatchVal;
		mat *bestFrequencies;
	};
}