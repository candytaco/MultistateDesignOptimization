#include "CuckooSearch.h"

namespace OPTIMIZER
{
	CuckooSearch::CuckooSearch()
	{
		CuckooSearch(0, NULL, NULL, 0, 1, 0.25, false);
	}

	CuckooSearch::CuckooSearch(int nMacrostates, map<int, Model> *models, SimilarityMeasure *similarityMeasure, int populationSize, double scaleParam, double elimination)
	{
		CuckooSearch(nMacrostates, models, similarityMeasure, 32, 1, 0.25, false);
	}
	CuckooSearch::CuckooSearch(int nMacrostates, map<int, Model> *models, SimilarityMeasure *similarityMeasure, int populationSize, double scaleParam, double elimination, bool continuousBoltzmann) : SearchAlgorithm(nMacrostates, models, similarityMeasure, continuousBoltzmann)
	{
		this->scaleParam = scaleParam;
		this->populationSize = populationSize;
		this->elimination = elimination;
		population = new list<Model>();
	}

	void CuckooSearch::initPopulation()
	{
		// TODO: implement this
		bestMatchVal = 0;
		population->clear();
		for (int i = 0; i < populationSize; i++)
		{
			double steepness = randDouble(randGen) * (steepnessRange[1] - steepnessRange[0]) + steepnessRange[0];
			double *weights = new double[nMacrostates];
			// assumes weight ranges are all [0, 1]
			for (int i = 0; i < nMacrostates; i++)
				weights[i] = searchWeights[i] ? randDouble(randGen) : weightMins[i];
			int ensembleSize = ensembleSizes[randGen() % nEnsembleSizes];
			int backrubTemp = backrubTemps[randGen() & nBackrubTemps];
			float boltzmann = continuousBoltzmann ? randDouble(randGen) * (boltzmannTemps[1] - boltzmannTemps[0]) + boltzmannTemps[0] : boltzmannTemps[randGen() % nBoltzmannTemps];

			Model *m;
			if (!continuousBoltzmann)
				m = new Model(*this->getModelByParams(backrubTemp, ensembleSize, boltzmann), ensembleSize, backrubTemp, boltzmann, weights, steepness);
			else
				m = new Model(*this->getModelByParams(backrubTemp, 0, 0), ensembleSize, backrubTemp, boltzmann, weights, steepness);
			m->recovery = similarityMeasure->getSimilarity(m->getFrequencies());
			population->push_back(*m);
		}
		population->sort();
		recordBestParams();
	}

	void CuckooSearch::iterate()
	{
		// TODO: finish implementing this
		initPopulation();
		// TODO: MPI things
		for (int iteration = 0; iteration < maxIterations; iteration++)
		for (int individual = 0; individual < populationSize; individual++)
		{
			// TODO: the meaty things here.
		}
	}

	void CuckooSearch::recordBestParams()
	{
		Model *best = &population->front();
		bestEnsembleSize = best->getEnsembleSize();
		bestBackrubTemp = best->getBackrubTemp();
		bestBoltzmannTemp = best->getBoltzmannTemp();
		bestSteepness = best->getSteepness();
		bestWeights = best->getWeights();
		bestMatchVal = best->recovery;
		bestFrequencies = best->getFrequencies();
	}

	string *CuckooSearch::toString()
	{
		return new string("Cuckoo search");
	}
}