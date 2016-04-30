#include "SearchAlgorithm.h"
#include "Optimizer.h"

namespace OPTIMIZER
{
	SearchAlgorithm::SearchAlgorithm()
	{
		SearchAlgorithm(0, NULL, NULL, false);
	}

	SearchAlgorithm::SearchAlgorithm(int nMacrostates, map<int, Model> *models, SimilarityMeasure *similarityMeasure)
	{
		SearchAlgorithm(nMacrostates, models, similarityMeasure, false);
	}

	SearchAlgorithm::SearchAlgorithm(int nMacrostates, map<int, Model> *models, SimilarityMeasure *similarityMeasure, bool continuousBoltzmann)
	{
		this->models = models;
		this->similarityMeasure = similarityMeasure;
		this->continuousBoltzmann = continuousBoltzmann;
		this->nMacrostates = nMacrostates;

		maxIterations = 0;
		ensembleSizes = NULL;
		backrubTemps = NULL;
		boltzmannTemps = NULL;
		steepnessRange = NULL;
		weightMins = NULL;
		weightMaxs = NULL;

		searchEnsemble =
			searchBackrub =
			searchBoltzmann =
			searchSteepness = true;

		searchWeights = false;

		bestEnsembleSize = -1;
		bestBackrubTemp = -1;
		bestBoltzmannTemp = -1;
		bestSteepness = -1;
		bestWeights = NULL;
		bestMatchVal = -1;
		bestFrequencies = NULL;

		suppressOutputs = true;
	}

	const int SearchAlgorithm::getNumMacrostates()
	{
		return this->nMacrostates;
	}

	void SearchAlgorithm::setSimilarityMeasure(SimilarityMeasure *similarityMeasure)
	{
		if (this->similarityMeasure != NULL)
			delete this->similarityMeasure;
		this->similarityMeasure = similarityMeasure;
	}

	void SearchAlgorithm::setMaxIterations(int iter)
	{
		this->maxIterations = iter;
	}

	void SearchAlgorithm::setSearchparameters(bool ensemble, bool backrub, bool boltzmann, bool steepness, bool *weights)
	{
		searchEnsemble = ensemble;
		searchBackrub = backrub;
		searchBoltzmann = boltzmann;
		if (this->searchWeights)
			delete[] searchWeights;
		searchWeights = weights;
	}

	void SearchAlgorithm::setParameterBounds(int *ensembleSizes, int nEnsembleSizes, double *backrubTemps, int nBackrubTemps, double *boltzmannTemps, int nBoltzmannTemps, double *steepnessRange, double *weightMins, double *weightMaxs)
	{
		if (this->ensembleSizes)
			delete[] this->ensembleSizes;
		this->ensembleSizes = ensembleSizes;
		this->nEnsembleSizes = nEnsembleSizes;
		if (this->backrubTemps)
			delete[] this->backrubTemps;
		this->nBackrubTemps = nBackrubTemps;
		if (this->boltzmannTemps)
			delete[] this->boltzmannTemps;
		this->boltzmannTemps = boltzmannTemps;
		this->nBoltzmannTemps = nBoltzmannTemps;
		if (this->steepnessRange)
			delete[] this->steepnessRange;
		this->steepnessRange = steepnessRange;
		if (this->weightMaxs)
			delete[] this->weightMaxs;
		this->weightMaxs = weightMaxs;
		if (this->weightMins)
			delete[] this->weightMins;
		this->weightMins = weightMins;
	}

	Model *SearchAlgorithm::getModelByParams(double param1, double param2, double param3)
	{
		return &models->at(calcParamsID(param1, param2, param3));
	}

	double *SearchAlgorithm::getBestParams()
	{
		double *out = new double[5 + nMacrostates];
		out[0] = bestEnsembleSize;
		out[1] = bestBackrubTemp;
		out[2] = bestBoltzmannTemp;
		out[3] = bestSteepness;
		for (int i = 4; i < 4 + nMacrostates; i++)
			out[i] = bestWeights[i - 4];
		out[4 + nMacrostates] = bestMatchVal;
		return out;
	}

	double **SearchAlgorithm::getBestFrequencies()
	{
		double **out = new double*[bestFrequencies->n_rows];
		for (int i = 0; i < bestFrequencies->n_rows; i++)
		{
			out[i] = new double[bestFrequencies->n_cols];
			for (int j = 0; j < bestFrequencies->n_cols; j++)
				out[i][j] = bestFrequencies->at(i, j);
		}
		return out;
	}

	void SearchAlgorithm::boundCheckBoltzmann(double *newBoltzmann)
	{
		if (!continuousBoltzmann)
			throw exception::exception("This search not on continuous Boltzmann range");

		if (!searchBoltzmann)
			*newBoltzmann = boltzmannTemps[0];
		else if (*newBoltzmann < -1)
			*newBoltzmann = boltzmannTemps[1];
		else if (*newBoltzmann == 0)
			*newBoltzmann = boltzmannTemps[0];
		else if (*newBoltzmann < boltzmannTemps[0])
			*newBoltzmann = 0;
		else if (*newBoltzmann > boltzmannTemps[1])
			*newBoltzmann = -1;
		return;
	}

	void SearchAlgorithm::boundCheckSteepness(double *newSteepness)
	{
		if (!searchSteepness)
			*newSteepness = steepnessRange[0];
		else if (*newSteepness < steepnessRange[0])
			*newSteepness = steepnessRange[0];
		else if (*newSteepness > steepnessRange[1])
			*newSteepness = steepnessRange[1];
		return;
	}

	void SearchAlgorithm::boundCheckWeights(double *weights)
	{
		// oh Jesus Christ this is ripe for segfaults if the arrays aren't indexed correctly
		for (int i = 0; i < nMacrostates; i++)
		{
			if (!searchWeights[i])
				weights[i] = weightMins[i];
			else if (weights[i] < weightMins[i])
				weights[i] = weightMins[i];
			else if (weights[i] > weightMaxs[i])
				weights[i] = weightMaxs[i];
			else continue;
		}
	}

	SearchAlgorithm::~SearchAlgorithm()
	{
		if (models)
			delete models;
		if (similarityMeasure)
			delete this->similarityMeasure;
		if (ensembleSizes)
			delete[] ensembleSizes;
		if (backrubTemps)
			delete[] backrubTemps;
		if (boltzmannTemps)
			delete[] boltzmannTemps;
		if (steepnessRange)
			delete[] steepnessRange;
		if (weightMins)
			delete[] weightMins;
		if (weightMaxs)
			delete[] weightMaxs;
		if (searchWeights)
			delete[] searchWeights;
		if (bestWeights)
			delete[] bestWeights;
		if (bestFrequencies)
			delete this->bestFrequencies;
	}
}