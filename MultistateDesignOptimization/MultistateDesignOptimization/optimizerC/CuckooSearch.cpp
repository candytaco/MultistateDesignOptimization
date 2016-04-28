#include "CuckooSearch.h"
// these should maybe go somewhere else.
#include <random>
#include <boost/math/distributions/normal.hpp>
#include <math>

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
		{
			list<Model>::iterator it = population->begin();
			for (int individual = 0; individual < populationSize; individual++)
			{
				// TODO: the meaty things here.
				double newSteep = nextLevyStep() + it->getSteepness();
				boundCheckSteepness(&newSteep);
				
			}

			//elimination of the worst individuals
			for (int i = 0; i < int(populationSize * elimination); i++)
				population->pop_back();
		}
	}
    
    void CukooSearch::nextLevyStep()
    {
        /* Generates a random number from this model's Levy distribution.
        The magnitude is drawn from a Levy distribution f(x; 0, scaleParam)
        While the direction is uniform random
        
        @param void
        @return float 
         
         Use the BOOST library quantile function.
         */
        
        // something like the following lines needs to be called at the beginning of the code, but I'm not sure where that is. we do not need to redefine the random number generation and the distributions EVERY TIME we run this function (that would be silly).
        default_random_engine e(time(NULL));
        boost::math::normal normal_dist(0.0, 1.0); // make the normal distribution
        boost::math::uniform_real_distribution uniform_dist05(0.5,1); // make the uniform distribution

        // draw from 0.5-1 because we are taking 1-r1_old/2, so we might as well save that computationtime.
        double r1 = uniform_dist05(e)
        double ppf = quantile(normal_dist,r1)
        // this is slightly different than the python code, but i think the python code should have been multiplying instead of dividing? double check!
        double randLevy = scaleParam * pow(ppf,-2)
        double signof = d(e)-0.75
        double randLevySigned = copysign(randLevy,signof) * 0.01 // Cuckoo search authors says to use 1/100 of the scale length

        return randLevySigned
        
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