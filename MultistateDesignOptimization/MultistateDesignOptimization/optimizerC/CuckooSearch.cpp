#include "CuckooSearch.h"
// these should maybe go somewhere else.
#include <random>

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
		e = new mt19937(time(NULL));
		normal_dist = new boost::math::normal(0.0, 1.0); // make the normal distribution
		uniform_dist05 = new boost::random::uniform_real_distribution<double>(0.5, 1); // make the uniform distribution
	}

	void CuckooSearch::initPopulation()
	{
		// TODO: implement this
		bestMatchVal = 0;
		population->clear();
		for (int i = 0; i < populationSize; i++)
		{
            double steepness = searchSteepness ? randDouble(randGen) * (steepnessRange[1] - steepnessRange[0]) + steepnessRange[0] : steepnessRange[0];
			double *weights = new double[nMacrostates];
            //TODO: add the range from the pyton code is going ton with the range?
			// assumes weight ranges are all [0, 1]
			for (int i = 0; i < nMacrostates; i++)
				weights[i] = searchWeights[i] ? randDouble(randGen) : weightMins[i]; // boolean ? <then this> : <else>
            int ensembleSize = searchEnsemble ? ensembleSizes[randGen() % nEnsembleSizes] : ensembleSizes[0];
            int backrubTemp = searchBackrub ? backrubTemps[randGen() % nBackrubTemps] : backrubTemps[0];
            double boltzmannTemp = searchBoltzmann ? boltzmannTemps[randGen() % nBoltzmannTemps] : boltzmannTemps[0];
            
            // should this next line be here?
            // double boltzmann = continuousBoltzmann ? randDouble(randGen) * (boltzmannTemps[1] - boltzmannTemps[0]) + boltzmannTemps[0] : boltzmannTemps[randGen() % nBoltzmannTemps];
            
			Model *m;
			if (!continuousBoltzmann)
				m = new Model(*this->getModelByParams(backrubTemp, ensembleSize, boltzmannTemp), ensembleSize, backrubTemp, boltzmannTemp, weights, steepness);
			else
				m = new Model(*this->getModelByParams(backrubTemp, 0, 0), ensembleSize, backrubTemp, boltzmannTemp, weights, steepness);
            //m->macrostatesUsed = searchWeights; // not sure if this line is necessary.
			m->recovery = similarityMeasure->getSimilarity(m->getFrequencies());
			population->push_back(*m);
		}
		population->sort();
		recordBestParams();
	}

	void CuckooSearch::iterate()
	{
		// TODO: finish implementing this
        // TODO: a bunch of these things can be set at the beginning of the code.
		initPopulation();
		// TODO: MPI things
		for (int iteration = 0; iteration < maxIterations; iteration++)
		{
			list<Model>::iterator it = population->begin();
			for (int individual = 0; individual < populationSize; individual++)
			{
				// TODO: the meaty things here.
                // TODO: can we do the search on a coarse grid first and then the smoother grid?
				double newSteep = nextLevyStep() + it->getSteepness();
				boundCheckSteepness(&newSteep);
                
				double *newWeights = nextLevySteps(nMacrostates);
				double *oldWeights = it->getWeights();
				for (int i = 0; i < nMacrostates; i++)
					newWeights[i] += oldWeights[i]; // is this correct? should this be 5 entries?

                boundCheckWeights(newWeights);
                
                double newEnsembleSize = searchEnsemble ? ensembleSizes[randGen() % nEnsembleSizes] : ensembleSizes[0];
				double newBackrubTemp = searchBackrub ? backrubTemps[randGen() % nBackrubTemps] : backrubTemps[0];
                
				Model *newModel;
                if (!continuousBoltzmann) {
                    double newBoltzmanTemp = searchBoltzmann ? boltzmannTemps[randGen() % nBoltzmannTemps] : boltzmannTemps[0];
                    newModel = new Model(*this->getModelByParams(newBackrubTemp, newEnsembleSize, newBoltzmanTemp), newEnsembleSize, newBackrubTemp, newBackrubTemp, newWeights, newSteep);
                }
                else {
                    double newBoltzmannTemp = nextLevyStep() + it->getBoltzmannTemp(); //TODO: check and figure out how we want to generate this.
                    // python code:
                    //	boltzmannStep = multiplier * (self.population[randParent1].getBoltzmannTemp() - self.population[randParent2].getBoltzmannTemp());
                    boundCheckBoltzmann(&newBoltzmannTemp);
                    newModel = new Model(*this->getModelByParams(newBackrubTemp, newEnsembleSize, newBoltzmannTemp), newEnsembleSize, newBackrubTemp, newBackrubTemp, newWeights, newSteep);
                }
                newModel->recovery = similarityMeasure->getSimilarity(newModel->getFrequencies);
                
                if (newModel->recovery > it->recovery)
                    
			}

			//elimination of the worst individuals
			for (int i = 0; i < int(populationSize * elimination); i++)
				population->pop_back();
		}
	}
    
    double CuckooSearch::nextLevyStep()
    {
        /* Generates a random number from this model's Levy distribution.
         The magnitude is drawn from a Levy distribution f(x; 0, scaleParam)
         While the direction is uniform random
         
         @param void
         @return float
         
         Use the BOOST library quantile function.
         */
        
        // something like the following lines needs to be called at the beginning of the code, but I'm not sure where that is. we do not need to redefine the random number generation and the distributions EVERY TIME we run this function (that would be silly).
           
            // draw from 0.5-1 because we are taking 1-r1_old/2, so we might as well save that computationtime.
        double r1 = uniform_dist05->operator()(*e);
        double ppf = quantile(*normal_dist, r1);
        // this is slightly different than the python code, but i think the python code should have been multiplying instead of dividing?double check! (CHECKED 4/30. This is correct).
        double randLevy = scaleParam * pow(ppf,-2);
        double signof = uniform_dist05->operator()(*e) - 0.75;
        float randLevySigned = copysign(randLevy,signof) * 0.01; // Cuckoo search authors says to use 1/100 of the scale length
            
            //     return randLevySigned;
            //printf("%f",randLevySigned);
        return randLevySigned;
    }
    
    double *CuckooSearch::nextLevySteps(int steps) {
        /*
         Generates an array of Levy steps
         
         @param steps		size of array
         @return float[]		an array of independent Levy steps
         */
		double *out = new double[steps];
        for (int i = 0; i< steps; i++)
            out[i] = nextLevyStep();
        return out;
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
