#include "stdafx.h"
#include "Model.h"
#include "SimilarityMeasure.h"
#include "CuckooSearch.h"

namespace OPTIMIZER
{
	CuckooSearch::CuckooSearch() : SearchAlgorithm(0, NULL, NULL, false),
		scaleParam(1),
		populationSize(0),
		elimination(0.25)
	{
		this->initilizeMembers();
	}

	CuckooSearch::CuckooSearch(int nMacrostates, map<int, Model> *models, SimilarityMeasure *similarityMeasure, int populationSize, double scaleParam, double elimination) : SearchAlgorithm(nMacrostates, models, similarityMeasure, false),
		scaleParam(scaleParam),
		populationSize(populationSize),
		elimination(elimination)
	{
		this->initilizeMembers();
	}
	CuckooSearch::CuckooSearch(int nMacrostates, map<int, Model> *models, SimilarityMeasure *similarityMeasure, int populationSize, double scaleParam, double elimination, bool continuousBoltzmann) : SearchAlgorithm(nMacrostates, models, similarityMeasure, continuousBoltzmann),
		scaleParam(scaleParam),
		populationSize(populationSize),
		elimination(elimination)
	{
		this->initilizeMembers();
	}

	void CuckooSearch::initilizeMembers()
	{
		population = vector<Model>();
        population.resize(populationSize);
		e = new mt19937(time(NULL));
		normal_dist = new boost::math::normal(0.0, 1.0); // make the normal distribution
		uniform_dist05 = new boost::random::uniform_real_distribution<double>(0.5, 1); // make the uniform distribution
		searchWeights = new bool[nMacrostates];
		for (int i = 0; i < nMacrostates; i++)
			searchWeights[i] = true;
	}

	bool CuckooSearch::sortCompModels(const Model&lhs, const Model &rhs)
	{
		return lhs > rhs;
	}

	void CuckooSearch::initPopulation()
	{
		// TODO: implement this
		bestMatchVal = 0;
        // TODO: commented the next line out when trying to debug the model deconstructor.
		//population.clear();
		for (int i = 0; i < populationSize; i++)
		{
            double steepness = searchSteepness ? randDouble(randGen) * (steepnessRange[1] - steepnessRange[0]) + steepnessRange[0] : steepnessRange[0];
			double *weights = new double[nMacrostates];
            //TODO: add the range from the pyton code is going ton with the range?
			// assumes weight ranges are all [0, 1]
			for (int j = 0; j < nMacrostates; j++)
				weights[j] = searchWeights[j] ? randDouble(randGen) : weightMins[j]; // boolean ? <then this> : <else>
            int ensembleSize = searchEnsemble ? ensembleSizes[randGen() % nEnsembleSizes] : ensembleSizes[0];
            double backrubTemp = searchBackrub ? backrubTemps[randGen() % nBackrubTemps] : backrubTemps[0];
            
            // should this next line be here?
            // double boltzmann = continuousBoltzmann ? randDouble(randGen) * (boltzmannTemps[1] - boltzmannTemps[0]) + boltzmannTemps[0] : boltzmannTemps[randGen() % nBoltzmannTemps];
            
			Model *m;
			if (!continuousBoltzmann)
			{
				double boltzmannTemp = searchBoltzmann ? boltzmannTemps[randGen() % nBoltzmannTemps] : boltzmannTemps[0];
				m = new Model(*this->getModelByParams(backrubTemp, ensembleSize, boltzmannTemp), ensembleSize, backrubTemp, boltzmannTemp, weights, steepness);
			}
			else
			{
				double boltzmannTemp = searchBoltzmann ? randDouble(randGen) * boltzmannTemps[1] : boltzmannTemps[0];
				m = new Model(*this->getModelByParams(backrubTemp, 0, 0), ensembleSize, backrubTemp, boltzmannTemp, weights, steepness);
			}
            //m->macrostatesUsed = searchWeights; // not sure if this line is necessary.
			m->recovery = similarityMeasure->getSimilarity(m->getFrequencies());
            printf("Here I am!\n");
            // this is a memory leak because we aren't deconstructing the models properly... but not sure what to do about that?
            population.at(i) = *m;
            //population.push_back(*m);
            //population[i] = *m;
            //population.push_back(*m);
		}
		sort(population.begin(), population.end(), &CuckooSearch::sortCompModels); // needed a 2-arg comparator
		//population.sort();
		recordBestParams();
	}

	void CuckooSearch::iterate()
	{
		// TODO: finish implementing this
        // TODO: a bunch of these things can be set at the beginning of the code.
		initPopulation();
		clock_t start = clock();
        omp_lock_t lock;
        omp_init_lock(&lock);
		// TODO: MPI things
        // TODO: make all of the variables private?
		for (int iteration = 0; iteration < maxIterations; iteration++)
		{
            printf("Iteration %d! \n",iteration);
            
			//list<Model>::iterator it = population.begin();
            int individual;
            #pragma omp parallel for
			for (individual = 0; individual < populationSize; individual++)
			{
                Model *newModel, oldModel;
                int numthreads = omp_get_num_threads();
                printf("number of openmp threads = %d\n",numthreads);
				Model it = population.at(individual); // because now population is a vector.
				// TODO: the meaty things here.
                // TODO: can we do the search on a coarse grid first and then the smoother grid?
                // TODO: all of these things should not be redefined within each iteration. pull them out because this is a memory leak?
				double newSteep = nextLevyStep() + it.getSteepness();
				boundCheckSteepness(&newSteep);
                
				double *newWeights = nextLevySteps(nMacrostates);
				double *oldWeights = it.getWeights();
				for (int i = 0; i < nMacrostates; i++)
					newWeights[i] += oldWeights[i]; // is this correct? should this be 5 entries?

                boundCheckWeights(newWeights);
                
                double newEnsembleSize = searchEnsemble ? ensembleSizes[randGen() % nEnsembleSizes] : ensembleSizes[0];
				double newBackrubTemp = searchBackrub ? backrubTemps[randGen() % nBackrubTemps] : backrubTemps[0];
                double newBoltzmannTemp;
                
                // TODO: in the next few lines we generate a new model but we don't clear it. this is also going to cause a memory leak. Do we have a temp model?
                
                if (!continuousBoltzmann) {
                    newBoltzmannTemp = searchBoltzmann ? boltzmannTemps[randGen() % nBoltzmannTemps] : boltzmannTemps[0];
                    newModel = new Model(*this->getModelByParams(newBackrubTemp, newEnsembleSize, newBoltzmannTemp), newEnsembleSize, newBackrubTemp, newBoltzmannTemp, newWeights, newSteep);
                }
                else {
                    newBoltzmannTemp = nextLevyStep() + it.getBoltzmannTemp(); //TODO: check and figure out how we want to generate this.
                    // python code:
                    //	boltzmannStep = multiplier * (self.population[randParent1].getBoltzmannTemp() - self.population[randParent2].getBoltzmannTemp());
                    boundCheckBoltzmann(&newBoltzmannTemp);
					newModel = new Model(*this->getModelByParams(newBackrubTemp, 0, 0), newEnsembleSize, newBackrubTemp, newBoltzmannTemp, newWeights, newSteep);
                }
                newModel->recovery = similarityMeasure->getSimilarity(newModel->getFrequencies());
                
                // TODO: can we find a way to do this without using *newModel? Or we just need to clear newModel after this happens I suppose
                if (*newModel > it)
                {
                    omp_set_lock(&lock);
                    population.at(individual).setParameters(newBackrubTemp,newBoltzmannTemp,newWeights,newSteep,newEnsembleSize);
                    omp_unset_lock(&lock);
                }
                    //population.at(individual) = *newModel;
                    
                //TODO: this is causing a segfault.
                //newModel->clearModel();
                //delete[] newModel;
                
                if (individual!=0) { // we don't want to randomly replace the best nest. this gives us some bias towards the best nest... (because it will always be there to get randomly selected).
                if (randDouble(randGen) < elimination) {
                    int randParent1 = randGen() % populationSize;
                    int randParent2 = randGen() % populationSize;
                    double multiplier = randDouble(randGen);
                    double steepnessStep = multiplier * (population.at(randParent1).getSteepness() - population.at(randParent2).getSteepness());
					double *parent1Weights = population.at(randParent1).getWeights();
					double *parent2Weights = population.at(randParent2).getWeights();
					newWeights = population.at(individual).getWeights();
					for (int i = 0; i < nMacrostates; i++)
                        newWeights[i] += multiplier * (parent1Weights[i]-parent2Weights[i]);
                    newSteep = population.at(individual).getSteepness() + steepnessStep;
                    boundCheckSteepness(&newSteep);
                    boundCheckWeights(newWeights);
                    
                    newEnsembleSize = searchEnsemble ? ensembleSizes[randGen() % nEnsembleSizes] : ensembleSizes[0];
                    newBackrubTemp = searchBackrub ? backrubTemps[randGen() % nBackrubTemps] : backrubTemps[0];
                    
                    // get rid of the new model here because we are just updating the old model.
					double newBoltzmannTemp;
                    if (!continuousBoltzmann) {
                        newBoltzmannTemp = searchBoltzmann ? boltzmannTemps[randGen() % nBoltzmannTemps] : boltzmannTemps[0];
                        //newModel = new Model(*this->getModelByParams(newBackrubTemp, newEnsembleSize, newBoltzmannTemp), newEnsembleSize, newBackrubTemp, newBackrubTemp, newWeights, newSteep);
                    }
                    else {
                        newBoltzmannTemp = multiplier * (population.at(randParent1).getBoltzmannTemp() - population.at(randParent2).getBoltzmannTemp()); //TODO: check and figure out how we want to generate this.
                        // python code:
                        //	boltzmannStep = multiplier * (self.population[randParent1].getBoltzmannTemp() - self.population[randParent2].getBoltzmannTemp());
                        newBoltzmannTemp += population.at(individual).getBoltzmannTemp();
                        boundCheckBoltzmann(&newBoltzmannTemp);
                        //newModel = new Model(*this->getModelByParams(newBackrubTemp, 0, 0), newEnsembleSize, newBackrubTemp, newBackrubTemp, newWeights, newSteep);
                    }
                    
                    newModel->recovery = similarityMeasure->getSimilarity(newModel->getFrequencies());
                    
					// this is currently a memory leak here - the old model is not deleted
					//population.at(individual) = *newModel; // this is where we will need to sync across the processes!
                    omp_set_lock(&lock);
                    population.at(individual).setParameters(newBackrubTemp,newBoltzmannTemp,newWeights,newSteep,newEnsembleSize);
                    omp_unset_lock(&lock);
                }
                }
                
			}
            
            //TODO: end of this... now quite sure what is supposed to happen
            // should only happen on one process? so we need to bring everything back together!
			sort(population.begin(), population.end(), &CuckooSearch::sortCompModels);
			recordBestParams();
			//elimination of the worst individuals
			/*for (int i = 0; i < int(populationSize * elimination); i++)
				population.pop_back();*/
		}
		clock_t end = clock();
		elapsedTime = double(end - start) / CLOCKS_PER_SEC;
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
		Model *best = &population.front();
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
		char out[200];
#ifndef _WIN32
		sprintf(out, "Cuckoo search running for %d generations", maxIterations);
#else
		sprintf_s(out, 200, "Cuckoo search running for %d generations", maxIterations);
#endif
		return new string(out);
	}
}
