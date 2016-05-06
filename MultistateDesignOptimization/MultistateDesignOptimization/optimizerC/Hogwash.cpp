#include "stdafx.h"
#include "Hogwash.h"

namespace OPTIMIZER
{
	Hogwash::Hogwash() : CuckooSearch()
	{}

	Hogwash::Hogwash(int nMacrostates, map<int, Model> *models, SimilarityMeasure *similarityMeasure, int populationSize, double scaleParam, double elimination)
		: CuckooSearch(nMacrostates, models, similarityMeasure, populationSize, scaleParam, elimination)
	{}

	Hogwash::Hogwash(int nMacrostates, map<int, Model> *models, SimilarityMeasure *similarityMeasure, int populationSize, double scaleParam, double elimination, bool continuousBoltzmann)
		: CuckooSearch(nMacrostates, models, similarityMeasure, populationSize, scaleParam, elimination, continuousBoltzmann)
	{}

	void Hogwash::iterate()
	{
		// TODO: finish implementing this
		// TODO: a bunch of these things can be set at the beginning of the code.
		initPopulation();
		clock_t start = clock();
		omp_lock_t lock;
		omp_init_lock(&lock);

		// TODO: MPI things
		// TODO: make all of the variables private?
#pragma omp parallel nowait
		for (int iteration = 0; iteration < maxIterations; iteration++)
		{
			omp_set_num_threads(12);

			//list<Model>::iterator it = population.begin();
			int individual;
#pragma omp parallel
			{
				int numthreads = omp_get_num_threads();
#pragma omp master
				{printf("Iteration %d! \n", iteration);
				printf("number of openmp threads = %d\n", numthreads); }

				Model *newModel = NULL, *temp = NULL;
				//bool createModel = true;
#pragma omp for private(newModel) private(temp) nowait	
				for (individual = 0; individual < populationSize; individual++)
				{
					Model it = *population.at(individual); // because now population is a vector.
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

					if (!continuousBoltzmann)
					{
						newBoltzmannTemp = searchBoltzmann ? boltzmannTemps[randGen() % nBoltzmannTemps] : boltzmannTemps[0];
						//if (createModel)
						{
							newModel = new Model(*this->getModelByParams(newBackrubTemp, newEnsembleSize, newBoltzmannTemp), newEnsembleSize, newBackrubTemp, newBoltzmannTemp, newWeights, newSteep);
							//	createModel = false;
						}
						//else
						//	newModel->setParameters(newBoltzmannTemp, newWeights, newSteep, newEnsembleSize);
					}
					else
					{
						newBoltzmannTemp = nextLevyStep() + it.getBoltzmannTemp(); //TODO: check and figure out how we want to generate this.
						// python code:
						//	boltzmannStep = multiplier * (self.population[randParent1]->getBoltzmannTemp() - self.population[randParent2]->getBoltzmannTemp());
						boundCheckBoltzmann(&newBoltzmannTemp);
						//if (createModel)
						//{
						newModel = new Model(*this->getModelByParams(newBackrubTemp, 0, 0), newEnsembleSize, newBackrubTemp, newBoltzmannTemp, newWeights, newSteep);
						//	createModel = false;
						//}
						//else
						//	newModel->setParameters(newBoltzmannTemp, newWeights, newSteep, newEnsembleSize);
					}
					newModel->recovery = similarityMeasure->getSimilarity(newModel->getFrequencies());

					// TODO: can we find a way to do this without using *newModel? Or we just need to clear newModel after this happens I suppose
					if (*newModel > it)
					{
						omp_set_lock(&lock);
						if (newBackrubTemp == population.at(individual)->getBackrubTemp())
						{
							population.at(individual)->setParameters(newBoltzmannTemp, newWeights, newSteep, newEnsembleSize);
							delete newModel;
						}
						else
						{
							temp = population.at(individual);
							population.at(individual) = newModel;
							delete temp;
							newModel = NULL;
							//	createModel = true;
						}
						omp_unset_lock(&lock);
					}
					else
						delete newModel;
					//population.at(individual) = *newModel;

					//TODO: this is causing a segfault.
					//newModel->clearModel();
					//delete[] newModel;

					if (individual != 0)
					{ // we don't want to randomly replace the best nest. this gives us some bias towards the best nest... (because it will always be there to get randomly selected).
						if (randDouble(randGen) < elimination)
						{
							int randParent1 = randGen() % populationSize;
							int randParent2 = randGen() % populationSize;
							double multiplier = randDouble(randGen);
							double steepnessStep = multiplier * (population.at(randParent1)->getSteepness() - population.at(randParent2)->getSteepness());
							double *parent1Weights = population.at(randParent1)->getWeights();
							double *parent2Weights = population.at(randParent2)->getWeights();
							newWeights = population.at(individual)->getWeights();
							for (int i = 0; i < nMacrostates; i++)
								newWeights[i] += multiplier * (parent1Weights[i] - parent2Weights[i]);
							newSteep = population.at(individual)->getSteepness() + steepnessStep;
							boundCheckSteepness(&newSteep);
							boundCheckWeights(newWeights);

							newEnsembleSize = searchEnsemble ? ensembleSizes[randGen() % nEnsembleSizes] : ensembleSizes[0];
							newBackrubTemp = searchBackrub ? backrubTemps[randGen() % nBackrubTemps] : backrubTemps[0];

							// get rid of the new model here because we are just updating the old model.
							double newBoltzmannTemp;
							if (!continuousBoltzmann)
							{
								newBoltzmannTemp = searchBoltzmann ? boltzmannTemps[randGen() % nBoltzmannTemps] : boltzmannTemps[0];
								newModel = new Model(*this->getModelByParams(newBackrubTemp, newEnsembleSize, newBoltzmannTemp), newEnsembleSize, newBackrubTemp, newBackrubTemp, newWeights, newSteep);
							}
							else
							{
								newBoltzmannTemp = multiplier * (population.at(randParent1)->getBoltzmannTemp() - population.at(randParent2)->getBoltzmannTemp()); //TODO: check and figure out how we want to generate this.
								// python code:
								//	boltzmannStep = multiplier * (self.population[randParent1]->getBoltzmannTemp() - self.population[randParent2]->getBoltzmannTemp());
								newBoltzmannTemp += population.at(individual)->getBoltzmannTemp();
								boundCheckBoltzmann(&newBoltzmannTemp);
								newModel = new Model(*this->getModelByParams(newBackrubTemp, 0, 0), newEnsembleSize, newBackrubTemp, newBackrubTemp, newWeights, newSteep);
							}

							newModel->recovery = similarityMeasure->getSimilarity(newModel->getFrequencies());

							// this is currently a memory leak here - the old model is not deleted
							//population.at(individual) = *newModel; // this is where we will need to sync across the processes!
							omp_set_lock(&lock);
							if (newBackrubTemp == population.at(individual)->getBackrubTemp())
							{
								population.at(individual)->setParameters(newBoltzmannTemp, newWeights, newSteep, newEnsembleSize);
								delete newModel;
							}
							else
							{
								temp = population.at(individual);
								population.at(individual) = newModel;
								delete temp;
								newModel = NULL;
								//	createModel = true;
							}
							omp_unset_lock(&lock);
						}
					}

				}
			}

			//TODO: end of this... now quite sure what is supposed to happen
			// should only happen on one process? so we need to bring everything back together!
			std::sort(population.begin(), population.end(), &CuckooSearch::sortCompModels);
			recordBestParams();
			//elimination of the worst individuals
			/*for (int i = 0; i < int(populationSize * elimination); i++)
			population.pop_back();*/
		}
		clock_t end = clock();
		elapsedTime = double(end - start) / CLOCKS_PER_SEC;
	}

}