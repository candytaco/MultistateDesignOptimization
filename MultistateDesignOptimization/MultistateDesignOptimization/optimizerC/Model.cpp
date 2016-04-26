#include "Model.h"
#include <cstdlib>

namespace OPTIMIZER
{
	Model::Model()
	{
		Model(0, 0, 0, 0, NULL, 0, 0, 0, false, false);
	}

	Model::Model(int nMacrostates, int ensembleSize, double backrubTemp, double boltzmannTemp, int nPositions, int positionOffset)
	{
		Model(nMacrostates, ensembleSize, backrubTemp, boltzmannTemp, NULL, 0, nPositions, positionOffset, false, false);
	}

	Model::Model(int nMacrostates, int ensembleSize, double backrubTemp, double boltzmannTemp, int nPositions, int positionOffset, bool useMicrostates)
	{
		Model(nMacrostates, ensembleSize, backrubTemp, boltzmannTemp, NULL, 0, nPositions, positionOffset, useMicrostates, false);
	}

	Model::Model(int nMacrostates, int ensembleSize, double backrubTemp, double boltzmannTemp, int nPositions, int positionOffset, bool useMicrostates, bool useAltAvgMethod)
	{
		Model(nMacrostates, ensembleSize, backrubTemp, boltzmannTemp, NULL, 0, nPositions, positionOffset, useMicrostates, useAltAvgMethod);
	}

	Model::Model(int nMacrostates, int ensembleSize, double backrubTemp, double boltzmannTemp, double *weights, double steepness, int nPositions, int positionOffset)
	{
		Model(nMacrostates, ensembleSize, backrubTemp, boltzmannTemp, weights, steepness, nPositions, positionOffset, false, false);
	}

	Model::Model(int nMacrostates, int ensembleSize, double backrubTemp, double boltzmannTemp, double *weights, double steepness, int nPositions, int positionOffset, bool useMicrostates)
	{
		Model(nMacrostates, ensembleSize, backrubTemp, boltzmannTemp, weights, steepness, nPositions, positionOffset, useMicrostates, false);
	}

	Model::Model(int nMacrostates, int ensembleSize, double backrubTemp, double boltzmannTemp, double *weights, double steepness, int nPositions, int positionOffset, bool useMicrostates, bool useAltAvgMethod)
	{
		Model(nMacrostates, ensembleSize, backrubTemp, boltzmannTemp, weights, steepness, nPositions, positionOffset, useMicrostates, useAltAvgMethod);
	}

	Model::Model(int nMacrostates, int ensembleSize, double backrubTemp, double boltzmannTemp, double *weights, double steepness, int nPositions, int positionOffset, bool useMicrostate, bool useAltaverageMethod)
	{
		this->recovery = 0;
		this->nMacrostates = nMacrostates;
		this->ensembleSize = ensembleSize;
		this->nPositions = nPositions;
		this->positionOffset = positionOffset;
		this->backrubTemp = backrubTemp;
		this->boltzmannTemp = boltzmannTemp;
		this->weights = weights;
		this->steepness = steepness;
		this->useAltAverageingMethod = useAltaverageMethod;
		this->useMicrostateData = useMicrostate;
		this->fitnessCalculated = false;
		this->fitnesses = new mat(nPositions, 20);
		this->frequencies = new mat(nPositions, 20);
		this->macrostateResidueEnergies = new cube(nPositions, 20, nMACROSTATES);
		if (useMicrostate)
		{
			areMicrostatesPicked = false;
			microstateResidueEnergies = new vector<cube>();
			for (int i = 0; i < nPositions; i++)
				microstateResidueEnergies->push_back(cube(20, nMACROSTATES, 1024));
			microstateCounts = new int*[nPositions];
			for (int i = 0; i < nPositions; i++)
				microstateCounts[i] = new int[nMacrostates];
			microstatesUsed = new int**[nPositions];
			for (int i = 0; i < nPositions; i++)
			{
				microstatesUsed[i] = new int*[nMacrostates];
				for (int j = 0; j < nMacrostates; j++)
					microstatesUsed[i][j] = new int[ensembleSize];
			}
			
		}
	}

	Model::Model(const Model &existing, int ensembleSize, double backrubTemp, double boltzmannTemp, double *weights, double steepness)
	{
		this->recovery = 0;
		this->nMacrostates = existing.nMacrostates;
		this->ensembleSize = ensembleSize;
		this->backrubTemp = backrubTemp;
		this->boltzmannTemp = boltzmannTemp;
		this->weights = weights;
		this->steepness = steepness;
		this->fitnesses = new mat(nPositions, 20);
		this->frequencies = new mat(nPositions, 20);
		this->useMicrostateData = existing.useMicrostateData;
		this->areMicrostatesPicked = false;
		this->useAltAverageingMethod = existing.useAltAverageingMethod;

		if (!useMicrostateData && ensembleSize == existing.ensembleSize)
			this->macrostateResidueEnergies = new cube(*(existing.macrostateResidueEnergies));
		else if (useMicrostateData && ensembleSize == existing.ensembleSize)
		{
			this->areMicrostatesPicked = true;
			this->selectedMicrostateEnergies = new vector<cube>(*existing.selectedMicrostateEnergies);
		}
		else
		{
			this->areMicrostatesPicked = false;
			this->microstateResidueEnergies = new vector<cube>(*existing.microstateResidueEnergies);
			this->microstateCounts = new int*[nPositions];
			for (int i = 0; i < nPositions; i++)
			{
				microstateCounts[i] = new int[nMacrostates];
				for (int j = 0; j < nMacrostates; j++)
					microstateCounts[i][j] = existing.microstateCounts[i][j];
			}
		}

		// TODO: shallow or deep copy?
		if (existing.microstatesUsed != NULL)
		{
			this->microstatesUsed = new int**[nPositions];
			for (int i = 0; i < nPositions; i++)
			{
				microstatesUsed[i] = new int*[nMacrostates];
				for (int j = 1; j < nMacrostates; j++)
				{
					microstatesUsed[i][j] = new int[ensembleSize];
					for (int k = 1; k < ensembleSize; k++)
						microstatesUsed[i][j][k] = existing.microstatesUsed[i][j][k];
				}
			}
		}

		else this->microstatesUsed = NULL;
	}

	void Model::addMacrostateData(MACROSTATES macrostate, int position, float *energies)
	{
		position -= this->positionOffset;
		for (int i = 0; i < 20; i++)
			(*macrostateResidueEnergies)(position, macrostate, i) = energies[i];
	}

	void Model::addMicrostateData(MACROSTATES macrostate, int position, float *energies)
	{
		position -= this->positionOffset;
		int microstateIndex = this->microstateCounts[position][macrostate];
		for (int i = 0; i < 20; i++)
			microstateResidueEnergies->at(position).at(i, macrostate, microstateIndex) = energies[i];
		microstateCounts[position][macrostate]++;
	}

	double *Model::getWeights()
	{
		double *w = new double[nMacrostates];
		for (int i = 0; i < nMacrostates; i++)
			w[i] = this->weights[i];
		return w;
	}

	mat *Model::getFrequencies()
	{
		if (this->ensembleSize == 0) return NULL;

		if (!isFrequenciesCalculated)
			calcFrequencies();
		
		return new mat(*frequencies);
	}

	void Model::calcFitness()
	{
		// TODO: check correctness
		if (useMicrostateData)
			averageMicrostates();

		mat minEnergies = min(*macrostateResidueEnergies, 1);
		mat offsets = minEnergies + log(99) / steepness;
		fitnesses = new mat(nPositions, 20);
		fitnesses->fill(1.0f);
		for (int i = 0; i < nPositions; i++)
		for (int j = 0; j < 20; j++)
		for (int k = 0; k < nMacrostates; k++)
		{
			double f = 1.0f / (1 + exp(steepness * (macrostateResidueEnergies->at(i, j, k) - offsets(i, k))));
			fitnesses->at(i, j) = fitnesses->at(i, j) * (1 - weights[k] + weights[k] * f);
		}
	}

	void Model::calcFrequencies()
	{
		// TODO: check that this actually works
		if (!isFrequenciesCalculated)
		{
			isFrequenciesCalculated = true;
			this->calcFitness();

			this->frequencies = &mat((*fitnesses) / (1 - *fitnesses));	// this works?
			mat sums = sum(*frequencies, 1);
			for (int i = 0; i < nPositions; i++)
				frequencies->at(i) = frequencies->at(i) / sums.at(i);
		}
	}

	void Model::averageMicrostates()
	{
		// TODO: implement this
		if (!areMicrostatesPicked)
		{
			// pick random indices
			for (int i = 0; i < nPositions; i++)
			for (int j = 0; j < nMacrostates; j++)
			for (int k = 0; k < ensembleSize; k++)
				microstatesUsed[i][j][k] = rand() % microstateCounts[i][j]; //TODO: make unique rand ints

			// pick energies
			for (int i = 0; i < nPositions; i++)
				selectedMicrostateEnergies->push_back(cube(20, nMacrostates, ensembleSize));	// unknown values? we shall see
			for (int i = 0; i < nPositions; i++)
			for (int j = 0; j < 20; j++)
			for (int k = 0; k < nMacrostates; k++)
			for (int l = 0; l < ensembleSize; l++)
				selectedMicrostateEnergies->at(i).at(j, k, l) = microstateResidueEnergies->at(i)(j, k, microstatesUsed[i][k][l]);
			areMicrostatesPicked = true;
		}

		//disregard alternative averaging method since it does not make any practical differences
		if (boltzmannTemp == 0)
			//TODO mathy thigns for slicing
			macrostateResidueEnergies = ;
	}

	bool Model::operator==(const Model &other) const
	{
		return this->recovery == other.recovery;
	}

	bool Model::operator <(const Model &other) const
	{
		return this->recovery < other.recovery;
	}

	Model::~Model()
	{
		fitnesses->~Mat();
		frequencies->~Mat();
		macrostateResidueEnergies->~Cube();

		if (microstateCounts)
		{
			for (int i = 0; i < nPositions; i++)
				delete[] microstateCounts[i];
			delete[] microstateCounts;
		}

		if (microstatesUsed)
		{
			for (int i = 0; i < nPositions; i++)
			{
				for (int j = 0; j < nMacrostates; j++)
					delete[] microstatesUsed[i][j];
				delete[] microstatesUsed[i];
			}
			delete[] microstatesUsed;
		}

		microstateResidueEnergies->~vector();
		selectedMicrostateEnergies->~vector();
	}
}
