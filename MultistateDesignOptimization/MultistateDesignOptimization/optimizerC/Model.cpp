#include "stdafx.h"
#include "Model.h"
#include <cstdlib>

namespace OPTIMIZER
{
	Model::Model()
	 :	recovery(0),
		nMacrostates(0),
		nPositions(0),
		positionOffset(0),
		backrubTemp(0),
		boltzmannTemp(0),
		weights(NULL),
		steepness(0),
		useAltAverageingMethod(false),
		useMicrostateData(false)
	{
		this->initializeMembers();
	}

	Model::Model(int nMacrostates, int ensembleSize, double backrubTemp, double boltzmannTemp, int nPositions, int positionOffset)
	 :	recovery(0),
		nMacrostates(nMacrostates),
		ensembleSize(ensembleSize),
		nPositions(nPositions),
		positionOffset(positionOffset),
		backrubTemp(backrubTemp),
		boltzmannTemp(boltzmannTemp),
		weights(NULL),
		steepness(0),
		useAltAverageingMethod(false),
		useMicrostateData(false)
	{
		this->initializeMembers();
	}

	Model::Model(int nMacrostates, int ensembleSize, double backrubTemp, double boltzmannTemp, int nPositions, int positionOffset, bool useMicrostates)
	 :	recovery(0),
	 nMacrostates(nMacrostates),
	 ensembleSize(ensembleSize),
		nPositions(nPositions),
		positionOffset(positionOffset),
		backrubTemp(backrubTemp),
		boltzmannTemp(boltzmannTemp),
		weights(NULL),
		steepness(0),
		useAltAverageingMethod(false),
		useMicrostateData(false)
	{
		this->initializeMembers();
	}

	Model::Model(int nMacrostates, int ensembleSize, double backrubTemp, double boltzmannTemp, int nPositions, int positionOffset, bool useMicrostates, bool useAltAvgMethod)
	 :	recovery(0),
		nMacrostates(nMacrostates),
		ensembleSize(ensembleSize),
		nPositions(nPositions),
		positionOffset(positionOffset),
		backrubTemp(backrubTemp),
		boltzmannTemp(boltzmannTemp),
		weights(NULL),
		steepness(0),
		useAltAverageingMethod(useAltAvgMethod),
		useMicrostateData(useMicrostates)
	{
		this->initializeMembers();
	}

	Model::Model(int nMacrostates, int ensembleSize, double backrubTemp, double boltzmannTemp, double *weights, double steepness, int nPositions, int positionOffset)
	 :	recovery(0),
		nMacrostates(nMacrostates),
		ensembleSize(ensembleSize),
		nPositions(nPositions),
		positionOffset(positionOffset),
		backrubTemp(backrubTemp),
		boltzmannTemp(boltzmannTemp),
		weights(weights),
		steepness(steepness),
		useAltAverageingMethod(false),
		useMicrostateData(false)
	{
		this->initializeMembers();
	}

	Model::Model(int nMacrostates, int ensembleSize, double backrubTemp, double boltzmannTemp, double *weights, double steepness, int nPositions, int positionOffset, bool useMicrostate)
		: recovery(0),
		nMacrostates(nMacrostates),
		ensembleSize(ensembleSize),
		nPositions(nPositions),
		positionOffset(positionOffset),
		backrubTemp(backrubTemp),
		boltzmannTemp(boltzmannTemp),
		weights(weights),
		steepness(steepness),
		useAltAverageingMethod(false),
		useMicrostateData(useMicrostate)
	{
			this->initializeMembers();
	}

	Model::Model(int nMacrostates, int ensembleSize, double backrubTemp, double boltzmannTemp, double *weights, double steepness, int nPositions, int positionOffset, bool useMicrostate, bool useAltaverageMethod)
		: recovery(0),
		nMacrostates(nMacrostates),
		ensembleSize(ensembleSize),
		nPositions(nPositions),
		positionOffset(positionOffset),
		backrubTemp(backrubTemp),
		boltzmannTemp(boltzmannTemp),
		weights(weights),
		steepness(steepness),
		useAltAverageingMethod(useAltaverageMethod),
		useMicrostateData(useMicrostate)
	{
		this->initializeMembers();
	}

	Model::Model(const Model &existing, int ensembleSize, double backrubTemp, double boltzmannTemp, double *weights, double steepness)
	{
		this->recovery = 0;
		this->nPositions = existing.nPositions;
		this->nMacrostates = existing.nMacrostates;
		this->ensembleSize = ensembleSize;
		this->backrubTemp = backrubTemp;
		this->boltzmannTemp = boltzmannTemp;
		this->weights = weights;
		this->steepness = steepness;
		this->fitnesses = mat(nPositions, 20, fill::zeros);
		this->frequencies = mat(nPositions, 20, fill::zeros);
		this->useMicrostateData = existing.useMicrostateData;
		this->areMicrostatesPicked = false;
		this->useAltAverageingMethod = existing.useAltAverageingMethod;

		if (!useMicrostateData && ensembleSize == existing.ensembleSize)
			this->macrostateResidueEnergies = cube((existing.macrostateResidueEnergies));
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

	void Model::initializeMembers()
	{
		this->fitnessCalculated = false;
		this->fitnesses = mat(nPositions, 20, fill::zeros);
		this->frequencies = mat(nPositions, 20, fill::zeros);
		this->macrostateResidueEnergies = cube(nPositions, nMacrostates, 20, fill::zeros);
		this->microstateCounts = NULL;
		this->microstateResidueEnergies = NULL;
		this->microstatesUsed = NULL;
		this->selectedMicrostateEnergies = NULL;
		if (useMicrostateData)
		{
			areMicrostatesPicked = false;
			microstateResidueEnergies = new vector<cube>();
			for (int i = 0; i < nPositions; i++)
				microstateResidueEnergies->push_back(cube(20, nMACROSTATES, 1024, fill::zeros));
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

	void Model::addMacrostateData(int macrostate, int position, double *energies)
	{
		//position -= this->positionOffset; // should already be done by Optimizer
		for (int i = 0; i < 20; i++)
			macrostateResidueEnergies(position, macrostate, i) = energies[i];	// this line gives a write access error...
	}

	void Model::addMicrostateData(int macrostate, int position, double *energies)
	{
		// position -= this->positionOffset; // should already be done by optimizer
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
		
		return &frequencies;
	}

	void Model::calcFitness()
	{
		// TODO: check correctness
		if (useMicrostateData)
			averageMicrostates();

		mat minEnergies = min(macrostateResidueEnergies, 2);
		mat offsets = minEnergies + log(99) / steepness;
		fitnesses = mat(nPositions, 20);
		fitnesses.fill(1.0f);
		for (int i = 0; i < nPositions; i++)
		for (int j = 0; j < 20; j++)
		for (int k = 0; k < nMacrostates; k++)
		{
			double f = 1.0f / (1 + exp(steepness * (macrostateResidueEnergies(i, k, j) - offsets(i, k))));
			fitnesses(i, j) = fitnesses(i, j) * (1 - weights[k] + weights[k] * f);
		}
	}

	void Model::calcFrequencies()
	{
		// TODO: check that this actually works
		if (!isFrequenciesCalculated)
		{
			isFrequenciesCalculated = true;
			this->calcFitness();

			this->frequencies = mat((fitnesses) / (1 - fitnesses));	// this works?
			mat sums = sum(frequencies, 1);
			for (int i = 0; i < nPositions; i++)
				frequencies(i) = frequencies(i) / sums.at(i);
		}
	}

	void Model::averageMicrostates()
	{
		cout << "Average microstates is not working" << endl;
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
		{
			for (int i = 0; i < nMacrostates; i++)
				// TODO: not correct
				macrostateResidueEnergies.slice(i) = min(selectedMicrostateEnergies->at(i), 2);
		}
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
		fitnesses.~Mat();
		frequencies.~Mat();
		macrostateResidueEnergies.~Cube();

		if (microstateCounts)
		{
			for (int i = 0; i < nPositions; i++)
			if (microstateCounts[i])
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

		if (microstateResidueEnergies)
			microstateResidueEnergies->~vector();
		if (selectedMicrostateEnergies)
			selectedMicrostateEnergies->~vector();
	}
}
