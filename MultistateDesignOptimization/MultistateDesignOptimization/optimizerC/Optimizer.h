#pragma once
#include "stdafx.h"
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include <list>
#include <string>

namespace OPTIMIZER
{
	int calcParamsID(double param1, double param2, double param3);

	class Optimizer
	{
	public:
		string *macrostates;	// list of macrostates
				
		/// <summary>
		/// Initializes a new instance of the <see cref="Optimizer"/> class.
		/// </summary>
		Optimizer();

		/// <summary>
		/// Initializes a new instance of the <see cref="Optimizer"/> class.
		/// </summary>
		/// <param name="existing">The existing.</param>

		Optimizer(const Optimizer &existing);		
		/// <summary>
		/// Copy constructor for <see cref="Optimizer"/> class.
		/// </summary>
		/// <param name="nMacrostates">The n macrostates.</param>
		/// <param name="macrostates">The macrostates.</param>
		Optimizer(int nMacrostates, string *macrostates);
		
		/// <summary>
		/// Initializes a new instance of the <see cref="Optimizer"/> class.
		/// </summary>
		/// <param name="nMacrostates">The n macrostates.</param>
		/// <param name="continuousBoltzmann">if set to <c>true</c> [continuous boltzmann].</param>
		Optimizer(int nMacrostates, bool continuousBoltzmann);
		
		/// <summary>
		/// Initializes a new instance of the <see cref="Optimizer"/> class.
		/// </summary>
		/// <param name="nMacrostates">The n macrostates.</param>
		/// <param name="macrostates">The macrostates.</param>
		/// <param name="continuousBoltzmann">if set to <c>true</c> [continuous boltzmann].</param>
		Optimizer(int nMacrostates, string *macrostates, bool continuousBoltzmann);
		
		/// <summary>
		/// Reads macrostate data from a .dat file that I created
		/// </summary>
		/// <param name="inFile">The in file.</param>
		void readData(string *inFile);
		
		/// <summary>
		/// Reads the microstate data from a .dat file that I created
		/// </summary>
		/// <param name="inFile">The in file.</param>
		void readMicrostateData(string *inFile);
		
		/// <summary>
		/// Reads the target frequencies.
		/// </summary>
		/// <param name="inFile">The in file.</param>
		void readTargetFrequencies(string *inFile);
		
		/// <summary>
		/// Writes the frequencies to fasta.
		/// </summary>
		/// <param name="outName">Name of the out.</param>
		void writeFrequenciesToFASTA(string *outName) { this->writeFrequenciesToFASTA(outName, 3, this->optimizationAlgorithm->getBestFrequencies()); }
		
		/// <summary>
		/// Writes the frequencies to fasta.
		/// </summary>
		/// <param name="outName">Name of the out.</param>
		/// <param name="precision">The precision.</param>
		void writeFrequenciesToFASTA(string *outName, int precision) { this->writeFrequenciesToFASTA(outName, precision, this->optimizationAlgorithm->getBestFrequencies()); }
		
		/// <summary>
		/// Writes the frequencies to fasta.
		/// </summary>
		/// <param name="outName">Name of the out.</param>
		/// <param name="precision">The precision.</param>
		/// <param name="frequencies">The frequencies.</param>
		void writeFrequenciesToFASTA(string *outName, int precision, double **frequencies); 
        		
		/// <summary>
		/// Writes the best parameters to text.
		/// </summary>
		/// <param name="outName">Name of the out.</param>
		void writeBestParamsToText(string *outName);
		
		/// <summary>
		/// Gets the model by parameters.
		/// </summary>
		/// <param name="">The .</param>
		/// <param name="">The .</param>
		/// <param name="">The .</param>
		/// <returns></returns>
		Model *getModelByParams(double, double, double);
				
		/// <summary>
		/// Uses the algorithm.
		/// </summary>
		/// <param name="">The .</param>
		void useAlgorithm(SearchAlgorithm *);
		
		/// <summary>
		/// Starts optimization process.
		/// </summary>
		void optimize();
		
		/// <summary>
		/// Gets the best parameters.
		/// </summary>
		/// <returns></returns>
		double *getBestParameters();
		// this needs to contain the following. i'm going to assume that they are each a float for now, but weights needs to be multiple floats!?
        // for the write functions, it is assumed these are stored in an array of floats.
        
        /* Keys:
        0: 'ensembleSize'
        1: 'backrubTemp'
        2: 'boltzmannTemp'
        3: 'steepness'
        4: 'weights' (actually an array of floats size nMacrostates)
        5: 'match'
        */
        
		/// <summary>
		/// Gets the best frequencies.
		/// </summary>
		/// <returns></returns>
		mat *getBestFrequencies();
		
		/// <summary>
		/// String representation of the optimizer.
		/// </summary>
		/// <returns></returns>
		string *toString();

		map<int, Model> *getModels() const;

		mat *getTargetFreqs() const;
		
		/// <summary>
		/// Finalizes an instance of the <see cref="Optimizer"/> class.
		/// </summary>
		~Optimizer();

	private:
		map<int, Model> *models;
		SearchAlgorithm *optimizationAlgorithm;
		int nPositions;
		int minPosition;
		mat *targetFrequencies;
		int nMacrostates;
		bool continuousBoltzmann;
	};
}