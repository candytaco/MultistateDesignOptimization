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
		/// Reads the data.
		/// </summary>
		/// <param name="inFile">The in file.</param>
		void readData(string *inFile);
		
		/// <summary>
		/// Reads the microstate data.
		/// </summary>
		/// <param name="inFile">The in file.</param>
		void readMicrostateData(string *inFile);
		
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
        void writeFrequenciesToFASTA(string *outName, int precision, double **frequencies) {
            /// <summary>
            /// Writes the best parameters to text.
            /// </summary>
            /// <param name="outName">Name of the out.</param>
            
            // 4/24 has issues ... what is the format for freqeuncies? frequencies must bme multiplied by nPositions (not currently happening).
            // what size is frequencies? is it a float by a float?
            // When do these get written? is everything writing to the same file or different files?
            // is this only ever called once or is it called multiple times?
            // nPositions : is this global?
            
            // TODO: check and make sure that the output file doesn't already end in .fasta!
            string outFileName = outName->append(".fasta"); // might not be necessary, but do it just in case.
            FILE *outputFile = fopen(outFileName.c_str(), "w");
            if(!outputFile) {
                perror("Error opening FASTA file");
            }
            
            int nEntries = pow(10, precision);
            
			int numbers;
            
            // this is going to create problems, but i'm confused as to what frequencies is. might need to loop through and multiple each element?
            double **numbers = frequencies;
            //int numbers = round(frequencies * nEntries); // uhhh what is frequencies? in the original code it is a numpy array. here it is...?
            
            vector<int> residueToWrite(nPositions,0); // i think this should allocate the vector of size nPositions filled with zeros.
            char *residues = "ACDEFGHIKLMNPQRSTVWY";
            
            // i think the following loop could be optimized, but i don't know how many times we use it.
            // why is there no "i" used within this loop? maybe i don't get what it is doing.
            for (int i = 0; i < nEntries; i++) {
                fprintf(outputFile,"> Null\n");
                for (int j = 0; j < nPositions; j++) {
                    while ((numbers[j][residueToWrite[j]] == 0) && residueToWrite[j] < 19)
                        residueToWrite[j]++;
                    numbers[j][residueToWrite[j]]--;
                    fprintf(outputFile, &residues[residueToWrite[j]]);
                }
                fprintf(outputFile, "\n");
            }
            fclose(outputFile);
            }
        

        void writeBestParamsToText(string *outName) {
		
		/// <summary>
		/// Calculates the parameters identifier.
		/// </summary>
		/// <param name="param1">The param1.</param>
		/// <param name="param2">The param2.</param>
		/// <param name="param3">The param3.</param>
		/// <returns></returns>
        
            char outFileName = outname + '.txt';
            FILE *outputFile = fopen(outFileName, "w");
			double *bestVals = getBestParameters();
            // for the write functions, it is assumed getBestParameters are stored in an array of floats.
            /* Keys:
             0: 'ensembleSize'
             1: 'backrubTemp'
             2: 'boltzmannTemp'
             3: 'steepness'
             4: 'weights'
             5: 'match'
             */
            
            fprintf(outputFile, "Ensemble Size: %.0f\n", bestVals[0]);
            fprintf(outputFile, "Backrub temperature: %.1f\n", bestVals[1]);
            if (bestVals[2] > 0)
                fprintf(outputFile, "Backrub temperature: %.9f\n", bestVals[2]);
            else if (bestVals[2] == 0)
                fprintf(outputFile, "Backrub temperature: mean\n ");
            else
                fprintf(outputFile, "Backrub temperature: inf\n ");
            fprintf(outputFile, "Steepness:  %.9f\n", bestVals[3]);
            fprintf(outputFile, "Weights:  ");
            for (int i=0; i < nMacrostates; i++) // nMacrostates must be defined earlier...
                fprintf(outputFile, "%.4f", bestVals[4][i]); // this is not the best structure, but change it later when necesary.
            fprintf(outputFile,"\nMatch: %.4f\n", bestVals[5]);
            
            //TODO: the following needs to be added, but not sure where these are stored. 
            /*
            outfile.write("Algorithm: {:s}\n".format(self.optimizationAlgorithm.__str__()));
            outfile.write("Similarity measure: {:s}\n".format(self.optimizationAlgorithm.similarityMeasure.__str__()));
            outfile.write("Elapsed time: {:s}\n".format(str(self.optimizationAlgorithm.elapsedTime)));
            outfile.close();
             */
            
            
        }
		int calcParamsID(double param1, double param2, double param3);
		
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
		void useAlgorithm(SearchAlgorithm);
		
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
		string *macrostates;
		int nMacrostates;
		bool continuousBoltzmann;
	};
}