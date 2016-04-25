#include "Optimizer.h"

using namespace OPTIMIZER;

Optimizer::Optimizer()
{
	Optimizer::Optimizer(0, NULL, false);
}

Optimizer::Optimizer(const Optimizer &existing)
{
	// TODO: implement this
}

Optimizer::Optimizer(int nMacrostates, string *macrostates)
{
	Optimizer::Optimizer(nMacrostates, macrostates, false);
}

Optimizer::Optimizer(int nMacrostates, bool continuousBoltzmann)
{
	Optimizer::Optimizer(nMacrostates, NULL, continuousBoltzmann);
}

Optimizer::Optimizer(int nMacrostates, string *macrostates, bool continuousBoltzmann)
{
	this->models = new map<int, Model>();
	this->optimizationAlgorithm = NULL;
	this->nPositions = 0;
	this->minPosition = 0;
	this->targetFrequencies = NULL;
	this->macrostates = macrostates;
	this->nMacrostates = nMacrostates;
	this->continuousBoltzmann = continuousBoltzmann;
}

void Optimizer::readData(string *inFile)
{

}

void Optimizer::readMicrostateData(string *inFile)
{

}

void Optimizer::writeFrequenciesToFASTA(string *outName, int precision, double **frequencies)
{
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
	if (!outputFile) {
		perror("Error opening FASTA file");
	}

	int nEntries = pow(10, precision);

	int numbers;

	// this is going to create problems, but i'm confused as to what frequencies is. might need to loop through and multiple each element?
	double **numbers = frequencies;
	//int numbers = round(frequencies * nEntries); // uhhh what is frequencies? in the original code it is a numpy array. here it is...?

	vector<int> residueToWrite(nPositions, 0); // i think this should allocate the vector of size nPositions filled with zeros.
	char *residues = "ACDEFGHIKLMNPQRSTVWY";

	// i think the following loop could be optimized, but i don't know how many times we use it.
	// why is there no "i" used within this loop? maybe i don't get what it is doing.
	for (int i = 0; i < nEntries; i++) {
		fprintf(outputFile, "> Null\n");
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

void Optimizer::writeBestParamsToText(string *outName)
{

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
	for (int i = 0; i < nMacrostates; i++) // nMacrostates must be defined earlier...
		fprintf(outputFile, "%.4f", bestVals[4][i]); // this is not the best structure, but change it later when necesary.
	fprintf(outputFile, "\nMatch: %.4f\n", bestVals[5]);

	//TODO: the following needs to be added, but not sure where these are stored. 
	/*
	outfile.write("Algorithm: {:s}\n".format(self.optimizationAlgorithm.__str__()));
	outfile.write("Similarity measure: {:s}\n".format(self.optimizationAlgorithm.similarityMeasure.__str__()));
	outfile.write("Elapsed time: {:s}\n".format(str(self.optimizationAlgorithm.elapsedTime)));
	outfile.close();
	*/


}

int Optimizer::calcParamsID(double param1, double param2, double param3)
{
	return NULL;
}
Model *Optimizer::getModelByParams(double, double, double)
{
	return NULL;
}
