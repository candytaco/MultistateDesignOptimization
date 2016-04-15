#include "stdafx.h"

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
		void writeFrequenciesToFASTA(string *outName, int precision, float **frequencies);
		
		/// <summary>
		/// Writes the best parameters to text.
		/// </summary>
		/// <param name="outName">Name of the out.</param>
		void writeBestParamsToText(string *outName);
		
		/// <summary>
		/// Calculates the parameters identifier.
		/// </summary>
		/// <param name="param1">The param1.</param>
		/// <param name="param2">The param2.</param>
		/// <param name="param3">The param3.</param>
		/// <returns></returns>
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
		float *getBestParameters();
		
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