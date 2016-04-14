#include "stdafx.h"

namespace OPTIMIZER
{
	class Optimizer
	{
	public:
		string *macrostates;	// list of macrostates

		Optimizer();
		Optimizer(Optimizer *existing);
		Optimizer(int nMacrostates, string *macrostates);
		Optimizer(int nMacrostates, bool continuousBoltzmann);
		Optimizer(int nMacrostates, string *macrostates, bool continuousBoltzmann);

		void readData(string *inFile);
		void readMicrostateData(string *inFile);
		void writeFrequenciesToFASTA(string *outName);
		void writeFrequenciesToFASTA(string *outName, int precision);
		void writeFrequenciesToFASTA(string *outName, int precision, mxArray *frequencies);
		void writeBestParamsToText(string *outName);

		int calcParamsID(double param1, double param2, double param3);
		Model *getModelByParams(double, double, double);
		
		void useAlgorithm(SearchAlgorithm);
		void optimize();

		float *getBestParameters();
		mxArray *getBestFrequencies();

		string *toString();

		~Optimizer();

	private:
		map<int, Model> models;
		SearchAlgorithm optimizationAlgorithm;
		int nPositions;
		int minPosition;
		mxArray *targetFrequencies;
		string *macrostates;
		int nMacrostates;
		bool continuousBoltzmann;
	};
}