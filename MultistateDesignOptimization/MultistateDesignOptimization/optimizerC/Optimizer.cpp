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

void Optimizer::writeFrequenciesToFASTA(string *outName, int precision, float **frequencies)
{

}

void Optimizer::writeBestParamsToText(string *outName)
{

}

int Optimizer::calcParamsID(double param1, double param2, double param3)
{
	return NULL;
}
Model *Optimizer::getModelByParams(double, double, double)
{
	return NULL;
}
