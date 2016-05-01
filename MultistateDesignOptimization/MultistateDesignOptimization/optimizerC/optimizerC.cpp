// optimizerC.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

using namespace OPTIMIZER;

int _tmain(int argc, _TCHAR* argv[])
{
	string fileName((char*)argv[1]);
	string targetFreqs((char*)argv[2]);
	Optimizer *optimizer = new Optimizer(6, false);
	optimizer->readTargetFrequencies(&targetFreqs);
	optimizer->readData(&fileName);
	CuckooSearch *cs = new CuckooSearch(6, optimizer->getModels(), new JensenShannonDistance(optimizer->getTargetFreqs), 32, 1, 0.2);
	optimizer->useAlgorithm(cs);
	optimizer->writeFrequenciesToFASTA(&string("test1.fasta"));
	optimizer->writeBestParamsToText(&string("test1.txt"));
	return 0;
}

