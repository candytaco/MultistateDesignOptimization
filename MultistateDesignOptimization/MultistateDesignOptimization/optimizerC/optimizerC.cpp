// optimizerC.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Optimizer.h"
#include "CuckooSearch.h"

using namespace OPTIMIZER;

int _tmain(int argc, _TCHAR* argv[])
{
	string fileName((char*)argv[1]);
	string targetFreqs((char*)argv[2]);
	Optimizer *optimizer = new Optimizer(6, false);
	optimizer->readTargetFrequencies(&targetFreqs);
	optimizer->readData(&fileName);
	CuckooSearch *cs = new CuckooSearch(6, optimizer->getModels(), );
	return 0;
}

