COMPILER=gcc
FLAGS= -Wall -O3

test : optimizerC.cpp stdafx.h Model.h SimilarityMeasure.h JensenShannonDistance.h SearchAlgorithm.h CuckooSearch.h Optimizer.h
	$(COMPILER) -o $(FLAGS) optimizerC.cpp stdafx.obj Model.obj SimilarityMeasure.obj JensenShannonDistance.obj SearchAlgorithm.obj CuckooSearch.obj Optimizer.obj

stdafx.obj : stdafx.h
	$(COMPILER) -c $(FLAGS) stdafx.cpp

Model.obj : stdafx.h Model.h
	$(COMPILER) -c $(FLAGS) stdafx.obj Model.cpp

SimilarityMeasure.obj : stdafx.h SimilarityMeasure.h
	$(COMPILER) -c $(FLAGS) stdafx.obj SimilarityMeasure.cpp

JensenShannonDistance.obj : stdafx.h SimilarityMeasure.h JensenShannonDistance.h
	$(COMPILER) -c $(FLAGS) stdafx.obj SimilarityMeassure.obj JensenShannonDistance.cpp

SearchAlgorithm.obj : stdafx.h Model.h SimilarityMeasure.h SearchAlgorithm.h
	$(COMPILER) -c $(FLAGS) stdafx.obj Model.obj SimilarityMeasure.obj SearchAlgorithm.cpp

CuckooSearch.obj : stdafx.h Model.h SimilarityMeasure.h SearchAlgorithm.h CuckooSearch.h
	$(COMPILER) -c $(FLAGS) stdafx.obj Model.obj SimilarityMeasure.obj SearchAlgorithm.obj CuckooSearch.cpp

Optimizer.obj: stdafx.h Model.h SimilarityMeasure.h SearchAlgorithm.h Optimizer.h
	$(COMPILER) -c $(FLAGS) stdafx.obj Model.obj SimilarityMeasure.obj SearchAlgorithm.obj Optimizer.cpp

clean :
	rm *.obj