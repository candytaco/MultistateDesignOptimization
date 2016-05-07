COMPILER=CC
FLAGS= -I/global/homes/j/jfeather/clibs/usr/include -I/usr/common/usg/boost/1.54.0-20160127/intel/include/ -L/global/homes/j/jfeather/clibs/usr/lib64 -Wall -O3 -std=c++11 -fpermissive -fopenmp
LIBRARY = -larmadillo -llapack -lblas

TARGETS =  stdafx.o Model.o SimilarityMeasure.o JensenShannonDistance.o SearchAlgorithm.o CuckooSearch.o   Optimizer.o test

all: $(TARGETS)

test : optimizerC.cpp stdafx.h Model.h SimilarityMeasure.h JensenShannonDistance.h SearchAlgorithm.h CuckooSearch.h Optimizer.h
	$(COMPILER) -o optimizerC.o $(FLAGS) optimizerC.cpp stdafx.o Model.o SimilarityMeasure.o JensenShannonDistance.o SearchAlgorithm.o CuckooSearch.o Optimizer.o

stdafx.o : stdafx.h
	$(COMPILER) -c $(FLAGS) stdafx.cpp

Model.o : stdafx.h Model.h Optimizer.h
	$(COMPILER) -c $(FLAGS) stdafx.o Model.cpp

SimilarityMeasure.o : stdafx.h SimilarityMeasure.h
	$(COMPILER) -c $(FLAGS) stdafx.o SimilarityMeasure.cpp

JensenShannonDistance.o : stdafx.h SimilarityMeasure.h JensenShannonDistance.h
	$(COMPILER) -c $(FLAGS) stdafx.o SimilarityMeasure.o JensenShannonDistance.cpp

SearchAlgorithm.o : stdafx.h Model.h SimilarityMeasure.h SearchAlgorithm.h
	$(COMPILER) -c $(FLAGS) stdafx.o Model.o SimilarityMeasure.o SearchAlgorithm.cpp

CuckooSearch.o : stdafx.h Model.h SimilarityMeasure.h SearchAlgorithm.h CuckooSearch.h
	$(COMPILER) -c $(FLAGS) stdafx.o Model.o SimilarityMeasure.o SearchAlgorithm.o CuckooSearch.cpp

Optimizer.o: stdafx.h Model.h SimilarityMeasure.h SearchAlgorithm.h Optimizer.h
	$(COMPILER) -c $(FLAGS) stdafx.o Model.o SimilarityMeasure.o SearchAlgorithm.o Optimizer.cpp

clean :
	rm *.o
