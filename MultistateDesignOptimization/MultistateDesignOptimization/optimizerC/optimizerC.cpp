// optimizerC.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Model.h"
#include "SimilarityMeasure.h"
#include "JensenShannonDistance.h"
#include "SearchAlgorithm.h"
#include "CuckooSearch.h"
#include "Optimizer.h"

using namespace OPTIMIZER;


#ifdef BOOP
int _tmain(int argc, _TCHAR* argv[])
{
    string defaultFileName("C:\\Users\\candy_000\\Source\\Repos\\MultistateDesignOptimization\\MultistateDesignOptimization\\MultistateDesignOptimization\\optimizerC\\macrostates.dat");
    string defaultTarget("C:\\Users\\candy_000\\Source\\Repos\\MultistateDesignOptimization\\MultistateDesignOptimization\\MultistateDesignOptimization\\optimizerC\\targetFreqs.dat");
    
    string fileName, targetFreqs;
    if (argc > 1)
    {
        fileName = string((char*)argv[1]);
        targetFreqs = string((char*)argv[2]);
    }
    else
    {
        fileName = defaultFileName;
        targetFreqs = defaultTarget;
    }
    
    Optimizer *optimizer = new Optimizer(6, false);
    optimizer->readTargetFrequencies(&targetFreqs);
    optimizer->readData(&fileName);
    CuckooSearch *cs = new CuckooSearch(6, optimizer->getModels(), new JensenShannonDistance(optimizer->getTargetFreqs()), 32, 1, 0.2);
    cs->setMaxIterations(1024);
    // search parameters
    int ensembleSizes[] = { 20, 50 };
    double backrubTemps[] = { 0.3, 0.6, 0.9, 1.2, 1.5, 1.8 };
    double boltzmannTemps[] = { 1, 5, -1, 0 };
    double steepness[] = { 1.0, 7.0 };
    double weightMins[] = { 0, 0, 0, 0, 0, 0 };
    double weightMaxs[] = { 1, 1, 1, 1, 1, 1 };
    cs->setParameterBounds(ensembleSizes, 2, backrubTemps, 6, boltzmannTemps, 4, steepness, weightMins, weightMaxs);
    optimizer->useAlgorithm(cs);
    optimizer->optimize();
    optimizer->writeFrequenciesToFASTA(&string("test1.fasta"));
    optimizer->writeBestParamsToText(&string("test1.txt"));
    return 0;
}

#else
#include <sys/time.h>
#include <time.h>

double read_timer( )
{
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

int main( int argc, char* argv[] )
{
    string defaultFileName("./microstates.dat");
    string defaultTarget("./targetFreqs.dat");
    double simulation_time = read_timer( );

    string fileName, targetFreqs;
    if (argc > 1)
    {
        fileName = string((char*)argv[1]);
        targetFreqs = string((char*)argv[2]);
    }
    else
    {
        fileName = defaultFileName;
        targetFreqs = defaultTarget;
    }
    
    Optimizer *optimizer = new Optimizer(5, true);
    optimizer->readTargetFrequencies(&targetFreqs);
    optimizer->readMicrostateData(&fileName);
    int thisPopulationSize = 128;
    CuckooSearch *cs = new CuckooSearch(6, optimizer->getModels(), new JensenShannonDistance(optimizer->getTargetFreqs()), thisPopulationSize, 1, 0.2, true);
    cs->setMaxIterations(1024);
    // search parameters
    int ensembleSizes[] = { 20, 50, 75, 100 };
    double backrubTemps[] = { 0.3, 0.6, 0.9, 1.2, 1.5, 1.8 };
    double boltzmannTemps[] = {-1, 5};
    double steepness[] = { 1.0, 7.0 };
    double weightMins[] = { 0, 0, 0, 0, 0 };
    double weightMaxs[] = { 1, 1, 1, 1, 1 };
    cs->setParameterBounds(ensembleSizes, 4, backrubTemps, 6, boltzmannTemps, 4, steepness, weightMins, weightMaxs);
    optimizer->useAlgorithm(cs);
    optimizer->optimize();
    string testf = "test2.fasta";
    string testt = "test2.txt";
    simulation_time = read_timer( ) - simulation_time;
    printf( "populationSize = %d,iterations = %d, simulation time = %g seconds", thisPopulationSize,1024,simulation_time);
    optimizer->writeFrequenciesToFASTA(&testf);
    optimizer->writeBestParamsToText(&testt);
    return 0;
}
#endif
