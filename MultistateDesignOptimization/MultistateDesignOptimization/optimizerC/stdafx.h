// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"

#include <stdio.h>
#include <tchar.h>
#include <string>
//#include "matrix.h"
#include <map>
#include <list>
#include <vector>
#include <cmath>
#include <random>
#include <armadillo>
#include <boost/math/distributions/normal.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <fstream>
#include <iostream>

using namespace std;
using namespace arma;

#include "Model.h"
#include "SimilarityMeasure.h"
#include "SearchAlgorithm.h"
#include "CuckooSearch.h"
#include "Optimizer.h"
#include "JensenShannonDistance.h"