#include "stdafx.h"

namespace OPTIMIZER
{
	const int nMACROSTATES = 6;	// TODO: needs to be changed

	enum RESIDUES
	{
		A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y
	};

	enum MACROSTATES	// TODO: needs to be a param
	{
		// TODO: fill this in
	};

	class Model
	{
	public:
		double recovery;				// META-FITNESS similarity to target, assigned by outside similarity measure

		Model();
		Model(int nMacrostates, int ensembleSize, double backrumTemp, double boltzmannTemp, int nPositions, int positionOffset);
		Model(int nMacrostates, int ensembleSize, double backrumTemp, double boltzmannTemp, int nPositions, int positionOffset, bool useMicrostate);
		Model(int nMacrostates, int ensembleSize, double backrumTemp, double boltzmannTemp, int nPositions, int positionOffset, bool useMicrostate, bool useAltAverageMethod);
		Model(int nMacrostates, int ensembleSize, double backrubTemp, double boltzmannTemp, double *weights, double steepness, int nPositions, int positionOffset);
		Model(int nMacrostates, int ensembleSize, double backrubTemp, double boltzmannTemp, double *weights, double steepness, int nPositions, int positionOffset, bool useMicrostate);
		Model(int nMacrostates, int ensembleSize, double backrubTemp, double boltzmannTemp, double *weights, double steepness, int nPositions, int positionOffset, bool useMicrostate, bool useAltaverageMethod);
		Model(const Model &existing, int ensembleSize, double backrubTemp, double boltzmannTemp, double *weights, double steepness);

		void addMacrostateData(MACROSTATES macrostate, int position, float *energies);
		void addMicrostateData(MACROSTATES macrostate, int position, float *energies);

		int getEnsembleSize() { return this->ensembleSize; }
		int getBackrubTemp() { return this->backrubTemp; }
		double *getWeights();
		double getSteepness() { return this->steepness; }
		double **getFrequencies();

		// comparators
		bool operator==(const Model &other) const;
		bool operator!=(const Model &other) const { return !(*this == other); }
		bool operator <(const Model &other) const;
		bool operator >(const Model &other) const { return other < *this; }
		bool operator<=(const Model &other) const { return !(*this > other); }
		bool operator>=(const Model &other) const { return !(*this < other); }

		~Model();

	private:
		int nMacrostates;				// number of macrostates
		int ensembleSize;				// ensemble size
		int nPositions;					// number of position examined
		int positionOffset;				// first position examines, for indexing purposes
		double backrubTemp;				// backrub temperature
		double boltzmannTemp;			// boltzmann averaging temperature
		double *weights;				// weights of each macrostate
		double steepness;				// sigmoid steepness

		// armadillo matrices/cubes
		mat *fitnesses;						// single item
		mat *frequencies;					// single item
		cube *macrostateResidueEnergies;	// single item

		// when microstates are used
		bool useMicrostateData;
		bool useAltAverageingMethod;
		bool isFrequenciesCalculated;
		bool areMicrostatesPicked;
		vector<cube> *microstateResidueEnergies;	// cube[]
		vector<cube> *selectedMicrostateEnergies;	// cube[]
		double **microstateCounts;
		double ***microstatesUsed;
		
		// TODO: fix Matrices - they're all 2D! The numbers are m,n sizes >.<
		//Matrix2d *fitnesses;			// Eigen matrix representations
		//Matrix2d *frequencies;
		//Matrix3d *residueEnergies;

		//double **fitnesses;				// double[position][residue fitness]
		//double **frequencies;			// double[position][residue frequency]
		//double ***residueEnergies;		// double[position][residue][macrostate energy]

		bool fitnessCalculated;

		void calcFitness();
		void calcFrequencies();
		void averageMicrostates();
	};
}