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
		
		/// <summary>
		/// Initializes a new instance of the <see cref="Model"/> class.
		/// </summary>
		/// <param name="nMacrostates">The n macrostates.</param>
		/// <param name="ensembleSize">Size of the ensemble.</param>
		/// <param name="backrubTemp">The backrub temporary.</param>
		/// <param name="boltzmannTemp">The boltzmann temporary.</param>
		/// <param name="weights">The weights.</param>
		/// <param name="steepness">The steepness.</param>
		/// <param name="nPositions">The n positions.</param>
		/// <param name="positionOffset">The position offset.</param>
		/// <param name="useMicrostate">if set to <c>true</c> [use microstate].</param>
		/// <param name="useAltaverageMethod">if set to <c>true</c> [use altaverage method].</param>
		Model(int nMacrostates, int ensembleSize, double backrubTemp, double boltzmannTemp, double *weights, double steepness, int nPositions, int positionOffset, bool useMicrostate, bool useAltaverageMethod);
		
		/// <summary>
		/// Initializes a new instance of the <see cref="Model"/> class.
		/// </summary>
		/// <param name="existing">The existing.</param>
		/// <param name="ensembleSize">Size of the ensemble.</param>
		/// <param name="backrubTemp">The backrub temporary.</param>
		/// <param name="boltzmannTemp">The boltzmann temporary.</param>
		/// <param name="weights">The weights.</param>
		/// <param name="steepness">The steepness.</param>
		Model(const Model &existing, int ensembleSize, double backrubTemp, double boltzmannTemp, double *weights, double steepness);
		
		/// <summary>
		/// Adds the macrostate data.
		/// </summary>
		/// <param name="macrostate">The macrostate.</param>
		/// <param name="position">The position.</param>
		/// <param name="energies">The energies.</param>
		void addMacrostateData(int macrostate, int position, double *energies);
		
		/// <summary>
		/// Adds the microstate data.
		/// </summary>
		/// <param name="macrostate">The macrostate.</param>
		/// <param name="position">The position.</param>
		/// <param name="energies">The energies.</param>
		void addMicrostateData(int macrostate, int position, double *energies);
		
		/// <summary>
		/// Gets the size of the ensemble.
		/// </summary>
		/// <returns></returns>
		const int getEnsembleSize() { return this->ensembleSize; }
		
		/// <summary>
		/// Gets the backrub temperature.
		/// </summary>
		/// <returns></returns>
		const double getBackrubTemp() { return this->backrubTemp; }
		
		/// <summary>
		/// Gets the boltzmann temperature.
		/// </summary>
		/// <returns></returns>
		const double getBoltzmannTemp() { return this->boltzmannTemp; }
		
		/// <summary>
		/// Gets the weights as a deep copy of the internal variables.
		/// </summary>
		/// <returns>double[] of weights</returns>
		double *getWeights();
		
		/// <summary>
		/// Gets the steepness.
		/// </summary>
		/// <returns></returns>
		double getSteepness() { return this->steepness; }
		
		/// <summary>
		/// Returns the frequencies with a deep copy of the object's internal variable
		/// </summary>
		/// <returns>double[position][residue] of frequency values</returns>
		mat *getFrequencies();

		// comparators based on similarity measure for sorting
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
		int **microstateCounts;
		int ***microstatesUsed;
		
		// TODO: fix Matrices - they're all 2D! The numbers are m,n sizes >.<
		//Matrix2d *fitnesses;			// Eigen matrix representations
		//Matrix2d *frequencies;
		//Matrix3d *residueEnergies;

		//double **fitnesses;				// double[position][residue fitness]
		//double **frequencies;			// double[position][residue frequency]
		//double ***residueEnergies;		// double[position][residue][macrostate energy]

		bool fitnessCalculated;
		
		/// <summary>
		/// Calculates fitnesses values.
		/// </summary>
		void calcFitness();
		
		/// <summary>
		/// Calculates the frequencies.
		/// </summary>
		void calcFrequencies();
		
		/// <summary>
		/// Averages the microstates.
		/// </summary>
		void averageMicrostates();
	};
}