#include "stdafx.h"

namespace OPTIMIZER
{
	class SearchAlgorithm
	{
	public:		
		/// <summary>
		/// Initializes a new instance of the <see cref="SearchAlgorithm"/> class.
		/// </summary>
		SearchAlgorithm();

		SearchAlgorithm(map<int, Model> models, SimilarityMeasure similarityMeasure);
		
		/// <summary>
		/// Initializes a new instance of the <see cref="SearchAlgorithm"/> class.
		/// </summary>
		/// <param name="models">The models.</param>
		/// <param name="similarityMeasure">The similarity measure.</param>
		/// <param name="continuousBoltzmann">if set to <c>true</c> [continuous boltzmann].</param>
		SearchAlgorithm(map<int, Model> models, SimilarityMeasure similarityMeasure, bool continuousBoltzmann);
		
		/// <summary>
		/// Sets the similarity measure.
		/// </summary>
		/// <param name="similarityMeasure">The similarity measure.</param>
		void setSimilarityMeasure(SimilarityMeasure similarityMeasure);
		
		/// <summary>
		/// Sets the maximum iterations.
		/// </summary>
		/// <param name="iter">The iter.</param>
		void setMaxIterations(int iter);
		
		/// <summary>
		/// Sets the continuous boltzmann.
		/// </summary>
		/// <param name="yes">if set to <c>true</c> [yes].</param>
		void setContinuousBoltzmann(bool yes);
		
		/// <summary>
		/// Gets the continuous boltzmann.
		/// </summary>
		/// <returns></returns>
		bool getContinuousBoltzmann();
		
		/// <summary>
		/// Sets the parameter bounds.
		/// </summary>
		/// <param name="ensembleSizes">The ensemble sizes.</param>
		/// <param name="backrubTemps">The backrub temps.</param>
		/// <param name="boltzmannTemps">The boltzmann temps.</param>
		/// <param name="steepnessRange">The steepness range.</param>
		/// <param name="weightMins">The weight mins.</param>
		/// <param name="weightMaxs">The weight maxs.</param>
		void setParameterBounds(int *ensembleSizes, double *backrubTemps, double *boltzmannTemps, double *steepnessRange, double *weightMins, double *weightMaxs);
		
		/// <summary>
		/// Sets the searchparameters.
		/// </summary>
		/// <param name="ensemble">if set to <c>true</c> [ensemble].</param>
		/// <param name="backrub">if set to <c>true</c> [backrub].</param>
		/// <param name="boltzmann">if set to <c>true</c> [boltzmann].</param>
		/// <param name="steepness">if set to <c>true</c> [steepness].</param>
		/// <param name="weights">The weights.</param>
		void setSearchparameters(bool ensemble, bool backrub, bool boltzmann, bool steepness, bool *weights);
		
		/// <summary>
		/// Starts the search.
		/// </summary>
		virtual void startSearch();		

		/// <summary>
		/// Iterates this instance.
		/// </summary>
		virtual void iterate();
		
		/// <summary>
		/// Gets the model by parameters.
		/// </summary>
		/// <param name="">The .</param>
		/// <param name="">The .</param>
		/// <param name="">The .</param>
		/// <returns></returns>
		Model *getModelByParams(double, double, double);
		
		/// <summary>
		/// Gets the best parameters.
		/// </summary>
		/// <returns></returns>
		double *getBestParams();
		
		/// <summary>
		/// Gets the best frequencies.
		/// </summary>
		/// <returns></returns>
		double **getBestFrequencies();
				
		/// <summary>
		/// The suppress outputs
		/// </summary>
		bool suppressOutputs;
		
		/// <summary>
		/// Finalizes an instance of the <see cref="SearchAlgorithm"/> class.
		/// </summary>
		~SearchAlgorithm();

	protected:
		
		/// <summary>
		/// Bounds the check boltzmann.
		/// </summary>
		/// <param name="newBoltzmann">The new boltzmann.</param>
		void boundCheckBoltzmann(double *newBoltzmann);
		
		/// <summary>
		/// Bounds the check steepness.
		/// </summary>
		/// <param name="newSteep">The new steep.</param>
		void boundCheckSteepness(double *newSteep);
		
		/// <summary>
		/// Bounds the check weights.
		/// </summary>
		/// <param name="weights">The weights.</param>
		void boundCheckWeights(double *weights);
		
		/// <summary>
		/// Records the best parameters.
		/// </summary>
		virtual void recordBestParams();

		map<int, Model> models;
		SimilarityMeasure *similarityMeasure;
		int maxIterations;

		int *ensembleSizes;
		double *backrubTemps;
		double *boltzmannTemps;
		bool continuousBoltzmann;
		double *steepnessRange;
		double *weightMins;
		double *weightMaxs;

		bool searchEnsemble;
		bool searchBackrub;
		bool searchBoltzmann;
		bool searchSteepness;
		bool *searchWeights;

		int bestEnsembleSize;
		double bestBackrubTemp;
		double bestBoltzmannTemp;
		double bestSteepness;
		double *bestWeights;
		double bestMatchVal;
		mat *bestFrequencies;
	};
}