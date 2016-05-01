#pragma once
#include "stdafx.h"

namespace OPTIMIZER
{
	class JensenShannonDistance : public SimilarityMeasure
	{
	public:
		JensenShannonDistance();
		JensenShannonDistance(const SimilarityMeasure &other);
		
		/// <summary>
		/// Initializes a new instance of the <see cref="JensenShannonDistance"/> class.
		/// </summary>
		/// <param name="targetFreqs">The target freqs.</param>
		JensenShannonDistance(mat *targetFreqs);
		
		/// <summary>
		/// Gets the similarity of expFrequencies to the target frequencies.
		/// Will not destroy expFrequencies in the process of doing so. Model's
		/// getFrequenceies was modified to return a pointer to its own copy.
		/// </summary>
		/// <param name="expFrequencies">The expected frequencies.</param>
		/// <returns></returns>
		double getSimilarity(mat *expFrequencies) override;

		~JensenShannonDistance();

		string *toString() override;

	private:		
		/// <summary>
		/// Kernel function for the kernalized implementation of the J-S divergence
		/// </summary>
		/// <param name="">The .</param>
		mat* h(mat*);

		mat *hTargetFreqs;
	};
}