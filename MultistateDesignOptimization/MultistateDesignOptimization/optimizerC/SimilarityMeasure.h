#pragma once
#include "stdafx.h"

namespace OPTIMIZER
{
	class SimilarityMeasure
	{
	public:		
		/// <summary>
		/// Initializes a new instance of the <see cref="SimilarityMeasure"/> class.
		/// </summary>
		SimilarityMeasure();
		
		/// <summary>
		/// Initializes a new instance of the <see cref="SimilarityMeasure"/> class.
		/// </summary>
		/// <param name="other">The other.</param>
		SimilarityMeasure(const SimilarityMeasure& other);
		
		/// <summary>
		/// Initializes a new instance of the <see cref="SimilarityMeasure"/> class.
		/// </summary>
		/// <param name="targetFreqs">The target freqs.</param>
		SimilarityMeasure(mat *targetFreqs);
		
		/// <summary>
		/// Sets the target frequencies.
		/// </summary>
		/// <param name="targetFreqs">The target freqs.</param>
		void setTargetFrequencies(mat *targetFreqs);

//		virtual SimilarityMeasure *clone();	// TODO: not needed because copy constructor
		
		/// <summary>
		/// Gets the similarity.
		/// </summary>
		/// <param name="expFrequencies">The exp frequencies.</param>
		/// <returns></returns>
		virtual double getSimilarity(mat *expFrequencies) = 0;
		
		/// <summary>
		/// Gets string descriptor
		/// </summary>
		/// <returns></returns>
		virtual string *toString() = 0;
		
		/// <summary>
		/// Finalizes an instance of the <see cref="SimilarityMeasure"/> class.
		/// </summary>
		~SimilarityMeasure();

	protected:
		mat *targetFrequencies;
	};
}