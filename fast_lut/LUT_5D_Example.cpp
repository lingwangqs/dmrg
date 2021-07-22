/*
 * LUT_5D_Example.cpp
 *
 * This shows how to create a 5D lookup table that computes and caches values
 * as they are requested.
 *
 * This is safe to use in OpenMP parallel regions.
 *
 *  Created on: Jul 20, 2021
 *      Author: bgregor
 */


#include "LUT_5D_Example.h"
#include <cmath>
#include <iostream>
using namespace std ;

// Call the parent initializer with an intializer_list, then store the phase.
LUT_5D_Example::LUT_5D_Example(const int dim0, const int dim1, const int dim2,
		const int dim3, const int dim4, const double phase) :
		LookupTable({dim0, dim1, dim2, dim3, dim4}), m_phase(phase) { }


// The inline flag is recommended.
inline double LUT_5D_Example::compute_value(const std::initializer_list<int> inds) {
	// For convenient access to the indices, extract them from the initializer_list.
	// The idea is that this compute_value is both computationally expensive and
	// that there are far fewer values needed in the program compared with what could
	// be computed for a simple lookup table with all values filled in. The convenience
	// is worth the miniscule computational cost here.
	int indices[m_n_dims] ;
	int i = 0 ;
	for (auto const &elem : inds) {
		indices[i] = elem ;
		++i ;
	}
	// Compute something with a sine function from the index values.  This is just a placeholder.
	double value =  indices[3] * sin(indices[1] * indices[4] + m_phase) + indices[2] - indices[0] ;
	if (m_verbose_computation) {
		cout << "LUT_5D_Example::compute_value -> Value computed at [" << indices[0] << "," << indices[1] << "," ;
		cout << indices[2] << "," << indices[3] << "," << indices[4] << "]: " << value << endl ;
	}
	return value ;
}
