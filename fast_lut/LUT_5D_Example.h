/*
 * LUT_5D_Example.h
 *
 *  Created on: Jul 20, 2021
 *      Author: bgregor
 */

#ifndef LUT_5D_EXAMPLE_H_
#define LUT_5D_EXAMPLE_H_

#include "LookupTable.h"

/* This is an example of a 5D lookup table.  As the values are
 * computed on-the-fly the compute_value method must be overridden.
 *
 * As an example, this stores a sine calculation.  A phase parameter
 * is provided as a constructor argument and is used in the compute_value
 * implementation.
 */

class LUT_5D_Example : public LookupTable {
public:
	// Initialize with the size of each dimension and an extra parameter.
	LUT_5D_Example(const int dim0, const int dim1, const int dim2, const int dim3, const int dim4, const double phase) ;

	// Default destructor is fine, everything is handled via containers.
	virtual ~LUT_5D_Example() = default ;

	// This variable is just used as a demo, when set to true the compute_value
	// function will print a message.
	bool m_verbose_computation = false ;

protected:
	double m_phase ; // Internal storage of a phase.

	// This is the override of the compute_value method. It MUST:
	// 	  **  take an initializer_list of ints (matrix index values) as the argument
	//    **  be virtual,
	//    **  be declared here.
	virtual double compute_value(const std::initializer_list<int> inds) ;

};


#endif /* LUT_5D_EXAMPLE_H_ */
