/*
 * LookupTable.h
 *
 *  Created on: Jul 13, 2021
 *      Author: bgregor
 */

#ifndef LOOKUPTABLE_H_
#define LOOKUPTABLE_H_

#include <initializer_list>
#include <string>
#include <cstdint>
#include <vector>
#include <cassert>
#include "tsl/hopscotch_map.h"


class LookupTable {
public:
	// No empty constructor
	LookupTable() = delete;

	// This takes any number of ints which are the dimensions of the
	// matrix to construct.  The initializer list handles any number
	// of arguments.
	LookupTable(const std::initializer_list<int> dims, bool reserve_size=false) ;

	virtual ~LookupTable() = default ;

	// This is the "real" operator() method. Its arguments are in the form
	// of an intializer list. The overrides below just call this.
	virtual double operator ()(const std::initializer_list<int> dims) ;

	// Provide overrides for operator() 1D thru 5D.
	virtual double operator() (const int i0) ;
	virtual double operator() (const int i0, const int i1) ;
	virtual double operator() (const int i0, const int i1, const int i2) ;
	virtual double operator() (const int i0, const int i1, const int i2, const int i3) ;
	virtual double operator() (const int i0, const int i1, const int i2, const int i3, const int i4) ;


	// How many values are actually stored in the table?
	uint64_t count() const ;

	// What % of the total potential storage is used in the table?
	double fill_percent() const ;

protected:
	std::vector<int> m_dims ; // Storage for the matrix dimensions
	// cumulative product of the dimensions. This is needed for the
	// sub2ind calculation, it's pre-calculated at initialization time
	// for performance.
	std::vector<uint64_t> m_dims_cp ;
	uint64_t m_num_elems ; // Total number of elements that could be stored.
	unsigned int m_n_dims ; // Number of dimensions

	// After testing, the Tessil hopscotch_map (open source with the MIT license) hash
	// table was essentially tied with the STL std::map for insert and lookup speed.
	// According to the Intel Advisor it made better use of the L2 cache than the std::map,
    //so let's use it here.
	// Source:  https://github.com/Tessil/hopscotch-map
	// benchmarks: https://tessil.github.io/2016/08/29/benchmark-hopscotch-map.html
	// The hopscotch_map is also pretty memory efficient.
	// Only the header files are needed.  All N-dimensional lookup tables
	// are stored as a lookup using a 64-bit integer as the lookup value.
	// This is computed from the indices, see sub2ind() below.
	// The hopscotch_map also has built-in serialization so it could be saved and read
	// from disk if desired.
	tsl::hopscotch_map<uint64_t, double> m_hash ;

	// Compute the linear vector index from a set of indices.
	// named sub2ind after the function that does this in Matlab.
    uint64_t sub2ind(const std::initializer_list<int> inds) ;

    // This is to be overridden in the child classes.
    virtual double compute_value(const std::initializer_list<int>) ;
};

#endif /* LOOKUPTABLE_H_ */
