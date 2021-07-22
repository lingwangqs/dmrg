/*
 * LookupTable.cpp
 *
 *  Created on: Jul 13, 2021
 *      Author: bgregor
 */

#include "LookupTable.h"


#include <algorithm>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <sstream>

LookupTable::LookupTable(const std::initializer_list<int> dims, bool reserve_size) {
	// Store the dimensions
	m_dims.insert(m_dims.begin(),dims) ;
	m_n_dims = m_dims.size() ;
	// Calculate the cumulative product, this is needed for the sub2ind method.
	m_dims_cp = std::vector<uint64_t>(m_dims.begin(), m_dims.end()) ;
	for (unsigned int i = 1 ; i < m_n_dims ; ++i) {
		m_dims_cp[i] *= m_dims_cp[i-1] ;
	}
	// Get the total possible number of elements to store
	m_num_elems = std::accumulate(m_dims.begin(), m_dims.end(), static_cast<uint64_t>(1), std::multiplies<uint64_t>());
	// Optional - if requested, the hash table can
}


inline double LookupTable::operator()(const std::initializer_list<int> inds) {
	// Look up a value. Compute it if necessary.
	uint64_t ind = sub2ind(inds) ;
	// For use with multiple threads - the "omp critical" forces serialization on
	// the read and possible update to m_hash.
	tsl::hopscotch_map<uint64_t,double>::iterator check ;
	double value ;
#pragma omp critical (hash_access)
	{
		// Check to see if it's in the table already.
		check = m_hash.find(ind);
		if (check != m_hash.end()) {
			value = check->second ;
		} else {
			// Value was not found. Compute it, insert it,
			// and return it. The "real" compute_value will
			// be implemented in a subclass.
			value = compute_value(inds) ;
			m_hash[ind] = value ;
			// m_hash updated, now return the value.
		}
	}
	return value ;
}

inline uint64_t LookupTable::sub2ind(const std::initializer_list<int> inds)
{
	// 1D offset.
	uint64_t offset = *inds.begin() ;

	// At least 2...calculate the 2D offset.
	if (m_n_dims >= 2) {
		offset += *(inds.begin() + 1)  * m_dims_cp[0] ;
	}
	// and if more than 2D keep adding to the offset.
	if (m_n_dims > 2) {
		// More than 2...calculate.
		for (uint64_t i = 2 ; i < m_n_dims ; ++i) {
			offset += m_dims_cp[i-1] *  *(inds.begin() + i) ;
		}
	}
	// Final check: if the offset is >= than the total number
	// of elements throw an exception!
	if (offset >= m_num_elems) {
		std::ostringstream msg ;
		msg << "Linear index " << offset << " is not within the dimensions of this " ;
		msg << m_n_dims << " dimensional table. " ;
		msg << "The dimensions of this table are: [ " ;
		for (auto const &elem: m_dims) {
			msg << elem << ", " ;
		}
		msg << "]. The requested index was: [ " ;
		for (auto const &elem: inds) {
			msg << elem << ", " ;
		}
		msg << "]." ;
		throw std::runtime_error(msg.str()) ;
	}
	return offset ;
}



inline double LookupTable::compute_value(const std::initializer_list<int>) {
	throw "LookupTable::compute_value() must be implemented in a child class." ;
}

uint64_t LookupTable::count() const {
	// Return the number of items in m_hash.
	uint64_t num_elems ;
	// Make sure no threads are updating the hash table
	// when this is called.
#pragma omp critical (hash_access)
	{
		num_elems = m_hash.size() ;
	}
	return num_elems ;
}

double LookupTable::fill_percent() const {
	return static_cast<double>(this->count()) / m_num_elems * 100.0 ;
}


// Below are the overrides for the operator() method.  This lets
// access to the table elements be used as:  table3d(i,j,k)

inline double LookupTable::operator()(const int i0)  {
	// Call the superclass operator() with the args as an initializer list.  This class's
	// compute_value will be called where needed to compute the desired value.
	return LookupTable::operator()({i0}) ;
}

inline double LookupTable::operator()(const int i0, const int i1)  {
	// Call the superclass operator() with the args as an initializer list.  This class's
	// compute_value will be called where needed to compute the desired value.
	return LookupTable::operator()({i0,i1}) ;
}

inline double LookupTable::operator()(const int i0, const int i1, const int i2)  {
	// Call the superclass operator() with the args as an initializer list.  This class's
	// compute_value will be called where needed to compute the desired value.
	return LookupTable::operator()({i0,i1,i2}) ;
}
inline double LookupTable::operator()(const int i0, const int i1, const int i2, const int i3)  {
	// Call the superclass operator() with the args as an initializer list.  This class's
	// compute_value will be called where needed to compute the desired value.
	return LookupTable::operator()({i0,i1,i2,i3}) ;
}

inline double LookupTable::operator()(const int i0, const int i1, const int i2, const int i3, const int i4)  {
	// Call the superclass operator() with the args as an initializer list.  This class's
	// compute_value will be called where needed to compute the desired value.
	return LookupTable::operator()({i0,i1,i2,i3,i4}) ;
}

