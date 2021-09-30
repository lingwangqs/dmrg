/*
 * LookupTable.cpp
 *
 *  Created on: Jul 13, 2021
 *      Author: bgregor
 */

#include <LookupTable.hpp>
#include <algorithm>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <sstream>

#include "backward.hpp"
using namespace backward;

LookupTable::LookupTable(const double default_value, const std::initializer_list<int> dims) :
        m_default_value(default_value) {
	// Store the dimensions
	m_dims.insert(m_dims.begin(),dims) ; // @suppress("Ambiguous problem")
	m_n_dims = m_dims.size() ;
	// Calculate the cumulative product, this is needed for the sub2ind method.
	m_dims_cp = std::vector<uint64_t>(m_dims.begin(), m_dims.end()) ;
	for (unsigned int i = 1 ; i < m_n_dims ; ++i) {
		m_dims_cp[i] *= m_dims_cp[i-1] ;
	}
	// Get the total possible number of elements to store
	m_num_elems = std::accumulate(m_dims.begin(), m_dims.end(), static_cast<uint64_t>(1), std::multiplies<uint64_t>());
}

inline double LookupTable::operator()(const std::initializer_list<int> inds) const {
	// Look up a value.
	uint64_t ind = sub2ind(inds) ;
	// For use with multiple threads - the "omp critical" forces serialization on
	// the read and possible update to m_hash. see the assign() method for an example.
	double value  ;
//#pragma omp critical (hash_access)
//	{
	// Check to see if it's in the table already.
	// The type of check is that of a const iterator into the map.
	auto check = m_hash.find(ind);
	if (check != m_hash.end()) {
		value = check->second ;  // Found it!
	} else {
                value  = m_default_value ; // Use the default value.
        }
//	} // can't return from inside an omp critical section.
	return value ;
}

inline double& LookupTable::operator ()(const std::initializer_list<int> inds) {
	// Set a value in the hash table.  
	uint64_t ind = sub2ind(inds) ;
        // Returning the value by reference is enough - if this index value
        // doesn't exist it will be created.
        return m_hash[ind] ;
}


inline uint64_t LookupTable::sub2ind(const std::initializer_list<int> inds) const
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
                
// Pretty print a stack trace
StackTrace st; st.load_here(32);

TraceResolver tr; tr.load_stacktrace(st);
for (size_t i = 0; i < st.size(); ++i) {
	ResolvedTrace trace = tr.resolve(st[i]);
	std::cout << "#" << i
		<< " " << trace.object_filename
		<< " " << trace.object_function
		<< " [" << trace.addr << "]"
	<< std::endl << std::flush;
}

                
                
		throw std::runtime_error(msg.str()) ;
	}
	return offset ;
}


uint64_t LookupTable::count() const {
	// Return the number of items in m_hash.
	uint64_t num_elems ;
	// Make sure no threads are updating the hash table
	// when this is called.
//#pragma omp critical (hash_access)
//	{
		num_elems = m_hash.size() ;
//	}
	return num_elems ;
}

double LookupTable::fill_percent() const {
	return static_cast<double>(this->count()) / m_num_elems * 100.0 ;
}

        
double LookupTable::get_default_value() const {
        return m_default_value ;
}
        
uint64_t LookupTable::overflow_size() {
        return m_hash.overflow_size() ;
}


LookupTable_5::LookupTable_5(const double default_value, const int dim0, const int dim1, const int dim2,
	const int dim3, const int dim4) : LookupTable(default_value, {dim0, dim1, dim2, dim3, dim4}) { }

inline double LookupTable_5::operator()(const int i0, const int i1, const int i2, const int i3, const int i4) const {
	// Call the superclass operator() with the args as an initializer list.  This class's
	// compute_value will be called where needed to compute the desired value.
	return LookupTable::operator()({i0,i1,i2,i3,i4}) ;
}

inline double& LookupTable_5::operator ()(const int i0, const int i1, const int i2,
		const int i3, const int i4) {
	return LookupTable::operator()({i0,i1,i2,i3,i4}) ;
}

LookupTable_4::LookupTable_4(const double default_value, const int dim0, const int dim1, const int dim2,
	const int dim3) : LookupTable(default_value, {dim0, dim1, dim2, dim3}) { }

inline double LookupTable_4::operator()(const int i0, const int i1, const int i2, const int i3) const {
	// Call the superclass operator() with the args as an initializer list.  This class's
	// compute_value will be called where needed to compute the desired value.
	return LookupTable::operator()({i0,i1,i2,i3}) ;
}

inline double& LookupTable_4::operator ()(const int i0, const int i1, const int i2,
		const int i3) {
	return LookupTable::operator()({i0,i1,i2,i3}) ;
}

LookupTable_3::LookupTable_3(const double default_value, const int dim0, const int dim1, const int dim2) : LookupTable(default_value, {dim0, dim1, dim2}){}

inline double LookupTable_3::operator ()(const int i0, const int i1,
		const int i2) const {
	return LookupTable::operator()({i0,i1,i2}) ;
}

inline double& LookupTable_3::operator ()(const int i0, const int i1, const int i2) {
	return LookupTable::operator()({i0,i1,i2}) ;
}

LookupTable_1::LookupTable_1(const double default_value, const int dim0) : LookupTable(default_value, {dim0}){}

inline double LookupTable_1::operator ()(const int i0) const {
	return LookupTable::operator()({i0}) ;
}

inline double& LookupTable_1::operator ()(const int i0) {
	return LookupTable::operator()({i0}) ;
}

LookupTable_2::LookupTable_2(const double default_value, const int dim0, const int dim1): LookupTable(default_value, {dim0, dim1}){}

inline double LookupTable_2::operator ()(const int i0, const int i1) const {
	return LookupTable::operator()({i0,i1}) ;
}

inline double& LookupTable_2::operator ()(const int i0, const int i1) {
	return LookupTable::operator()({i0,i1}) ;
}
