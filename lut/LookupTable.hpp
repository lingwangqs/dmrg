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


// This implements a sparse N-dimensional lookup table.
// It stores double precision values. It could be templated
// to store any value, of course. Values placed in the table
// are stored in a hash table, and any index that's requested
// returns the default value.  The default value is set when
// the table object is constructed.

// Updating the table is not thread safe, although there is 
// commented code to implement OpenMP locks that make updating
// the table thread-safe. It is commented out to avoid the minor
// performance hit when using locks as the table updates in this
// program are done in a single thread.


class LookupTable {
public:
        LookupTable() = default;

        // This takes any number of ints which are the dimensions of the
        // matrix to construct.  The initializer list handles any number
        // of arguments.
        LookupTable(const double default_value, const std::initializer_list<int> dims) ;

        virtual ~LookupTable() = default ;

        // This is the "real" operator() method. Its arguments are in the form
        // of an initializer list. The overrides below just call this.
        virtual double operator()(const std::initializer_list<int> inds) const;

        // This returns the location of the stored value in the map, so that
        // it can be assigned, i.e. tbl_obj(2,3,4)=32.0 ;
        virtual double& operator() (const std::initializer_list<int> inds);

        // How many values are actually stored in the table?
        uint64_t count() const ;

        // What % of the total potential storage is used in the table?
        double fill_percent() const ;

        //  getter for the default value.
        double get_default_value() const ;
        
protected:
        // What's the default value of an unfilled index?
        double m_default_value  ;
        
        std::vector<int> m_dims ; // Storage for the matrix dimensions
        // cumulative product of the dimensions. This is needed for the
        // sub2ind calculation, it's pre-calculated at initialization time
        // for performance.
        std::vector<uint64_t> m_dims_cp ;
        uint64_t m_num_elems = 0 ; // Total number of elements that could be stored.
        unsigned int m_n_dims = 0 ; // Number of dimensions

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
    uint64_t sub2ind(const std::initializer_list<int> inds) const ;
};

// Specific implementations for 5D, 4D, 3D, 2D, and 1D tables.
// Just to be complete...


class LookupTable_1 : public LookupTable {
public:
        LookupTable_1() {}
        LookupTable_1(const double default_value, const int dim0) ;

        virtual double operator() (const int i0) const ;
        virtual double& operator() (const int i0) ;
};

class LookupTable_2 : public LookupTable {
public:
        LookupTable_2() {}
        LookupTable_2(const double default_value, const int dim0, const int dim1) ;

        virtual double operator() (const int i0, const int i1) const ;
        virtual double& operator() (const int i0, const int i1) ;
};

class LookupTable_3 : public LookupTable {
public:
        LookupTable_3() {}
        LookupTable_3(const double default_value, const int dim0, const int dim1, const int dim2) ;

        virtual double operator() (const int i0, const int i1, const int i2) const ;
        virtual double& operator() (const int i0, const int i1, const int i2) ;
};

class LookupTable_4 : public LookupTable {
public:
        LookupTable_4() {}
        LookupTable_4(const double default_value, const int dim0, const int dim1, const int dim2, const int dim3) ;

        virtual double operator() (const int i0, const int i1, const int i2, const int i3) const ;
        virtual double& operator() (const int i0, const int i1, const int i2, const int i3) ;

};

class LookupTable_5 : public LookupTable {
public:
        LookupTable_5() {}
        LookupTable_5(const double default_value, const int dim0, const int dim1, const int dim2, const int dim3, const int dim4) ;

        virtual double operator() (const int i0, const int i1, const int i2, const int i3, const int i4) const ;
        virtual double& operator() (const int i0, const int i1, const int i2, const int i3, const int i4)  ;

};

#endif /* LOOKUPTABLE_H_ */
