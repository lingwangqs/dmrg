/*
 *  Fast lookup table example.
 *
 *  Created on: Jul 20, 2021
 *      Author: bgregor
 */


#include "LUT_5D_Example.h"

#include <iostream>
using namespace std ;

int main(int argc, char **argv) {

	// Make a 5D lookup table.
	double phase = 0.2 ;
	LUT_5D_Example table(8,3,4,5,6, phase) ;

	// Turn on verbose table output.
	table.m_verbose_computation = true ;

	// Fetch a value from the table. This should print a computation message.
	cout << "Fetching a value from the table. This fetch should print a computation message." << endl ;
	double value = table(1,1,2,3,4) ;
	cout << "table(1,1,2,3,4) = " << value << endl ;

	cout << "Fetching the same value from the table. This fetch should NOT print a computation message." << endl ;
	cout << "table(1,1,2,3,4) = " << value << endl ;

	// Turn off verbosity.
	table.m_verbose_computation = false ;

	// Fetch 240 values.  Do this twice.  After the
	// count should be 241 as we've looked up 241 unique values.
	value = 0 ;
	// Fetching again just retrieves values, it doesn't recompute them.
#pragma omp parallel for shared(table) reduction(+:value)
	for (int i = 0 ; i < 8 ; ++i) {
		for (int j = 0 ; j < 2 ; ++j) {
			for (int k = 0 ; k < 4 ; ++k) {
				for (int l = 0 ; l < 5 ; ++l) {
					for (int m = 0 ; m < 6 ; ++m) {
						value += table(i,j,k,l,m) ;
					}
					value += table(0,0,0,0,0) ;
				}
			}
		}
	}
	cout << "VALUE " << value << endl ;
	value = 0 ;
	// Fetching again just retrieves values, it doesn't recompute them.
#pragma omp parallel for shared(table) reduction(+:value)
	for (int i = 0 ; i < 8 ; ++i) {
		for (int j = 0 ; j < 2 ; ++j) {
			for (int k = 0 ; k < 4 ; ++k) {
				for (int l = 0 ; l < 5 ; ++l) {
					for (int m = 0 ; m < 6 ; ++m) {
						value += table(i,j,k,l,m) ;
					}
					value += table(0,0,0,0,0) ;
				}
			}
		}
	}

	cout << "VALUE " << value << endl ;

	value = 0 ;
#pragma omp parallel for shared(table) reduction(+:value)
	for (int i = 0 ; i < 100000 ; ++i)
		value += table(0,0,0,0,0) ;

	// Print the number of elements in the table.
	cout << "The table has stored " << table.count() << " values." << endl ;
	// Print the percentage of the table that has been used.
	cout << "This is " << table.fill_percent() << "% of the total possible storage." << endl ;



	return 0;
}
