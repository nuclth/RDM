/***************************************************

Header to define custom tensor data structures.

***************************************************/

#ifndef MATRIX_DEFINE_H
#define MATRIX_DEFINE_H

#include "boost/multi_array.hpp"
#include "boost/array.hpp"

typedef boost::multi_array<double, 1> one_array;
typedef boost::multi_array<double, 2> two_array;
typedef boost::multi_array<double, 3> three_array;
typedef boost::multi_array<double, 4> four_array;
typedef boost::multi_array<double, 5> five_array;
typedef boost::multi_array<double, 6> six_array;
typedef boost::multi_array<double, 7> seven_array;
typedef boost::multi_array<double, 8> eight_array;

#endif
