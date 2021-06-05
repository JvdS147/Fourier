#ifndef MATHSFUNCTIONS_H
#define MATHSFUNCTIONS_H

/* *********************************************
Copyright (c) 2013-2021, Cornelis Jan (Jacco) van de Streek
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of my employers nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL CORNELIS JAN VAN DE STREEK BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
********************************************* */

class Vector3D;
class Angle;

#include <cstddef> // For definition of size_t
#include <vector>

// Adds a list of doubles trying to avoid adding very small to very large numbers
double add_doubles( const std::vector< double > & values );

// Adds a list of doubles trying to avoid adding very small to very large numbers
double add_absolute_doubles( const std::vector< double > & values );

// Adds a list of doubles trying to avoid adding very small to very large numbers
double add_squared_doubles( const std::vector< double > & values );

double calculate_average( const std::vector< double > & values );
double calculate_minimum( const std::vector< double > & values );
double calculate_maximum( const std::vector< double > & values );
size_t calculate_minimum( const std::vector< size_t > & values );
size_t calculate_maximum( const std::vector< size_t > & values );

int greatest_common_divisor( const int lhs, const int rhs );

int round_to_int( const double x );

size_t round_to_size_t( const double x );

// ********************************************************************************

double absolute_relative_difference( const double lhs, const double rhs );

// It is more efficient to calculate sine and cosine of the same angle simultaneously
// The algorithm that is used is an APPROXIMATION
void sincos( Angle angle, double & sine, double & cosine );

// Calculates sqrt( x^2 + y^2 ) trying to prevent underflow or overflow
double hypothenuse( const double x, const double y );

inline double square( const double x )
{
    return ( x * x );
}

// returns 0 if x = 0.0
inline int sign( const double x )
{
    return ( 0.0 < x ) - ( x < 0.0 );
}

inline bool is_even( const int i )
{
    return ( i % 2 ) == 0;
}

inline bool is_odd( const int i )
{
    return ! is_even( i );
}

Vector3D cylindrical2Cartesian( const double r, Angle phi, const double z );

// Returns a random number between 0.0 and 1.0, inclusive.
double uniform_distribution_1();

// Returns a random number between 0.0 and 1.0, exclusive.
double uniform_distribution_2();

double normal_distribution( const double mean = 0.0, const double sigma = 1.0 );

size_t Poisson_distribution( const double mean );

// For testing
size_t Poisson_distribution_exact( const double mean );
size_t Poisson_distribution_Gaussian_0( const double mean );
size_t Poisson_distribution_Gaussian_1( const double mean );
size_t Poisson_distribution_Gaussian_2( const double mean );

double calculate_sample_RMSD( const std::vector< double > & lhs, const std::vector< double > & rhs );

double ln( const double x );

#endif // MATHSFUNCTIONS_H

