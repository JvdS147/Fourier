#ifndef BASICMATHSFUNCTIONS_H
#define BASICMATHSFUNCTIONS_H

/* *********************************************
Copyright (c) 2013-2023, Cornelis Jan (Jacco) van de Streek
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

#include <cstddef> // For definition of size_t

const double TOLERANCE = 0.000001;

// To calculate the average of four values:
// double average = average( value_1, value_2 );
// average = average( value_3, average, 2.0 );
// average = average( value_4, average, 3.0 );
// Even this works:
//    double prev_estimate; // No need to initialise...
//    double next_estimate; // No need to initialise...
//    size_t iStep( 0 );
//    while ( iStep < 1000000 )
//    {
//        prev_estimate = next_estimate;
//        double current_value = some_function();
//        next_estimate = average( current_value, prev_estimate, iStep );
//        ++iStep;
//    }
double average( const double lhs, const double rhs, const double weight = 1.0 );

int greatest_common_divisor( const int lhs, const int rhs );

int round_to_int( const double x );

size_t round_to_size_t( const double x );

inline double absolute( const double x )
{
    return ( x < 0.0 ) ? -x : x;
}

inline size_t absolute( const int x )
{
    return ( x < 0 ) ? -x : x;
}

inline bool nearly_equal( const double lhs, const double rhs, const double tolerance = TOLERANCE )
{
    return ( absolute( rhs - lhs ) < tolerance );
}

bool triquality( const double x1, const double x2, const double x3, const double tolerance = TOLERANCE );

inline bool nearly_zero( const double lhs, const double tolerance = TOLERANCE )
{
    return ( absolute( lhs ) < tolerance );
}

// This has to be here because it uses functions defined in Utilities.h and in MathsFunctions.h
// and otherwise there would be circular references. Well, that is why we started the current file.
bool nearly_integer( const double value, const double tolerance = TOLERANCE );

double absolute_relative_difference( const double lhs, const double rhs );

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

#endif // BASICMATHSFUNCTIONS_H

