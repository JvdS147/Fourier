/* *********************************************
Copyright (c) 2013-2022, Cornelis Jan (Jacco) van de Streek
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

#include "BasicMathsFunctions.h"

#include <stdexcept>
#include <iostream> // For debugging

// ********************************************************************************

double average( const double lhs, const double rhs, const double weight )
{
    if ( ! nearly_integer( weight ) )
        std::cout << "average(): warning: weight is not an integer, is that intended?" << std::endl;
    return ( lhs + weight * rhs ) / ( 1.0 + weight );
}

// ********************************************************************************

int greatest_common_divisor( const int lhs, const int rhs )
{
    return ( ( rhs == 0 ) ? lhs : greatest_common_divisor( rhs, lhs % rhs ) );
}

// ********************************************************************************

int round_to_int( const double x )
{
    return ( x < 0 ) ? static_cast<int>( x - 0.5 ) : static_cast<int>( x + 0.5 );
}

// ********************************************************************************

size_t round_to_size_t( const double x )
{
    if ( x < -0.5 )
        throw std::runtime_error( "round_to_size_t(): value is negative." );
    return static_cast<size_t>( x + 0.5 );
}

// ********************************************************************************

bool triquality( const double x1, const double x2, const double x3, const double tolerance )
{
    double average = ( x1 + x2 + x3 ) / 3.0;
    if ( ! nearly_equal( x1, average, tolerance ) )
        return false;
    if ( ! nearly_equal( x2, average, tolerance ) )
        return false;
    if ( ! nearly_equal( x3, average, tolerance ) )
        return false;
    return true;
}

// ********************************************************************************

bool nearly_integer( const double value, const double tolerance )
{
    return nearly_equal( value, round_to_int( value ), tolerance );
}

// ********************************************************************************

double absolute_relative_difference( const double lhs, const double rhs )
{
    return absolute( lhs - rhs ) / ( 0.5 * ( lhs + rhs ) );
}

// ********************************************************************************

