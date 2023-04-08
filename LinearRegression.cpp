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

#include "LinearRegression.h"

#include "BasicMathsFunctions.h"

// ********************************************************************************

void linear_regression( const std::vector< double > & x, const std::vector< double > & y, const std::vector< double > & s, double & a, double & b )
{
    size_t N = x.size();
    if ( ( y.size() != N ) || ( s.size() != N ) )
        throw std::runtime_error( "linear_regression(): Error: sizes must be equal." );
    if ( N < 2 )
        throw std::runtime_error( "linear_regression(): Error: need at least two points." );
    double S   = 0.0;
    double Sx  = 0.0;
    double Sy  = 0.0;
    double Sxx = 0.0;
    double Sxy = 0.0;
    for ( size_t i( 0 ); i != N; ++i )
    {
        if ( nearly_zero( s[i] ) )
            throw std::runtime_error( "linear_regression(): Error: weight close to 0.0." );
        double ss = 1.0 / square( s[i] );
        S   += ss;
        Sx  += x[i] * ss;
        Sy  += y[i] * ss;
        Sxx += square( x[i] ) * ss;
        Sxy += x[i] * y[i] * ss;
    }
    double Delta = S * Sxx - square( Sx );
    a = ( Sxx*Sy - Sx*Sxy ) / Delta;
    b = ( S*Sxy - Sx*Sy ) / Delta;
}

// ********************************************************************************

void linear_regression( const std::vector< double > & x, const std::vector< double > & y, double & a, double & b )
{
    size_t N = x.size();
    if ( y.size() != N )
        throw std::runtime_error( "linear_regression(): Error: sizes must be equal." );
    if ( N < 2 )
        throw std::runtime_error( "linear_regression(): Error: need at least two points." );
    double Sx  = 0.0;
    double Sy  = 0.0;
    double Sxx = 0.0;
    double Sxy = 0.0;
    for ( size_t i( 0 ); i != N; ++i )
    {
        Sx  += x[i];
        Sy  += y[i];
        Sxx += square( x[i] );
        Sxy += x[i] * y[i];
    }
    double Delta = N * Sxx - square( Sx );
    a = ( Sxx*Sy - Sx*Sxy ) / Delta;
    b = ( N*Sxy - Sx*Sy ) / Delta;
}

// ********************************************************************************

void fit_exponential( const std::vector< double > & x, const std::vector< double > & y, double & a, double & b )
{
    std::vector< double > lny;
    lny.reserve( y.size() );
    for ( size_t i( 0 ); i != y.size(); ++i )
        lny.push_back( ln( y[i] ) );
    linear_regression( x, lny, a, b );
    a = exp( a );
}

// ********************************************************************************

LinearFunction linear_regression( const std::vector< double > & x, const std::vector< double > & y, const std::vector< double > & s )
{
    double a;
    double b;
    linear_regression( x, y, s, a, b );
    return LinearFunction( a, b );
}

// ********************************************************************************

LinearFunction linear_regression( const std::vector< double > & x, const std::vector< double > & y )
{
    double a;
    double b;
    linear_regression( x, y, a, b );
    return LinearFunction( a, b );
}

// ********************************************************************************

ExponentialFunction fit_exponential( const std::vector< double > & x, const std::vector< double > & y )
{
    double a;
    double b;
    fit_exponential( x, y, a, b );
    return ExponentialFunction( a, b );
}

// ********************************************************************************

