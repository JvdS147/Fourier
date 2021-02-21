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

#include "MathFunctions.h"
#include "Angle.h" // This is bound to give circular references later on
#include "MathConstants.h"
#include "Sort.h"
#include "Vector3D.h" // This is bound to give circular references later on

#include <cmath>
#include <stdexcept>
#include <cstdlib>
#include <iostream> // For debugging

// ********************************************************************************

// @@ Not sophisticated enough, but a start.
// Two problems: should sort by magnitude (absolute value) and should use a better algorithm (adding pairwise, for example)
double add_doubles( const std::vector< double > & values )
{
    double result( 0.0 );
    std::vector< size_t > sorted_map = sort( values );
    for ( size_t i( 0 ); i != sorted_map.size(); ++i )
        result += values[ sorted_map[i] ];
    return result;
}

// ********************************************************************************

// @@ Not sophisticated enough, but a start.
// Should use a better algorithm (adding pairwise, for example)
double add_absolute_doubles( const std::vector< double > & values )
{
    double result( 0.0 );
    std::vector< double > absolute_values;
    absolute_values.reserve( values.size() );
    for ( size_t i( 0 ); i != values.size(); ++i )
        absolute_values.push_back( std::abs( values[ i ] ) );
    std::vector< size_t > sorted_map = sort( absolute_values );
    for ( size_t i( 0 ); i != sorted_map.size(); ++i )
        result += std::abs( absolute_values[ sorted_map[i] ] );
    return result;
}

// ********************************************************************************

// @@ Not sophisticated enough, but a start.
// Should use a better algorithm (adding pairwise, for example)
double add_squared_doubles( const std::vector< double > & values )
{
    double result( 0.0 );
    std::vector< double > squared_values;
    squared_values.reserve( values.size() );
    for ( size_t i( 0 ); i != values.size(); ++i )
        squared_values.push_back( std::abs( values[ i ] ) );
    std::vector< size_t > sorted_map = sort( squared_values );
    for ( size_t i( 0 ); i != sorted_map.size(); ++i )
        result += std::abs( squared_values[ sorted_map[i] ] );
    return result;
}

// ********************************************************************************

double calculate_average( const std::vector< double > & values )
{
    if ( values.empty() )
        throw std::runtime_error( "calculate_average(): no values." );
    return add_doubles( values ) / values.size();
}

// ********************************************************************************

double calculate_minimum( const std::vector< double > & values )
{
    if ( values.empty() )
        throw std::runtime_error( "calculate_minimum(): empty" );
    double result = values[0];
    for ( size_t i( 1 ); i != values.size(); ++i )
        result = std::min( result, values[i] );
    return result;
}

// ********************************************************************************

double calculate_maximum( const std::vector< double > & values )
{
    if ( values.empty() )
        throw std::runtime_error( "calculate_maximum(): empty" );
    double result = values[0];
    for ( size_t i( 1 ); i != values.size(); ++i )
        result = std::max( result, values[i] );
    return result;
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

void sincos( Angle angle, double & sine, double & cosine )
{
//    sine = sin( x );

// Sine
    double x = angle.value_in_radians();
    while ( x > CONSTANT_PI )
        x -= 2.0*CONSTANT_PI;
    while ( x < -CONSTANT_PI )
        x += 2.0*CONSTANT_PI;
    double B = 4.0/CONSTANT_PI;
    double C = -4.0/(CONSTANT_PI*CONSTANT_PI);
    double y = B * x + C * x * fabs(x);
    double P = 0.225;
    sine = P * (y * fabs(y) - y) + y; // = Q * y + P * y * fabs(y), Q = 0.775 ( P + Q = 1.0 )

//    cosine = cos( x );
    x += CONSTANT_PI/2.0; // cos(x) = sin(x + pi/2)
    if ( x > CONSTANT_PI ) // Original x > pi/2
        x -= 2.0*CONSTANT_PI; // Wrap: cos(x) = cos(x - 2 pi)
    y = B * x + C * x * fabs(x);
    cosine = P * (y * fabs(y) - y) + y; // = Q * y + P * y * fabs(y), Q = 0.775 ( P + Q = 1.0 )
}

// ********************************************************************************

// hypothenuse
double hypothenuse( const double x, const double y )
{
    double ax = fabs( x );
    double ay = fabs( y );
    double t = std::min( ax, ay );
    ax = std::max( ax, ay );
    if ( t == 0.0 )
        return ax;
    t /= ax;
    return ax * sqrt( 1.0 + square( t ) );
}

// ********************************************************************************

Vector3D cylindrical2Cartesian( const double r, Angle phi, const double z )
{
    return Vector3D( r * phi.cosine(), r * phi.sine(), z );
}

// ********************************************************************************

double uniform_distribution_1()
{
    return static_cast<double>( rand() ) / static_cast<double>( RAND_MAX );
}

// ********************************************************************************

double uniform_distribution_2()
{
    int value;
    do
    {
        value = rand();
    }
    while ( ( value == 0 ) || ( value == RAND_MAX ) );
//    std::cout << "value = " << value << std::endl;
    return static_cast<double>( value ) / static_cast<double>( RAND_MAX );
//    return ( static_cast<double>( rand() ) + 1.0 ) / ( static_cast<double>( RAND_MAX ) + 2.0 );
}

// ********************************************************************************

double normal_distribution( const double mean, const double sigma )
{
    double u = uniform_distribution_2();
    double v = uniform_distribution_2();
    // Box-Muller method, we discard the second random number
    double x = sqrt( -2.0 * log( u ) ) * cos( 2.0 * CONSTANT_PI * v );
//    y = sqrt( -2.0 * log( u ) ) * sin( 2.0 * CONSTANT_PI * v );
    if ( x > 10000.0 )
        std::cout << "x = " << x << ", u = " << u << ", v = " << v << std::endl;
    return mean + sigma*x;
}

// ********************************************************************************

size_t Poisson_distribution( const double mean )
{
    if ( mean < 0.0 )
        throw std::runtime_error( "poisson_distribution(): mean cannot be negative." );
    if ( mean > 100.0 )
        return round_to_int( normal_distribution( mean, sqrt( mean ) ) );
    // The following is exact but numerically unstable for values greater than about 800
    return Poisson_distribution_exact( mean );
}

// ********************************************************************************

size_t Poisson_distribution_exact( const double mean )
{
    if ( mean < 0.0 )
        throw std::runtime_error( "poisson_distribution(): mean cannot be negative." );
    // The following is exact but numerically unstable for values greater than about 800
    // because exp( -mean ) sooner or later becomes 0.0
    double L = exp( -mean );
    size_t k = 0;
    double p = 1.0;
    do
    {
        ++k;
        double u = uniform_distribution_2();
        p *= u;
    }
    while ( p > L );
    return k - 1;
}

// ********************************************************************************

size_t Poisson_distribution_Gaussian_0( const double mean )
{
    if ( mean < 0.0 )
        throw std::runtime_error( "poisson_distribution(): mean cannot be negative." );
    return round_to_size_t( std::abs( normal_distribution( mean, sqrt( mean ) ) ) );
}

// ********************************************************************************

size_t Poisson_distribution_Gaussian_1( const double mean )
{
    if ( mean < 0.0 )
        throw std::runtime_error( "poisson_distribution(): mean cannot be negative." );
    return round_to_size_t( std::abs( normal_distribution( mean-0.5, sqrt( mean ) ) ) );
}

// ********************************************************************************

size_t Poisson_distribution_Gaussian_2( const double mean )
{
    if ( mean < 0.0 )
        throw std::runtime_error( "poisson_distribution(): mean cannot be negative." );
    return round_to_size_t( std::abs( normal_distribution( mean+0.5, sqrt( mean ) ) ) );
}

// ********************************************************************************

double calculate_sample_RMSD( const std::vector< double > & lhs, const std::vector< double > & rhs )
{
    // @@ Poor design, should use data structure that forces the two to be the same
    if ( lhs.size() != rhs.size() ) 
        throw std::runtime_error( "calculate_RMSD(): left-hand side and right-hand side should contain the same number of values." );
    if ( lhs.empty() ) 
        throw std::runtime_error( "calculate_RMSD(): lists are empty." );
    if ( lhs.size() == 1 ) 
        throw std::runtime_error( "calculate_RMSD(): lists must contain at least two values." );
    double average_lhs = calculate_average( lhs );
    double average_rhs = calculate_average( rhs );
    double result( 0.0 );
    for ( size_t i( 0 ); i != lhs.size(); ++i )
        result += square( ( lhs[i] - average_lhs ) - ( rhs[i] - average_rhs ) );
    result /= ( lhs.size() - 1.0 );
    return std::sqrt( result );
}

// ********************************************************************************

