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

#include "MathsFunctions.h"
#include "Angle.h"
#include "Sort.h"

#include <stdexcept>
#include <cstdlib>
#include <iostream> // For debugging

// ********************************************************************************

double bisection( const Function f, const double target_y_value, const double initial_x_value, const double tolerance )
{
    double increment( 0.01 );
    double left_bracket  = initial_x_value - increment;
    double right_bracket = initial_x_value + increment;
    while ( ! ( sign( f( left_bracket ) - target_y_value ) * sign( f( right_bracket ) - target_y_value ) < 0 ) )
    {
        left_bracket  -= increment;
        right_bracket += increment;
        //increment *= 2.0;
    }
    double result = initial_x_value;
    while ( std::abs( f( result ) - target_y_value ) > tolerance )
    {
        if ( sign( f( result ) - target_y_value ) * sign( f( right_bracket ) - target_y_value ) < 0 )
            left_bracket = result;
        else
            right_bracket = result;
        result = ( right_bracket + left_bracket ) / 2.0;
    }
    return result;
}

// ********************************************************************************

double integral( const Function f, const double start, const double end, const double step_size )
{
    double result( 0.0 );
    size_t nintervals = round_to_int( ( end - start ) / step_size );
    if ( nintervals == 0 )
        throw std::runtime_error( "integral(): step_size too big." );
    double interval = ( end - start ) / nintervals;
    double lhs = f( start );
    for ( size_t i( 0 ); i != nintervals; ++i )
    {
        double rhs = f( start + ( i + 1 ) * interval );
        result += ( ( lhs + rhs ) / 2.0 ) * interval;
        lhs = rhs;
    }
    return result;
}

// ********************************************************************************

double integral_Monte_Carlo( const Function f, const double start, const double end, const size_t npoints )
{
    double result( 0.0 );
    if ( npoints == 0 )
        throw std::runtime_error( "integral_Monte_Carlo(): npoints == 0." );
    double interval = ( end - start ) / npoints;
    for ( size_t i( 0 ); i != npoints; ++i )
    {
        result += f( start + uniform_distribution_1() * ( end - start ) ) * interval;
    }
    return result;
}

// ********************************************************************************

double Legendre_polynomial( const size_t order, const double cosine_theta )
{
    switch ( order )
    {
        case 0 : return 1.0;
        case 1 : return cosine_theta;
        case 2 : return (   3.0 * square( cosine_theta ) - 1.0 ) / 2.0;
        case 3 : return (   5.0 * std::pow( cosine_theta, 3 ) -   3.0 * cosine_theta ) / 2.0;
        case 4 : return (  35.0 * std::pow( cosine_theta, 4 ) -  30.0 * square( cosine_theta ) + 3.0 ) / 8.0;
        case 5 : return (  63.0 * std::pow( cosine_theta, 5 ) -  70.0 * std::pow( cosine_theta, 3 ) + 15.0 * cosine_theta ) / 8.0;
        case 6 : return ( 231.0 * std::pow( cosine_theta, 6 ) - 315.0 * std::pow( cosine_theta, 4 ) + 105.0 * square( cosine_theta ) - 5.0 ) / 16.0;
        default : throw std::runtime_error( "Legendre_polynomial(): order not yet implemented." );
    }
}

// ********************************************************************************

// Must return size_t, because we use recursion
// Simplistic algorithm, will overflow very quickly.
size_t factorial( const size_t n )
{
    if ( n == 0 )
        return 1;
    return n * factorial( n-1 );
}

// ********************************************************************************

// @@ Not sophisticated enough, but a start.
// Two problems: should sort by magnitude (absolute value) and should use a better algorithm (adding pairwise, for example)
double add_doubles( const std::vector< double > & values )
{
    double result( 0.0 );
    Mapping sorted_map = sort( values );
    for ( size_t i( 0 ); i != sorted_map.size(); ++i )
        result += values[ sorted_map[ i ] ];
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
    Mapping sorted_map = sort( absolute_values );
    for ( size_t i( 0 ); i != sorted_map.size(); ++i )
        result += std::abs( absolute_values[ sorted_map[ i ] ] );
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
    Mapping sorted_map = sort( squared_values );
    for ( size_t i( 0 ); i != sorted_map.size(); ++i )
        result += std::abs( squared_values[ sorted_map[ i ] ] );
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
        throw std::runtime_error( "calculate_minimum(double): empty" );
    double result = values[0];
    for ( size_t i( 1 ); i != values.size(); ++i )
        result = std::min( result, values[i] );
    return result;
}

// ********************************************************************************

double calculate_maximum( const std::vector< double > & values )
{
    if ( values.empty() )
        throw std::runtime_error( "calculate_maximum(double): empty" );
    double result = values[0];
    for ( size_t i( 1 ); i != values.size(); ++i )
        result = std::max( result, values[i] );
    return result;
}

// ********************************************************************************

size_t calculate_minimum( const std::vector< size_t > & values )
{
    if ( values.empty() )
        throw std::runtime_error( "calculate_minimum(size_t): empty" );
    size_t result = values[0];
    for ( size_t i( 1 ); i != values.size(); ++i )
        result = std::min( result, values[i] );
    return result;
}

// ********************************************************************************

size_t calculate_maximum( const std::vector< size_t > & values )
{
    if ( values.empty() )
        throw std::runtime_error( "calculate_maximum(size_t): empty" );
    size_t result = values[0];
    for ( size_t i( 1 ); i != values.size(); ++i )
        result = std::max( result, values[i] );
    return result;
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

// Centered around 0.0, area normalised to 1.0.
double Lorentzian( const double x, const double FWHM )
{
    return (1.0/CONSTANT_PI) * ( ( 0.5*FWHM ) / ( square( x ) + square( 0.5*FWHM ) ) );
}

// ********************************************************************************

// Centered around 0.0, area normalised to 1.0.
double Gaussian( const double x, const double FWHM )
{
    double sigma = FWHM / 2.35482; // 2.35482 = 2*sqrt(2*ln(2))
    return exp( -square(x) / (2.0*square( sigma )) ) / ( sigma*sqrt(2.0*CONSTANT_PI) );
}

// ********************************************************************************

// Centered around 0.0, area normalised to 1.0.
// Needs: FWHM (in degrees 2theta), eta, 2theta w.r.t. 0.0
// eta should probably be 0.68 for the FWHM of the pseudo-Voigt to be the same as the
// FWHM of the individual Lorentzian and Gaussian.
double pseudo_Voigt( const double x, const double FWHM )
{
    double eta( 0.9 );
// For a *full* Voigt, when the FWHM is set to the same value for the Gaussian and the Lorentzian part, the FWHM of the
// resulting full Voigt is also that same FWHM. This is no longer true for our pseudo-Voigt, so we
// must calculate the correct FWHM.
    return eta * Lorentzian( x, FWHM ) + (1.0-eta) * Gaussian( x, FWHM );
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

double ln( const double x )
{
    if ( nearly_zero( x ) )
        return 1.0;
    return log( x );
}

// ********************************************************************************

