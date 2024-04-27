/* *********************************************
Copyright (c) 2013-2024, Cornelis Jan (Jacco) van de Streek
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

#include "FingerCoxJephcoat_functions.h"

#include "Angle.h"
#include "MathsFunctions.h"

#include <stdexcept>
#include <iostream> // For debugging

// ********************************************************************************

double h( const double two_phi, const double two_theta, const double L )
{
    return L * sqrt( ( square( cos( two_phi * degrees2radians ) ) / square( cos( two_theta * degrees2radians ) ) ) - 1.0 );
}

// ********************************************************************************

double W( const double two_phi, const double two_theta, const double two_phi_min, const double two_phi_inflection, const double H, const double S, const double L )
{
    if ( two_phi < two_phi_min )
        return 0.0;
    // h( 2phi ) goes to 0.0 as 2phi goes to 2theta.
    if ( two_phi > two_theta )
        return 0.0;
    if ( two_phi < two_phi_inflection )
        return H + S - h( two_phi, two_theta, L );
    return 2.0 * std::min( H, S );
}

// ********************************************************************************

double D( const double two_phi, const double two_theta, const double two_phi_min, const double two_phi_inflection, const double H, const double S, const double L )
{
    return ( L / ( 2.0 * H * h( two_phi, two_theta, L ) * cos( two_phi * degrees2radians ) ) ) * W( two_phi, two_theta, two_phi_min, two_phi_inflection, H, S, L );
}

// ********************************************************************************

std::vector< double > asymmetric_peak( const double two_theta, const std::vector< double > & two_phi_values, const double FWHM, const double H, const double S, const double L )
{
    if ( two_theta > 90.0 )
        throw std::runtime_error( "asymmetric_peak(): two_theta > 90.0." );
    size_t N( 30 ); // Number of points for the Gauss-Legendre quadrature.
    // N can be made to depend on the distance from two_theta. In the Finger-Cox-Jephcoat paper they
    // recommend N = 30 for the region within six FWHMs from two_theta, N = 10 for everywhere else.
    double two_phi_min = 0.0;
    double cosine_argument = cos( two_theta * degrees2radians ) * sqrt( square( ( H + S ) / L ) + 1.0 );
    if ( cosine_argument < 1.0 )
        two_phi_min = radians2degrees * acos( cosine_argument );
    double two_phi_inflection = 0.0;
    cosine_argument = cos( two_theta * degrees2radians ) * sqrt( square( ( H - S ) / L ) + 1.0 );
    if ( cosine_argument < 1.0 )
        two_phi_inflection = radians2degrees * acos( cosine_argument );
    std::vector< double > result( two_phi_values.size(), 0.0 );
    std::vector< double > x_i;
    std::vector< double > w_i;
    Gauss_Legendre_quadrature( -1.0, 1.0, N, x_i, w_i );
    for ( size_t i( 0 ); i != two_phi_values.size(); ++i )
    {
        if ( two_phi_values[i] > 90.0 )
            throw std::runtime_error( "asymmetric_peak(): two_phi_values[i] > 90.0." );
    //    std::cout << D( two_phi_values[i], two_theta, two_phi_min, two_phi_inflection, H, S, L ) << std::endl;
        double numerator = 0.0;
        double denominator = 0.0;
        for ( size_t j( 0 ); j != N; ++j )
        {
            double delta_n = ( two_theta + two_phi_min ) / 2.0 + ( two_theta - two_phi_min ) * x_i[j] / 2.0;
//std::cout << "delta_n " << delta_n << std::endl;

//std::cout << "w " << w_i[j] << std::endl;
//std::cout << "W " << W( delta_n, two_theta, two_phi_min, two_phi_inflection, H, S, L ) << std::endl;
//std::cout << "h " << h( delta_n, two_theta, L ) << std::endl;
//std::cout << "cos " << cos( delta_n * degrees2radians ) << std::endl;
//std::cout << " D " << D( delta_n, two_theta, two_phi_min, two_phi_inflection, H, S, L ) << std::endl;

            double term = ( W( delta_n, two_theta, two_phi_min, two_phi_inflection, H, S, L ) ) / h( delta_n, two_theta, L );
            term *= w_i[j] / cos( delta_n * degrees2radians );
            numerator   += term * pseudo_Voigt( two_phi_values[i] - delta_n, FWHM, 0.9 );
            denominator += term;
        }
        result[i] = numerator / denominator;
    }
//std::cout << "two_phi_min = " << two_phi_min << std::endl;
//std::cout << "two_phi_inflection = " << two_phi_inflection << std::endl;
    return result;
}

// ********************************************************************************

