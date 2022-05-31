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
    if ( ( two_phi + 0.001 ) > two_theta )
        return 0.0;
    if ( two_phi < two_phi_inflection )
        return H + S - h( two_phi, two_theta, L );
    return 2.0 * std::min( H, S );
}

// ********************************************************************************

double R( const double x, const double two_theta, const double FWHM, const double eta )
{
    return pseudo_Voigt( x - two_theta, FWHM, eta );
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
    double two_phi_min = 0.0;
    double cosine_argument = cos( two_theta * degrees2radians ) * sqrt( square( ( H + S ) / L ) + 1.0 );
//std::cout << "cosine_argument = " << cosine_argument << std::endl;
    if ( cosine_argument < -0.999 )
        two_phi_min = 180.0;
    else if ( cosine_argument > 0.999 )
        two_phi_min = 0.0;
    else
        two_phi_min = radians2degrees * acos( cosine_argument );
    double two_phi_inflection = 0.0;
    cosine_argument = cos( two_theta * degrees2radians ) * sqrt( square( ( H - S ) / L ) + 1.0 );
//std::cout << "cosine_argument = " << cosine_argument << std::endl;
    if ( cosine_argument < -0.999 )
        two_phi_inflection = 180.0;
    else if ( cosine_argument > 0.999 )
        two_phi_inflection = 0.0;
    else
        two_phi_inflection = radians2degrees * acos( cosine_argument );
//std::cout << "two_phi_min = " << two_phi_min << std::endl;
//std::cout << "two_phi_inflection = " << two_phi_inflection << std::endl;
    std::vector< double > result( two_phi_values.size(), 0.0 );
    std::vector< double > x_i;
    std::vector< double > w_i;
    Gauss_Legendre_quadrature( -1.0, 1.0, N, x_i, w_i );
      //  std::cout << "####################    DDDDD    ######################" << std::endl;
    
    for ( size_t i( 0 ); i != two_phi_values.size(); ++i )
    {
        if ( two_phi_values[i] > 90.0 )
            throw std::runtime_error( "asymmetric_peak(): two_phi_values[i] > 90.0." );
    //    std::cout << D( two_phi_values[i], two_theta, two_phi_min, two_phi_inflection, H, S, L ) << std::endl;
        if ( two_phi_values[i] < two_phi_min )
            continue;
        // h( 2phi ) goes to 0.0 as 2phi goes to 2theta.
        if ( ( two_phi_values[i] + 0.001 ) > two_theta )
            continue;
        double numerator = 0.0;
        double denominator = 0.0;
        for ( size_t j( 0 ); j != N; ++j )
        {
            double delta_n = ( two_theta + two_phi_min ) / 2.0 + ( two_theta - two_phi_min ) * x_i[j] / 2.0;
//            numerator   += ( w_i[j] * W( delta_n, two_theta, two_phi_min, two_phi_inflection, H, S, L ) * R( two_phi_values[i] - delta_n, two_theta, FWHM, 0.9 ) ) / ( h( delta_n, two_theta, L ) * cos( delta_n * degrees2radians ) );
//            denominator += ( w_i[j] * W( delta_n, two_theta, two_phi_min, two_phi_inflection, H, S, L )                                                          ) / ( h( delta_n, two_theta, L ) * cos( delta_n * degrees2radians ) );
            double term = ( w_i[j] * W( delta_n, two_theta, two_phi_min, two_phi_inflection, H, S, L ) ) / ( h( delta_n, two_theta, L ) * cos( delta_n * degrees2radians ) );
            numerator   += term * R( two_phi_values[i] - delta_n, two_theta, FWHM, 0.9 );
            denominator += term;
        }
        result[i] = numerator / denominator;
    }
  //      std::cout << "####################    DDDDD    ######################" << std::endl;
    return result;
}

// ********************************************************************************

