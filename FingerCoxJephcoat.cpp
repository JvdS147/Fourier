/* *********************************************
Copyright (c) 2013-2025, Cornelis Jan (Jacco) van de Streek
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

#include "FingerCoxJephcoat.h"

#include "Angle.h"
#include "BasicMathsFunctions.h"
#include "MathsFunctions.h"

#include <iostream> // for debugging
#include <stdexcept>

// ********************************************************************************

FingerCoxJephcoat::FingerCoxJephcoat( const double H, const double S, const double L ):eta_(0.9),N_(250)
{
    set_HSL( H, S, L );
    Gauss_Legendre_quadrature( -1.0, 1.0, N_, x_i_, w_i_ );
}

// ********************************************************************************

FingerCoxJephcoat::FingerCoxJephcoat( const double A, const double B ):eta_(0.9),N_(250)
{
    set_A( A );
    set_B( B );
    Gauss_Legendre_quadrature( -1.0, 1.0, N_, x_i_, w_i_ );
}

// ********************************************************************************

FingerCoxJephcoat FingerCoxJephcoat::from_DASH( const double HPSL, const double HMSL )
{
    return FingerCoxJephcoat( 0.5 * ( HPSL + HMSL ), 0.5 * ( HPSL - HMSL ) );
}

// ********************************************************************************

void FingerCoxJephcoat::set_HSL( const double H, const double S, const double L )
{
    if ( H < 0.01 )
        throw std::runtime_error( "FingerCoxJephcoat::set_HSL(): H < 0.01." );
    if ( S < 0.01 )
        throw std::runtime_error( "FingerCoxJephcoat::set_HSL(): S < 0.01." );
    if ( L < 0.01 )
        throw std::runtime_error( "FingerCoxJephcoat::set_HSL(): L < 0.01." );
    A_ = H/L;
    B_ = S/L;
}

// ********************************************************************************

void FingerCoxJephcoat::set_A( const double A )
{
    if ( A < 0.0001 )
        throw std::runtime_error( "FingerCoxJephcoat::set_A(): A < 0.0001." );
    A_ = A;
}

// ********************************************************************************

void FingerCoxJephcoat::set_B( const double B )
{
    if ( B < 0.0001 )
        throw std::runtime_error( "FingerCoxJephcoat::set_B(): B < 0.0001." );
    B_ = B;
}

// ********************************************************************************

std::vector< double > FingerCoxJephcoat::asymmetric_peak( const Angle two_theta, const std::vector< Angle > & two_phi_values, const double FWHM ) const
{
    if ( two_theta > Angle::angle_90_degrees() )
        throw std::runtime_error( "FingerCoxJephcoat::asymmetric_peak(): two_theta > 90.0." );
    std::vector< double > result( two_phi_values.size(), 0.0 );
    Angle two_phi_min;
    Angle two_phi_infl;
    double cosine_argument = two_theta.cosine() * sqrt( square( A_ + B_ ) + 1.0 );
    if ( cosine_argument < 1.0 )
        two_phi_min = arccosine( cosine_argument );
    cosine_argument = two_theta.cosine() * sqrt( square( A_ - B_ ) + 1.0 );
    if ( cosine_argument < 1.0 )
        two_phi_infl = arccosine( cosine_argument );
    for ( size_t i( 0 ); i != two_phi_values.size(); ++i )
    {
        if ( two_phi_values[i] > Angle::angle_90_degrees() )
            throw std::runtime_error( "FingerCoxJephcoat::asymmetric_peak_H_is_S(): two_phi_values[i] > 90.0." );
        double numerator = 0.0;
        double denominator = 0.0;
        for ( size_t j( 0 ); j != N_; ++j )
        {
            Angle delta_n = ( two_theta + two_phi_min ) / 2.0 + ( two_theta - two_phi_min ) * x_i_[j] / 2.0;
            double C = sqrt( ( square( delta_n.cosine() ) / square( two_theta.cosine() ) ) - 1.0 );
            double term;
            if ( delta_n < two_phi_infl )
                term = ( ( A_ + B_ ) / C ) - 1.0;
            else
            {
                if ( A_ < B_ )
                    term = 2.0 * A_ / C;
                else
                    term = 2.0 * B_ / C;
            }
            term *= w_i_[j];
            term /= delta_n.cosine();
            numerator   += term * pseudo_Voigt( ( two_phi_values[i] - delta_n ).value_in_degrees(), FWHM, eta_ );
            denominator += term;
        }
        result[i] = numerator / denominator;
    }
    return result;
}

// ********************************************************************************

std::vector< double > FingerCoxJephcoat::asymmetric_peak_H_is_S( const Angle two_theta, const std::vector< Angle > & two_phi_values, const double FWHM ) const
{
    if ( two_theta > Angle::angle_90_degrees() )
        throw std::runtime_error( "FingerCoxJephcoat::asymmetric_peak(): two_theta > 90.0." );
    std::vector< double > result( two_phi_values.size(), 0.0 );
    Angle two_phi_min;
    double cosine_argument = two_theta.cosine() * sqrt( 4.0*square( A_ ) + 1.0 );
    if ( cosine_argument < 1.0 )
        two_phi_min = arccosine( cosine_argument );
    for ( size_t i( 0 ); i != two_phi_values.size(); ++i )
    {
        if ( two_phi_values[i] > Angle::angle_90_degrees() )
            throw std::runtime_error( "FingerCoxJephcoat::asymmetric_peak_H_is_S(): two_phi_values[i] > 90.0." );
        double numerator = 0.0;
        double denominator = 0.0;
        for ( size_t j( 0 ); j != N_; ++j )
        {
            Angle delta_n = ( two_theta + two_phi_min ) / 2.0 + ( two_theta - two_phi_min ) * x_i_[j] / 2.0;
            double C = sqrt( ( square( delta_n.cosine() ) / square( two_theta.cosine() ) ) - 1.0 );
            double term;
            if ( delta_n < two_theta )
                term = ( 2.0 * A_ / C ) - 1.0;
            else
                term = 2.0 * A_ / C;
            term *= w_i_[j] / delta_n.cosine();
            numerator   += term * pseudo_Voigt( ( two_phi_values[i] - delta_n ).value_in_degrees(), FWHM, eta_ );
            denominator += term;
        }
        result[i] = numerator / denominator;
    }
    return result;
}

// ********************************************************************************

