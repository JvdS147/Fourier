#ifndef FINGERCOXJEPHCOAT_H
#define FINGERCOXJEPHCOAT_H

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

class Angle;

#include <cstddef> // For definition of size_t
#include <vector>

/*
  An asymmetric peak.
  Currently only correctly implemented for two_theta and two_phi_values lower than 45 degrees.
  Generally speaking peak asymmetry is only visible up to about 30 degrees 2theta and beyond 150 degrees 2theta.
  we assume that we have not measured beyond 90 degrees 2theta.
  We do things properly here: express everything in terms of A = H/L and B = S/L and phi and theta are Angle objects.
  The peak shape is fixed to be pseudo-Voigt, but this is the case in all implemenetations I have seen and it is always mentioned how good the results are.
*/
class FingerCoxJephcoat
{
public:

    // Note that internally, this class only works with H/L and S/L.
    FingerCoxJephcoat( const double H, const double S, const double L );

    FingerCoxJephcoat( const double A, const double B );

    void set_HSL( const double H, const double S, const double L );

    // A = H/L
    double A() const { return A_; }
    void set_A( const double A );

    // B = S/L
    double B() const { return B_; }
    void set_B( const double B );

    // Returns the peak shape, normalised to an area of 1.0.
    // two_theta is the peak position.
    std::vector< double > asymmetric_peak( const Angle two_theta, const std::vector< Angle > & two_phi_values, const double FWHM ) const;

    // Sets B = A.
    std::vector< double > asymmetric_peak_H_is_S( const Angle two_theta, const std::vector< Angle > & two_phi_values, const double FWHM ) const;

private:
    double A_;
    double B_;
    double eta_;
    size_t N_;
    std::vector< double > x_i_;
    std::vector< double > w_i_;

};

#endif // FINGERCOXJEPHCOAT_H

