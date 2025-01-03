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

#include "SphericalHarmonics.h"
#include "MathsFunctions.h"

#include <stdexcept>

// ********************************************************************************

Complex spherical_harmonics( const size_t l, const int m, const Angle alpha, const Angle beta )
{
    return sqrt( ( (2.0*l+1.0) / (4.0*CONSTANT_PI) ) ) * Racah_spherical_harmonics( l, m, alpha, beta );
}

// ********************************************************************************

Complex Racah_spherical_harmonics( const size_t l, const int m, const Angle alpha, const Angle beta )
{
    if ( absolute( m ) > l )
        throw std::runtime_error( "Racah_spherical_harmonics(): |m| > l." );
    if ( m < 0 )
    {
        Complex result = Racah_spherical_harmonics( l, -m, alpha, beta );
        result.conjugate();
        return is_even( -m ) ? result : -result;
    }
    double factorial_lmm = factorial(l-m);
    double factorial_lpm = factorial(l+m);
    double result = sqrt( factorial_lmm / factorial_lpm );
    result *= associated_Legendre_polynomial( l, m, beta.cosine() );
    Complex z( 0.0, m * alpha.value_in_radians() );
    return result * exponential( z );
}

// ********************************************************************************

