/* *********************************************
Copyright (c) 2013-2020, Cornelis Jan (Jacco) van de Streek
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

#include "RandomNumberGenerator.h"

#include <cstddef> // For definition of size_t

// ********************************************************************************

RandomNumberGenerator_integer::RandomNumberGenerator_integer( const int idum )
{
    IM1 = 2147483563; // 2,147,483,563
    IM2 = 2147483399;
    IA1 = 40014;
    IA2 = 40692;
    IQ1 = 53668;
    IQ2 = 52774;
    IR1 = 12211;
    IR2 = 3791;
    NTAB = 32;
    IMM1 = IM1 - 1;
    NDIV = 1 + IMM1/NTAB;
    iv_ = std::vector< int >(32,0);
    
    idum_ = idum;
    if ( idum_ < 0 )
        idum_ = - idum_;
    if ( idum_ == 0 )
        idum_ = 1539;
    idum2_ = idum_;
    for ( size_t j( 0 ); j != 8; ++j )
    {
        int k = idum_/IQ1;
        idum_ = IA1*(idum_-k*IQ1)-k*IR1;
        if ( idum_ < 0 )
            idum_ += IM1;
    }
    for ( size_t j( 0 ); j != 32; ++j )
    {
        int k = idum_/IQ1;
        idum_ = IA1*(idum_-k*IQ1)-k*IR1;
        if ( idum_ < 0 )
            idum_ += IM1;
        iv_[j] = idum_;
    }
    iy_ = iv_[0];
}

// ********************************************************************************

int RandomNumberGenerator_integer::next_number( const int start, const int end )
{
    int k = idum_ / IQ1;
    idum_ = IA1*(idum_-k*IQ1)-k*IR1;
    if ( idum_ < 0 )
        idum_ += IM1;
    k = idum2_ / IQ2;
    idum2_ = IA2*(idum2_-k*IQ2)-k*IR2;
    if ( idum2_ < 0 )
        idum2_ += IM2;
    int j = iy_ / NDIV;
    iy_ = iv_[j] - idum2_;
    iv_[j] = idum_;
    if ( iy_ < 0 )
        iy_ += IM1;
    return start + ( iy_ % ( 1 + end - start ) );
}

// ********************************************************************************

RandomNumberGenerator_double::RandomNumberGenerator_double( const int idum )
{
    IM1 = 2147483563;
    IM2 = 2147483399;
    IA1 = 40014;
    IA2 = 40692;
    IQ1 = 53668;
    IQ2 = 52774;
    IR1 = 12211;
    IR2 = 3791;
    NTAB = 32;
    IMM1 = IM1 - 1;
    NDIV = 1 + IMM1/NTAB;
    iv_ = std::vector< int >(32,0);
    AM = 1.0 / static_cast<double>( IM1 );
    
    idum_ = idum;
    if ( idum_ < 0 )
        idum_ = - idum_;
    if ( idum_ == 0 )
        idum_ = 1539;
    idum2_ = idum_;
    for ( size_t j( 0 ); j != 8; ++j )
    {
        int k = idum_/IQ1;
        idum_ = IA1*(idum_-k*IQ1)-k*IR1;
        if ( idum_ < 0 )
            idum_ += IM1;
    }
    for ( size_t j( 0 ); j != 32; ++j )
    {
        int k = idum_/IQ1;
        idum_ = IA1*(idum_-k*IQ1)-k*IR1;
        if ( idum_ < 0 )
            idum_ += IM1;
        iv_[j] = idum_;
    }
    iy_ = iv_[0];
}

// ********************************************************************************

double RandomNumberGenerator_double::next_number()
{
    int k = idum_ / IQ1;
    idum_ = IA1*(idum_-k*IQ1)-k*IR1;
    if ( idum_ < 0 )
        idum_ += IM1;
    k = idum2_ / IQ2;
    idum2_ = IA2*(idum2_-k*IQ2)-k*IR2;
    if ( idum2_ < 0 )
        idum2_ += IM2;
    int j = iy_ / NDIV;
    iy_ = iv_[j] - idum2_;
    iv_[j] = idum_;
    if ( iy_ < 0 )
        iy_ += IM1;
    return AM*iy_;
}

// ********************************************************************************

