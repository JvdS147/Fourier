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

#include "MillerIndices.h"
#include "BasicMathsFunctions.h"
#include "Utilities.h"

#include <iostream>

// ********************************************************************************

MillerIndices::MillerIndices( const int h, const int k, const int l ): h_(h), k_(k), l_(l)
{
}

// ********************************************************************************

void MillerIndices::make_relative_prime()
{
    int gcd = greatest_common_divisor( greatest_common_divisor( h_, k_), l_ );
    if ( gcd == 0 )
        return;
    h_ /= gcd;
    k_ /= gcd;
    l_ /= gcd;
}

// ********************************************************************************

bool MillerIndices::is_000() const
{
    return ( ( h_ == 0 ) && ( k_ == 0 ) && ( l_ == 0 ) );
}

// ********************************************************************************

std::string MillerIndices::to_string() const
{
    return int2string( h() ) + " " + int2string( k() ) + " " + int2string( l() );
}

// ********************************************************************************

std::ostream & operator<<( std::ostream & os, const MillerIndices & miller_indices )
{
    os << miller_indices.h() << " " << miller_indices.k() << " " << miller_indices.l();
    return os;
}

// ********************************************************************************

bool operator<( const MillerIndices & lhs, const MillerIndices & rhs )
{
    if ( rhs.h() < lhs.h() )
        return true;
    if ( lhs.h() < rhs.h() )
        return false;
    if ( rhs.k() < lhs.k() )
        return true;
    if ( lhs.k() < rhs.k() )
        return false;
    return ( rhs.l() < lhs.l() );
}

// ********************************************************************************

bool operator==( const MillerIndices & lhs, const MillerIndices & rhs )
{
    return ( ( lhs.h() == rhs.h() ) &&
             ( lhs.k() == rhs.k() ) &&
             ( lhs.l() == rhs.l() ));
}

// ********************************************************************************

int operator*( const MillerIndices & lhs, const MillerIndices & rhs )
{
    return ( lhs.h() * rhs.h() + lhs.k() * rhs.k() + lhs.l() * rhs.l() );
}

// ********************************************************************************
