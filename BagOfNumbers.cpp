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

#include "BagOfNumbers.h"

#include <cstdlib>
#include <stdexcept>

// ********************************************************************************

BagOfNumbers::BagOfNumbers( const int idum ):RNG_int_(idum)
{
}
// ********************************************************************************

BagOfNumbers::BagOfNumbers( const size_t nnumbers, const int idum ):RNG_int_(idum)
{
    values_.reserve( nnumbers );
    for ( size_t i( 0 ); i != nnumbers; ++i )
        values_.push_back( i );
}

// ********************************************************************************

void BagOfNumbers::remove( const int value )
{
    for ( size_t i(0); i != size(); ++i )
    {
        if ( values_[i] == value )
        {
            remove_position( i );
            return;
        }
    }
}

// ********************************************************************************

void BagOfNumbers::add( const int value )
{
    values_.push_back( value );
}

// ********************************************************************************

size_t BagOfNumbers::draw()
{
    if ( values_.empty() )
        throw std::runtime_error( "BagOfNumbers::draw(): bag is empty." );
    size_t iPos = RNG_int_.next_number( 0, values_.size() - 1 );
    size_t result = values_[iPos];
    remove_position( iPos );
    return result;
}

// ********************************************************************************

size_t BagOfNumbers::draw_with_replace()
{
    if ( values_.empty() )
        throw std::runtime_error( "BagOfNumbers::draw_with_replace(): bag is empty." );
    size_t iPos = RNG_int_.next_number( 0, values_.size() - 1 );
    return values_[iPos];
}

// ********************************************************************************

void BagOfNumbers::remove_position( const size_t iPos )
{
    if ( values_.empty() )
        return;
    values_[iPos] = values_[size()-1];
    values_.pop_back();
}

// ********************************************************************************

