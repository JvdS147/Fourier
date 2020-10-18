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

BagOfNumbers::BagOfNumbers( const size_t nvalues, const int idum ):
set_of_numbers_(nvalues),
RNG_int_(idum)
{
}

// ********************************************************************************

void BagOfNumbers::set_duplicates_policy( const SetOfNumbers::DuplicatesPolicy duplicates_policy )
{
    set_of_numbers_.set_duplicates_policy( duplicates_policy );
}

// ********************************************************************************

void BagOfNumbers::remove( const size_t value )
{
    set_of_numbers_.remove( value );
}

// ********************************************************************************

void BagOfNumbers::add( const size_t value )
{
    set_of_numbers_.add( value );
}

// ********************************************************************************

size_t BagOfNumbers::draw()
{
    size_t result = draw_with_replace();
    set_of_numbers_.remove( result );
    return result;
}

// ********************************************************************************

size_t BagOfNumbers::draw_with_replace() const
{
    if ( set_of_numbers_.empty() )
        throw std::runtime_error( "BagOfNumbers::draw_with_replace(): bag is empty." );
    size_t iPos = RNG_int_.next_number( 0, set_of_numbers_.size() - 1 );
    return set_of_numbers_.value( iPos );
}

// ********************************************************************************

