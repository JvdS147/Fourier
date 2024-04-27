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

#include "Mapping.h"
#include "Utilities.h"

#include <stdexcept>

// ********************************************************************************

Mapping::Mapping( const size_t nvalues )
{
    mapping_ = initialise_with_sequential_values( nvalues );
}

// ********************************************************************************

Mapping::Mapping( const std::vector< size_t > & values ):
mapping_(values)
{
    check_consistency();
}

// ********************************************************************************

void Mapping::swap( const size_t i, const size_t j )
{
    if ( ( i < mapping_.size() ) && ( j < mapping_.size() ) )
        std::swap( mapping_[i], mapping_[j] );
    else
        throw std::runtime_error( "Mapping::swap(): error: index out of range." );
}

// ********************************************************************************

void Mapping::push_back()
{
    mapping_.push_back( size() );
}

// ********************************************************************************

size_t Mapping::operator[]( const size_t i ) const
{
    if ( i < mapping_.size() )
        return mapping_[i];
    throw std::runtime_error( "Mapping::operator[]: error: index out of range." );
}

// ********************************************************************************

void Mapping::invert()
{
    std::vector< size_t > new_mapping( mapping_.size() );
    for ( size_t i( 0 ); i != mapping_.size(); ++i )
        new_mapping[ mapping_[i] ] = i;
    mapping_ = new_mapping;
}

// ********************************************************************************

void Mapping::check_consistency() const
{
    std::vector< bool > done( mapping_.size(), false );
    for ( size_t i( 0 ); i != mapping_.size(); ++i )
    {
        if ( ! ( mapping_[i] < mapping_.size() ) )
            throw std::runtime_error( "Mapping::check_consistency(): error: value out of range." );
        if ( done[ mapping_[i] ] )
            throw std::runtime_error( "Mapping::check_consistency(): error: value found more than once." );
        done[ mapping_[i] ] = true;
    }
}

// ********************************************************************************

