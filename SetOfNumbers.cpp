/* *********************************************
Copyright (c) 2013-2023, Cornelis Jan (Jacco) van de Streek
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

#include "SetOfNumbers.h"
#include "Sort.h"

#include <iostream>
#include <stdexcept>

// ********************************************************************************

SetOfNumbers::SetOfNumbers( const DuplicatesPolicy duplicates_policy ):
duplicates_policy_(duplicates_policy),
empty_is_allowed_(true)
{
}

// ********************************************************************************

SetOfNumbers::SetOfNumbers( const size_t nvalues, const DuplicatesPolicy duplicates_policy ):
duplicates_policy_(duplicates_policy),
empty_is_allowed_(true)
{
    values_.reserve( nvalues );
    for ( size_t i( 0 ); i != nvalues; ++i )
        values_.push_back( i );
    sorted_map_ = Mapping( nvalues );
    is_sorted_ = true;
}

// ********************************************************************************

// Fills the set with the numbers [begin,end]
SetOfNumbers::SetOfNumbers( const size_t begin, const size_t end, const DuplicatesPolicy duplicates_policy ):
duplicates_policy_(duplicates_policy),
empty_is_allowed_(true)
{
    if ( end < begin )
        throw std::runtime_error( "SetOfNumbers::SetOfNumbers(): end < begin." );
    size_t nvalues( end-begin+1 );
    values_.reserve( nvalues );
    for ( size_t i( 0 ); i != nvalues; ++i )
        values_.push_back( begin+i );
    sorted_map_ = Mapping( nvalues );
    is_sorted_ = true;
}

// ********************************************************************************

SetOfNumbers::SetOfNumbers( const std::vector< size_t > & values, const DuplicatesPolicy duplicates_policy ):
duplicates_policy_(duplicates_policy),
empty_is_allowed_(true)
{
    values_ = values;
    is_sorted_ = false;
    check_for_duplicates();
}

// ********************************************************************************

void SetOfNumbers::set_duplicates_policy( const DuplicatesPolicy duplicates_policy )
{
    duplicates_policy_ = duplicates_policy;
    check_for_duplicates();
}

// ********************************************************************************

void SetOfNumbers::set_empty_is_allowed( const bool empty_is_allowed )
{
    empty_is_allowed_ = empty_is_allowed;
    check_if_empty();
}

// ********************************************************************************

size_t SetOfNumbers::value( const size_t index ) const
{
    if ( ! ( index < values_.size() ) )
        throw std::runtime_error( "SetOfNumbers::value(): index out of range." );
    sort();
    return values_[ sorted_map_[ index ] ];
}

// ********************************************************************************

std::vector< size_t > SetOfNumbers::values() const
{
    std::vector< size_t > result;
    result.reserve( size() );
    for ( size_t i( 0 ); i != size(); ++i )
        result.push_back( value( i ) );
    return result;
}

// ********************************************************************************

SetOfNumbers SetOfNumbers::unique_values() const
{
    return SetOfNumbers( values(), AUTO_REMOVE );
}

// ********************************************************************************

// Returns how often value occurs in the set.
size_t SetOfNumbers::frequency( const size_t value ) const
{
    size_t result( 0 );
    for ( size_t i( 0 ); i != size(); ++i )
    {
        if ( values_[i] == value )
            ++result;
    }
    return result;
}

// ********************************************************************************

void SetOfNumbers::add( const size_t value )
{
    values_.push_back( value );
    is_sorted_ = false;
    check_for_duplicates();
}

// ********************************************************************************

void SetOfNumbers::remove( const size_t value )
{
    for ( size_t i( 0 ); i != size(); ++i )
    {
        if ( values_[i] == value )
        {
            remove_position( i );
            check_if_empty();
            return;
        }
    }
}

// ********************************************************************************

bool SetOfNumbers::contains( const size_t value ) const
{
    for ( size_t i( 0 ); i != size(); ++i )
    {
        if ( values_[i] == value )
            return true;
    }
    return false;
}

// ********************************************************************************

bool SetOfNumbers::contains_duplicates() const
{
    if ( duplicates_policy_ != ALLOWED )
        return false;
    sort();
    for ( size_t i( 1 ); i < size(); ++i )
    {
        if ( values_[ sorted_map_[ i ] ] == values_[ sorted_map_[ i-1 ] ] )
            return true;
    }
    return false;
}

// ********************************************************************************

// The attributes duplicates_allowed_ and empty_is_allowed_ of the
// argument values are ignored, the values of *this are kept.
void SetOfNumbers::add( const SetOfNumbers & set_of_numbers )
{
    for ( size_t i( 0 ); i != set_of_numbers.size(); ++i )
        this->add( set_of_numbers.value( i ) );
}

// ********************************************************************************

// The attributes duplicates_allowed_ and empty_is_allowed_ of the
// argument values are ignored, the values of *this are kept.
void SetOfNumbers::remove( const SetOfNumbers & set_of_numbers )
{
    for ( size_t i( 0 ); i != set_of_numbers.size(); ++i )
        this->remove( set_of_numbers.value( i ) );
}

// ********************************************************************************

SetOfNumbers SetOfNumbers::in_common( const SetOfNumbers & set_of_numbers )
{
    SetOfNumbers result( duplicates_policy() );
    for ( size_t i( 0 ); i != size(); ++i )
    {
        if ( result.contains( value( i ) ) )
            continue;
        if ( ! set_of_numbers.contains( value( i ) ) )
            continue;
        size_t nvalues = std::min( frequency( value( i ) ), set_of_numbers.frequency( value( i ) ) );
        for ( size_t j( 0 ); j != nvalues; ++j )
            result.add( value( i ) );
    }
    return result;
}

// ********************************************************************************

bool SetOfNumbers::operator==( const SetOfNumbers & rhs ) const
{
    if ( this->size() != rhs.size() )
        return false;
    if ( this->duplicates_policy() != rhs.duplicates_policy() )
        return false;
    if ( this->empty_is_allowed() != rhs.empty_is_allowed() )
        return false;
    for ( size_t i( 0 ); i != this->size(); ++i )
    {
        if ( this->value( i ) != rhs.value( i ) )
            return false;
    }
    return true;
}

// ********************************************************************************

bool SetOfNumbers::operator!=( const SetOfNumbers & rhs ) const
{
    return ( ! ( this->values_ == rhs.values_ ) ); 
}
    
// ********************************************************************************

void SetOfNumbers::show() const
{
    for ( size_t i( 0 ); i != size(); ++i )
        std::cout << value( i ) << " ";
    std::cout << std::endl;
}

// ********************************************************************************

void SetOfNumbers::remove_position( const size_t index )
{
    if ( index < values_.size() )
    {
        is_sorted_ = false;
        values_[index] = values_[size()-1];
        values_.pop_back();
    }
    else
        throw std::runtime_error( "SetOfNumbers::remove_position(): index out of range." );
    
}

// ********************************************************************************

void SetOfNumbers::sort() const
{
    if ( is_sorted_ )
        return;
    sorted_map_ = ::sort( values_ );
    is_sorted_ = true;
}

// ********************************************************************************

void SetOfNumbers::check_for_duplicates()
{
    if ( duplicates_policy_ == ALLOWED )
        return;
    sort();
    for ( size_t i( 1 ); i < size(); ++i )
    {
        if ( values_[ sorted_map_[ i ] ] == values_[ sorted_map_[ i-1 ] ] )
        {
            if ( duplicates_policy_ == THROW )
                throw std::runtime_error( "SetOfNumbers: duplicates found." );
            else
            {
                remove_position( sorted_map_[ i ] );
                check_for_duplicates();
                return;
            }
        }
    }
}

// ********************************************************************************

void SetOfNumbers::check_if_empty() const
{
    if ( empty_is_allowed_ )
        return;
    if ( empty() )
        throw std::runtime_error( "SetOfNumbers: set is empty." );
}

// ********************************************************************************

// The attributes duplicates_allowed_ and empty_is_allowed_ of the
// arguments are ignored
// duplicates_allowed_ is set to true, empty_is_allowed_ is set to true
SetOfNumbers merge( const SetOfNumbers & lhs, const SetOfNumbers & rhs )
{
    SetOfNumbers result;
    result.add( lhs );
    result.add( rhs );
    return result;
}

// ********************************************************************************

