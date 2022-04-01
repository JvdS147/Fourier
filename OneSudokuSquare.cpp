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

#include "OneSudokuSquare.h"
#include "Utilities.h"

#include <stdexcept>
#include <iostream>

// ********************************************************************************

OneSudokuSquare::OneSudokuSquare():
values_( 1, 9, SetOfNumbers::THROW )
{
    values_.set_empty_is_allowed( false );
}

// ********************************************************************************

OneSudokuSquare::OneSudokuSquare( const size_t value )
{
    if ( value == 0 )
        *this = OneSudokuSquare();
    else
    {
        if ( ( value < 1 ) || ( value > 9 ) )
            throw std::runtime_error( "OneSudokuSquare::OneSudokuSquare( size_t ): incorrect value: " + int2string( value ) );
        values_.add( value );
        values_.set_duplicates_policy( SetOfNumbers::THROW );
        values_.set_empty_is_allowed( false );
    }
}

// ********************************************************************************

OneSudokuSquare::OneSudokuSquare( const SetOfNumbers & values ):
values_(values)
{
    values_.set_duplicates_policy( SetOfNumbers::THROW );
    values_.set_empty_is_allowed( false );
    if ( values_.value( 0 ) < 1 )
        throw std::runtime_error( "OneSudokuSquare::OneSudokuSquare( std::vector< size_t > ): incorrect value: " + int2string( values_.value( 0 ) ) );
    if ( values_.value( values_.size() - 1 ) > 9 )
        throw std::runtime_error( "OneSudokuSquare::OneSudokuSquare( std::vector< size_t > ): incorrect value: " + int2string( values_.value( values_.size() - 1 ) ) );    
}

// ********************************************************************************

size_t OneSudokuSquare::value( const size_t i ) const
{
    return values_.value( i );
}

// ********************************************************************************

size_t OneSudokuSquare::value() const
{
    if ( ! this->solved() )
        throw std::runtime_error( "OneSudokuSquare::value(): square is not solved." );
    return values_.value( 0 );
}

// ********************************************************************************

bool OneSudokuSquare::set( const size_t value )
{
    if ( this->contains( value ) )
        return false;
    values_.add( value );
    return true;
}

// ********************************************************************************

bool OneSudokuSquare::unset( const size_t value )
{
    if ( ! this->contains( value ) )
        return false;
    values_.remove( value );
    return true;
}

// ********************************************************************************

bool OneSudokuSquare::unset( const OneSudokuSquare & sudoku_square )
{
    SetOfNumbers new_values( SetOfNumbers::THROW );
    for ( size_t i( 0 ); i != values_.size(); ++i )
    {
        if ( ! sudoku_square.contains( values_.value( i ) ) )
            new_values.add( values_.value( i ) );
    }
    new_values.set_empty_is_allowed( false );
    bool result = ( values_ != new_values );
    values_ = new_values;
    return result;
}

// ********************************************************************************

bool OneSudokuSquare::unset( const SetOfNumbers & values )
{
    if ( values.empty() )
        return false;
    return this->unset( OneSudokuSquare( values ) );
}

// ********************************************************************************

bool OneSudokuSquare::contains( const size_t value ) const
{
    return values_.contains( value );
}

// ********************************************************************************

size_t OneSudokuSquare::size() const
{
    return values_.size();
}

// ********************************************************************************

bool OneSudokuSquare::solved() const
{
    return ( values_.size() == 1 );
}

// ********************************************************************************

void OneSudokuSquare::show() const
{
    values_.show();
}

// ********************************************************************************

OneSudokuSquare merge( const OneSudokuSquare & lhs, const OneSudokuSquare & rhs )
{
    SetOfNumbers set_of_numbers( lhs.values() );
    set_of_numbers.set_duplicates_policy( SetOfNumbers::AUTO_REMOVE );
    set_of_numbers.add( rhs.values() );
    set_of_numbers.set_duplicates_policy( SetOfNumbers::THROW );
    return OneSudokuSquare( set_of_numbers );
}

// ********************************************************************************

