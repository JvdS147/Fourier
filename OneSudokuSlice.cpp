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

#include "OneSudokuSlice.h"

#include <stdexcept>
#include <iostream>

// ********************************************************************************

//OneSudokuSlice::OneSudokuSlice():
//id_(0)
//{
//    for ( size_t i( 0 ); i != 9; ++i )
//        values_.push_back( OneSudokuSquare() );
//}

// ********************************************************************************

//OneSudokuSlice::OneSudokuSlice( const size_t id, const SudokuSliceType type, const std::string & one_row ):
//id_(id),
//type_(type)
//{
//    for ( size_t i( 0 ); i != 9; ++i )
//        values_.push_back( OneSudokuSquare( string2integer( one_row.substr( i , 1 ) ) ) );
//}

// ********************************************************************************

OneSudokuSlice::OneSudokuSlice( const size_t id, const SudokuSliceType type, const std::vector< OneSudokuSquare > & values ):
id_(id),
type_(type)
{
    if ( ( values.size() != 6 ) && ( values.size() != 9 ) )
        throw std::runtime_error( "Sudoku::Sudoku( std::vector< std::string > ): incorrect number of squares." );
    values_ = values;
}

// ********************************************************************************

OneSudokuSquare OneSudokuSlice::square( const size_t i ) const
{
    return values_[i];
}

// ********************************************************************************

void OneSudokuSlice::set_square( const size_t i, const OneSudokuSquare & value )
{
    values_[i] = value;
}

// ********************************************************************************

bool OneSudokuSlice::solved() const
{
    if ( size() == 6 )
        std::cout << "OneSudokuSlice::solved(): we should never be here for size() == 6." << std::endl;
     for ( size_t i( 0 ); i != size(); ++i )
     {
         if ( ! values_[i].solved() )
            return false;
     }
     return true;
}

// ********************************************************************************

SetOfNumbers OneSudokuSlice::collect_solved_numbers() const
{
    SetOfNumbers result;
    for ( size_t i( 0 ); i != size(); ++i )
    {
        if ( values_[i].solved() )
            result.add( values_[i].value() );
    }
    return result;
}

// ********************************************************************************

SetOfNumbers OneSudokuSlice::collect_all_possible_numbers() const
{
    OneSudokuSquare merged( square( 0 ) );
    for ( size_t i( 1 ); i != size(); ++i )
        merged = merge( merged, square( i ) );
    return merged.values();
}

// ********************************************************************************

std::vector< size_t > OneSudokuSlice::indices_of_unsolved_squares() const
{
    std::vector< size_t > result;
    for ( size_t i( 0 ); i != size(); ++i )
    {
        if ( ! values_[i].solved() )
            result.push_back( i );
    }
    return result;
}

// ********************************************************************************

std::vector< size_t > OneSudokuSlice::indices_of_solved_squares() const
{
    std::vector< size_t > result;
    for ( size_t i( 0 ); i != size(); ++i )
    {
        if ( values_[i].solved() )
            result.push_back( i );
    }
    return result;
}

// ********************************************************************************

void OneSudokuSlice::show() const
{
    for ( size_t i( 0 ); i != size(); ++i )
    {
        values_[i].show();
        std::cout << " | ";
    }
    std::cout << std::endl;
}

// ********************************************************************************

