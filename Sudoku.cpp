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

#include "Sudoku.h"
#include "OneSudokuSlice.h"
#include "SetOfNumbers.h"
#include "Sort.h"
#include "TransSquareDependency.h"
#include "Utilities.h"

#include <stdexcept>
#include <iostream>

static const size_t sudoku_slice_mapping[27][9] =
                   {
                       {  0,  1,  2,  3,  4,  5,  6,  7,  8 },
                       {  9, 10, 11, 12, 13, 14, 15, 16, 17 },
                       { 18, 19, 20, 21, 22, 23, 24, 25, 26 },
                       { 27, 28, 29, 30, 31, 32, 33, 34, 35 },
                       { 36, 37, 38, 39, 40, 41, 42, 43, 44 },
                       { 45, 46, 47, 48, 49, 50, 51, 52, 53 },
                       { 54, 55, 56, 57, 58, 59, 60, 61, 62 },
                       { 63, 64, 65, 66, 67, 68, 69, 70, 71 },
                       { 72, 73, 74, 75, 76, 77, 78, 79, 80 },
                       
                       {  0,  9, 18, 27, 36, 45, 54, 63, 72 },
                       {  1, 10, 19, 28, 37, 46, 55, 64, 73 },
                       {  2, 11, 20, 29, 38, 47, 56, 65, 74 },
                       {  3, 12, 21, 30, 39, 48, 57, 66, 75 },
                       {  4, 13, 22, 31, 40, 49, 58, 67, 76 },
                       {  5, 14, 23, 32, 41, 50, 59, 68, 77 },
                       {  6, 15, 24, 33, 42, 51, 60, 69, 78 },
                       {  7, 16, 25, 34, 43, 52, 61, 70, 79 },
                       {  8, 17, 26, 35, 44, 53, 62, 71, 80 },
                       
                       {  0,  1,  2,  9, 10, 11, 18, 19, 20 },
                       {  3,  4,  5, 12, 13, 14, 21, 22, 23 },
                       {  6,  7,  8, 15, 16, 17, 24, 25, 26 },
                       { 27, 28, 29, 36, 37, 38, 45, 46, 47 },
                       { 30, 31, 32, 39, 40, 41, 48, 49, 50 },
                       { 33, 34, 35, 42, 43, 44, 51, 52, 53 },
                       { 54, 55, 56, 63, 64, 65, 72, 73, 74 },
                       { 57, 58, 59, 66, 67, 68, 75, 76, 77 },
                       { 60, 61, 62, 69, 70, 71, 78, 79, 80 }
                   };

static size_t sudoku_determining_values_mapping[54][6] =
                   {
                       {   9, 10, 11, 18, 19, 20 },
                       {  12, 13, 14, 21, 22, 23 },
                       {  15, 16, 17, 24, 25, 26 },
                       {   0,  1,  2, 18, 19, 20 },
                       {   3,  4,  5, 21, 22, 23 },
                       {   6,  7,  8, 24, 25, 26 },
                       {   0,  1,  2,  9, 10, 11 },
                       {   3,  4,  5, 12, 13, 14 },
                       {   6,  7,  8, 15, 16, 17 },
                       {  36, 37, 38, 45, 46, 47 },
                       {  39, 40, 41, 48, 49, 50 },
                       {  42, 43, 44, 51, 52, 53 },
                       {  27, 28, 29, 45, 46, 47 },
                       {  30, 31, 32, 48, 49, 50 },
                       {  33, 34, 35, 51, 52, 53 },
                       {  27, 28, 29, 36, 37, 38 },
                       {  30, 31, 32, 39, 40, 41 },
                       {  33, 34, 35, 42, 43, 44 },
                       {  63, 64, 65, 72, 73, 74 },
                       {  66, 67, 68, 75, 76, 77 },
                       {  69, 70, 71, 78, 79, 80 },
                       {  54, 55, 56, 72, 73, 74 },
                       {  57, 58, 59, 75, 76, 77 },
                       {  60, 61, 62, 78, 79, 80 },
                       {  54, 55, 56, 63, 64, 65 },
                       {  57, 58, 59, 66, 67, 68 },
                       {  60, 61, 62, 69, 70, 71 },
                       {   1,  2, 10, 11, 19, 20 },
                       {   0,  2,  9, 11, 18, 20 },
                       {   0,  1,  9, 10, 18, 19 },
                       {   4,  5, 13, 14, 22, 23 },
                       {   3,  5, 12, 14, 21, 23 },
                       {   3,  4, 12, 13, 21, 22 },
                       {   7,  8, 16, 17, 25, 26 },
                       {   6,  8, 15, 17, 24, 26 },
                       {   6,  7, 15, 16, 24, 25 },
                       {  28, 29, 37, 38, 46, 47 },
                       {  27, 29, 36, 38, 45, 47 },
                       {  27, 28, 36, 37, 45, 46 },
                       {  31, 32, 40, 41, 49, 50 },
                       {  30, 32, 39, 41, 48, 50 },
                       {  30, 31, 39, 40, 48, 49 },
                       {  34, 35, 43, 44, 52, 53 },
                       {  33, 35, 42, 44, 51, 53 },
                       {  33, 34, 42, 43, 51, 52 },
                       {  55, 56, 64, 65, 73, 74 },
                       {  54, 56, 63, 65, 72, 74 },
                       {  54, 55, 63, 64, 72, 73 },
                       {  58, 59, 67, 68, 76, 77 },
                       {  57, 59, 66, 68, 75, 77 },
                       {  57, 58, 66, 67, 75, 76 },
                       {  61, 62, 70, 71, 79, 80 },
                       {  60, 62, 69, 71, 78, 80 },
                       {  60, 61, 69, 70, 78, 79 }
                   };

static size_t sudoku_values_to_be_changed_mapping[54][6] =
                   {
                       {   3,  4,  5,  6,  7,  8 },
                       {   0,  1,  2,  6,  7,  8 },
                       {   0,  1,  2,  3,  4,  5 },
                       {  12, 13, 14, 15, 16, 17 },
                       {   9, 10, 11, 15, 16, 17 },
                       {   9, 10, 11, 12, 13, 14 },
                       {  21, 22, 23, 24, 25, 26 },
                       {  18, 19, 20, 24, 25, 26 },
                       {  18, 19, 20, 21, 22, 23 },
                       {  30, 31, 32, 33, 34, 35 },
                       {  27, 28, 29, 33, 34, 35 },
                       {  27, 28, 29, 30, 31, 32 },
                       {  39, 40, 41, 42, 43, 44 },
                       {  36, 37, 38, 42, 43, 44 },
                       {  36, 37, 38, 39, 40, 41 },
                       {  48, 49, 50, 51, 52, 53 },
                       {  45, 46, 47, 51, 52, 53 },
                       {  45, 46, 47, 48, 49, 50 },
                       {  57, 58, 59, 60, 61, 62 },
                       {  54, 55, 56, 60, 61, 62 },
                       {  54, 55, 56, 57, 58, 59 },
                       {  66, 67, 68, 69, 70, 71 },
                       {  63, 64, 65, 69, 70, 71 },
                       {  63, 64, 65, 66, 67, 68 },
                       {  75, 76, 77, 78, 79, 80 },
                       {  72, 73, 74, 78, 79, 80 },
                       {  72, 73, 74, 75, 76, 77 },
                       {  27, 36, 45, 54, 63, 72 },
                       {  28, 37, 46, 55, 64, 73 },
                       {  29, 38, 47, 56, 65, 74 },
                       {  30, 39, 48, 57, 66, 75 },
                       {  31, 40, 49, 58, 67, 76 },
                       {  32, 41, 50, 59, 68, 77 },
                       {  33, 42, 51, 60, 69, 78 },
                       {  34, 43, 52, 61, 70, 79 },
                       {  35, 44, 53, 62, 71, 80 },
                       {   0,  9, 18, 54, 63, 72 },
                       {   1, 10, 19, 55, 64, 73 },
                       {   2, 11, 20, 56, 65, 74 },
                       {   3, 12, 21, 57, 66, 75 },
                       {   4, 13, 22, 58, 67, 76 },
                       {   5, 14, 23, 59, 68, 77 },
                       {   6, 15, 24, 60, 69, 78 },
                       {   7, 16, 25, 61, 70, 79 },
                       {   8, 17, 26, 62, 71, 80 },
                       {   0,  9, 18, 27, 36, 45 },
                       {   1, 10, 19, 28, 37, 46 },
                       {   2, 11, 20, 29, 38, 47 },
                       {   3, 12, 21, 30, 39, 48 },
                       {   4, 13, 22, 31, 40, 49 },
                       {   5, 14, 23, 32, 41, 50 },
                       {   6, 15, 24, 33, 42, 51 },
                       {   7, 16, 25, 34, 43, 52 },
                       {   8, 17, 26, 35, 44, 53 }
                   };

// ********************************************************************************

Sudoku::Sudoku()
{
    current_slice_ = CyclicInteger( 0, number_of_slices()-1, number_of_slices()-1 );
    current_trans_square_dependency_ = CyclicInteger( 0, 53, 53 );
    values_.reserve( 81 );
    for ( size_t i( 0 ); i != 81; ++i )
        values_.push_back( OneSudokuSquare() );
}

// ********************************************************************************

Sudoku::Sudoku( const std::vector< std::string > & rows )
{
    current_slice_ = CyclicInteger( 0, number_of_slices()-1, number_of_slices()-1 );
    current_trans_square_dependency_ = CyclicInteger( 0, 53, 53 );
    values_.reserve( 81 );
    for ( size_t i( 0 ); i != 81; ++i )
        values_.push_back( OneSudokuSquare() );
    if ( rows.size() != 9 )
        throw std::runtime_error( "Sudoku::Sudoku( std::vector< std::string > ): incorrect number of rows." );
    size_t index( 0 );
    for ( size_t i( 0 ); i != rows.size(); ++i )
    {
        if ( rows[i].size() != 9 )
            throw std::runtime_error( "Sudoku::Sudoku( std::vector< std::string > ): incorrect number of squares." );
        for ( size_t j( 0 ); j != rows.size(); ++j )
        {
            int value = string2integer( rows[i].substr( j, 1 ) );
            if ( value != 0 )
                this->update_square( index, OneSudokuSquare( value ) );
            ++index;
        }
    }
}

// ********************************************************************************

OneSudokuSlice Sudoku::slice( const size_t i ) const
{
    if ( ! ( i < number_of_slices() ) )
        throw std::runtime_error( "Sudoku::slice( size_t ): invalid slice number." );
    return OneSudokuSlice( i, slice_id2type( i ), get_slice_values( i ) );
}

// ********************************************************************************

OneSudokuSlice Sudoku::row( const size_t i ) const
{
    if ( ! ( i < 9 ) )
        throw std::runtime_error( "Sudoku::row( size_t ): invalid row number." );
    return OneSudokuSlice( i, slice_id2type( i ), get_slice_values( i ) );
}

// ********************************************************************************

OneSudokuSlice Sudoku::column( const size_t i ) const
{
    if ( ! ( i < 9 ) )
        throw std::runtime_error( "Sudoku::column( size_t ): invalid column number." );
    return OneSudokuSlice( i+9, slice_id2type( i+9 ), get_slice_values( i+9 ) );
}

// ********************************************************************************

OneSudokuSlice Sudoku::block( const size_t i ) const
{
    if ( ! ( i < 9 ) )
        throw std::runtime_error( "Sudoku::block( size_t ): invalid block number." );
    return OneSudokuSlice( i+18, slice_id2type( i+18 ), get_slice_values( i+18 ) );
}

// ********************************************************************************

OneSudokuSlice Sudoku::next_slice()
{
    ++current_slice_;
    OneSudokuSlice result = slice( current_slice_.current_value() );
    // This can go wrong, partially because sudoku.solved() returns false
    // when the sudoku has been solved but contains a contradiction,
    // but this is not true for a slice. So sudoku.solved() may return false when all slices return true
    // We may call this while the sudoku has been solved (due to the recursion in update_square()), in which case
    // we would end up with an infinite loop.
    size_t infinite_loop_guard( 0 );
    while ( result.solved() )
    {
        ++current_slice_;
        result = slice( current_slice_.current_value() );
        ++infinite_loop_guard;
        if ( infinite_loop_guard == number_of_slices() )
            break;
    }
    return result;
}

// ********************************************************************************

void Sudoku::update_slice( const OneSudokuSlice & slice )
{
    if ( slice.size() != 9 )
        throw std::runtime_error( "Sudoku::update_slice( OneSudokuSlice ): slice should have size() == 9." );    
    for ( size_t i( 0 ); i != slice.size(); ++i )
        update_square( sudoku_slice_mapping[slice.id()][i], slice.square( i ) );
}

// ********************************************************************************

bool Sudoku::update_square( const size_t i, OneSudokuSquare square )
{
    // Due to the recursion, it is possible that a square that we were in the middle of dealing with
    // now becomes active again, but has meanwhile been solved
    if ( values_[i].solved() )
        return false;
    if ( values_[i] == square )
        return false;
    SetOfNumbers new_values = square.values();
    // Check that we do not add values that before were impossible
    new_values = new_values.in_common( values_[i].values() );
    square = OneSudokuSquare( new_values );
    values_[i] = square;
    if ( square.solved() )
    {
        // Collect the indices of the unsolved squares affected by this square.
        std::vector< size_t > affected_squares_indices;
        std::vector< size_t > slice_indices = square2slices( i );
        for ( size_t iSlice( 0 ); iSlice != slice_indices.size(); ++iSlice )
        {
            for ( size_t j( 0 ); j != 9; ++j )
            {
                if ( ! values_[ sudoku_slice_mapping[slice_indices[iSlice]][j] ].solved() )
                    affected_squares_indices.push_back( sudoku_slice_mapping[slice_indices[iSlice]][j] );
            }
        }
        // Loop over all squares and if not solved, remove the square we just solved by updating the square
        // (so that if that quare is now solved, it triggers another update)
        for ( size_t j( 0 ); j != affected_squares_indices.size(); ++j )
        {
            OneSudokuSquare square_2 = values_[ affected_squares_indices[j] ];
            if ( square_2.solved() )
                continue;
            if ( ! square_2.contains( square.value() ) )
                continue;
            this->unset( affected_squares_indices[j], square.value() );
        }
    }
    return true;
}

// ********************************************************************************

bool Sudoku::unset( const size_t i, const size_t value )
{
    OneSudokuSquare square = values_[i];
    square.unset( value );
    return update_square( i, square );
}

// ********************************************************************************

TransSquareDependency Sudoku::next_trans_square_dependency()
{
    ++current_trans_square_dependency_;
    std::vector< OneSudokuSquare > determining_values;
    for ( size_t i( 0 ); i != 6; ++i )
        determining_values.push_back( values_[ sudoku_determining_values_mapping[current_trans_square_dependency_.current_value()][i] ] );
    std::vector< OneSudokuSquare > values_to_be_changed;
    for ( size_t i( 0 ); i != 6; ++i )
        values_to_be_changed.push_back( values_[ sudoku_values_to_be_changed_mapping[current_trans_square_dependency_.current_value()][i] ] );
    return TransSquareDependency( current_trans_square_dependency_.current_value(), determining_values, values_to_be_changed );
}

// ********************************************************************************

// To be called immediately after next_trans_square_dependency().
void Sudoku::update_trans_square_dependency( const TransSquareDependency & tsd )
{
    OneSudokuSlice values_to_be_changed = tsd.values_to_be_changed();
    for ( size_t i( 0 ); i != values_to_be_changed.size(); ++i )
        update_square( sudoku_values_to_be_changed_mapping[tsd.id()][i], values_to_be_changed.square( i ) );
}

// ********************************************************************************

bool Sudoku::solved() const
{
    for ( size_t i( 0 ); i != values_.size(); ++i )
    {
        if ( ! values_[i].solved() )
            return false;
    }
    return ( ! this->there_are_contradictions() );
}

// ********************************************************************************

bool Sudoku::there_are_contradictions() const
{
    for ( size_t i( 0 ); i != Sudoku::number_of_slices(); ++i )
    {
        OneSudokuSlice slice = this->slice( i );
        SetOfNumbers values = slice.collect_solved_numbers();
        if ( values.contains_duplicates() )
            return true;
    }
    return false;
}

// ********************************************************************************

// Returns the three slices that contain this quare
std::vector< size_t > Sudoku::square2slices( const size_t square_index )
{
    std::vector< size_t > result;
    for ( size_t i( 0 ); i != Sudoku::number_of_slices(); ++i )
    {
        for ( size_t j( 0 ); j != 9; ++j )
        {
            if ( sudoku_slice_mapping[i][j] == square_index )
            {
                result.push_back( i );
                break;
            }
        }
    }
    return result;
}

// ********************************************************************************

void Sudoku::show() const
{
    for ( size_t i( 0 ); i != 9; ++i )
    {
        std::cout << "  -------------------------------------------------------------------------" << std::endl;
        for ( size_t iRow( 0 ); iRow != 3; ++iRow )
        {
            std::cout << "  |  ";
            for ( size_t j( 0 ); j != 9; ++j )
            {
                if ( values_[ i*9 + j ].contains( iRow*3 + 1 ) )
                    std::cout << iRow*3 + 1;
                else
                    std::cout << " ";
                if ( values_[ i*9 + j ].contains( iRow*3 + 2 ) )
                    std::cout << iRow*3 + 2;
                else
                    std::cout << " ";
                if ( values_[ i*9 + j ].contains( iRow*3 + 3 ) )
                    std::cout << iRow*3 + 3;
                else
                    std::cout << " ";
                std::cout << "  |  ";
            }
            std::cout << std::endl;
        }
    }
    std::cout << "  -------------------------------------------------------------------------" << std::endl;
    this->show_statistics();
}

// ********************************************************************************

void Sudoku::show_statistics() const
{
    size_t nsolved( 0 );
    size_t npossibilities( 0 );
    for ( size_t i( 0 ); i != values_.size(); ++i )
    {
        if ( values_[i].solved() )
            ++nsolved;
        else
            npossibilities += values_[i].size();
    }
    std::cout << "# Solved = " << nsolved << std::endl;
    std::cout << "# Possibilities = " << npossibilities << std::endl;
}

// ********************************************************************************

OneSudokuSlice::SudokuSliceType Sudoku::slice_id2type( const size_t i ) const
{
    if ( i < 9 )
        return OneSudokuSlice::ROW;
    if ( i < 18 )
        return OneSudokuSlice::COLUMN;
    return OneSudokuSlice::SQUARE;
}

// ********************************************************************************

std::vector< OneSudokuSquare > Sudoku::get_slice_values( const size_t j ) const
{
    std::vector< OneSudokuSquare > result;
    for ( size_t i( 0 ); i != 9; ++i )
        result.push_back( values_[ sudoku_slice_mapping[j][i] ] );
    return result;
}

// ********************************************************************************

void Sudoku::initialise_trans_square_dependency_mappings()
{
    return;
    size_t index_1( 0 );
    size_t index_2( 0 );
    // Rows
    for ( size_t i( 0 ); i != 81; i = i+3 )
    {
        // index 0 is row, index 1 is column and index 2 is square
        std::vector< size_t > slices = square2slices( i );
        // determining_values is square minus { i, i+1, i+2 }
        // values_to_be_changed is row minus  { i, i+1, i+2 }
        index_2 = 0;
        for ( size_t j( 0 ); j != 9; ++j )
        {
            size_t value = sudoku_slice_mapping[ slices[2] ][j];
            if ( ( value != i     ) &&
                 ( value != i + 1 ) &&
                 ( value != i + 2 ) )
            {
                sudoku_determining_values_mapping[index_1][index_2] = value;
                ++index_2;
            }
        }
        std::vector< size_t > values_to_be_changed;
        index_2 = 0;
        for ( size_t j( 0 ); j != 9; ++j )
        {
            size_t value = sudoku_slice_mapping[ slices[0] ][j];
            if ( ( value != i     ) &&
                 ( value != i + 1 ) &&
                 ( value != i + 2 ) )
            {
                sudoku_values_to_be_changed_mapping[index_1][index_2] = value;
                ++index_2;
            }
        }
        ++index_1;
    }
    // Columns
    for ( size_t i( 0 ); i != 81; ++i )
    {
        if ( i > 62 )
            continue;
        if ( ( i > 35 ) && ( i < 54 ) )
            continue;
        if ( ( i >  8 ) && ( i < 27 ) )
            continue;
        // index 0 is row, index 1 is column and index 2 is square
        std::vector< size_t > slices = square2slices( i );
        // determining_values is square minus { i, i+1, i+2 }
        // values_to_be_changed is row minus  { i, i+1, i+2 }
        index_2 = 0;
        for ( size_t j( 0 ); j != 9; ++j )
        {
            size_t value = sudoku_slice_mapping[ slices[2] ][j];
            if ( ( value != i      ) &&
                 ( value != i +  9 ) &&
                 ( value != i + 18 ) )
            {
                sudoku_determining_values_mapping[index_1][index_2] = value;
                ++index_2;
            }
        }
        std::vector< size_t > values_to_be_changed;
        index_2 = 0;
        for ( size_t j( 0 ); j != 9; ++j )
        {
            size_t value = sudoku_slice_mapping[ slices[1] ][j];
            if ( ( value != i      ) &&
                 ( value != i +  9 ) &&
                 ( value != i + 18 ) )
            {
                sudoku_values_to_be_changed_mapping[index_1][index_2] = value;
                ++index_2;
            }
        }
        ++index_1;
    }
    for ( size_t i( 0 ); i != 54; ++i )
    {
        std::cout << "                       { ";
        for ( size_t j( 0 ); j != 6; ++j )
        {
            std::cout << " " << sudoku_determining_values_mapping[i][j];
            if ( j != 5 )
                std::cout << ",";
        }
        std::cout << " }";
        if ( i != 53 )
            std::cout << ",";
        std::cout << std::endl;
    }
    std::cout << std::endl;
    for ( size_t i( 0 ); i != 54; ++i )
    {
        std::cout << "                       { ";
        for ( size_t j( 0 ); j != 6; ++j )
        {
            std::cout << " " << sudoku_values_to_be_changed_mapping[i][j];
            if ( j != 5 )
                std::cout << ",";
        }
        std::cout << " }";
        if ( i != 53 )
            std::cout << ",";
        std::cout << std::endl;
    }
}

// ********************************************************************************

