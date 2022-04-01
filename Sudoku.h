#ifndef SUDOKU_H
#define SUDOKU_H

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

//class OneSudokuSlice;
class TransSquareDependency;

#include "OneSudokuSquare.h"
#include "CyclicInteger.h"
#include "OneSudokuSlice.h"

#include <string>
#include <vector>

/*

*/
class Sudoku
{
public:

    // Default constructor
    Sudoku();
    
    Sudoku( const std::vector< std::string > & rows );

    static size_t number_of_slices() { return 27; }
    
    OneSudokuSlice slice( const size_t i ) const;

    // 0 based
    OneSudokuSlice row( const size_t i ) const;

    // 0 based
    OneSudokuSlice column( const size_t i ) const;

    // 0 based
    OneSudokuSlice block( const size_t i ) const;

    size_t nsquares() const { return values_.size(); }
    
    OneSudokuSquare square( const size_t index ) const { return values_[index]; }
    
    // Could be called set_square(), I guess    
    bool update_square( const size_t i, OneSudokuSquare square );

    // Does not return solved slices
    OneSudokuSlice next_slice();

    // Consistency checking is performed
    void update_slice( const OneSudokuSlice & slice );
    
    // i is the index of the square
    bool unset( const size_t i, const size_t value );
    
    size_t number_of_trans_square_dependencies() const { return 54; }

    TransSquareDependency next_trans_square_dependency();

    // To be called immediately after next_trans_square_dependency().
    void update_trans_square_dependency( const TransSquareDependency & tsd );

    bool solved() const;
    
    bool there_are_contradictions() const;
    
    // Returns the three slices that contain this quare
    // index 0 is ROW, index 1 is COLUMN and index 2 is SQUARE
    static std::vector< size_t > square2slices( const size_t square_index );

    void show() const;

    void show_statistics() const;
    
    static void initialise_trans_square_dependency_mappings();

private:
    
    std::vector< OneSudokuSquare > values_;
    CyclicInteger current_slice_;
    CyclicInteger current_trans_square_dependency_;
    
    OneSudokuSlice::SudokuSliceType slice_id2type( const size_t i ) const;
    std::vector< OneSudokuSquare > get_slice_values( const size_t i ) const;
};

#endif // SUDOKU_H

