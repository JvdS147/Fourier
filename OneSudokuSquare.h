#ifndef ONESUDOKUSQUARE_H
#define ONESUDOKUSQUARE_H

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

#include "SetOfNumbers.h"

#include <cstddef> // For definition of size_t
#include <vector>

/*
Cannot be empty.
No duplicates.
Only 1 through 9.
*/
class OneSudokuSquare
{
public:

    // Default constructor, initialises to 1 though 9
    OneSudokuSquare();

    // If the value is 0, initialises to 1 though 9
    explicit OneSudokuSquare( const size_t value );

    explicit OneSudokuSquare( const SetOfNumbers & values );
    
    SetOfNumbers values() const { return values_; }
    
    size_t value( const size_t i ) const;
    
    // Returns its value when solved, throws otherwise
    size_t value() const;
    
    // Returns true if something has changed
    // The value is added to the other values, i.e. the bit for the value is set
    bool set( const size_t value );

//    // Returns true if something has changed
//    // The value is added to the other values, i.e. the bit for the value is set
//    bool set( const OneSudokuSquare & sudoku_square );
    // Returns true if something has changed
    bool unset( const size_t value );
    
    // The version with the OneSudokuSquare argument is less versatile (it does not allow for the vector to be empty)
    // and the two can probably be merged and replaced by the std::vector< int > version.
    // Returns true if something has changed
    bool unset( const OneSudokuSquare & sudoku_square );
    
    // Returns true if something has changed
    bool unset( const SetOfNumbers & values );

    bool contains( const size_t value ) const;

    size_t size() const;
    
    bool solved() const;

    void show() const;

    bool operator==( const OneSudokuSquare & rhs ) const { return ( this->values_ == rhs.values_ ); }

private:
    // Should probably store a SetOfNumbers
    SetOfNumbers values_; // Always stored in ascending order
};

// @@ Should have been a member function?
// @@ Should simply have used std::vector< size_t > ?
OneSudokuSquare merge( const OneSudokuSquare & lhs, const OneSudokuSquare & rhs );

// Duplicates are removed.
//std::vector< size_t > merge( const std::vector< size_t > & lhs, const std::vector< size_t > & rhs );

#endif // ONESUDOKUSQUARE_H

