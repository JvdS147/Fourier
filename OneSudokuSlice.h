#ifndef ONESUDOKUSLICE_H
#define ONESUDOKUSLICE_H

/* *********************************************
Copyright (c) 2013-2021, Cornelis Jan (Jacco) van de Streek
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

#include <cstddef> // For definition of size_t
#include <vector>
#include <string>

/*
    Each square is initialised to all nine values. In other words, there is no marker for "unknown", neither is there a marker for "value is known",
    every square simply has between 1 and 9 values.
*/
class OneSudokuSlice
{
public:

    enum SudokuSliceType { OTHER, ROW, COLUMN, SQUARE };

    // Default constructor
 //   OneSudokuSlice();

//    OneSudokuSlice( const size_t id, const std::string & one_row );

    OneSudokuSlice( const size_t id, const SudokuSliceType type, const std::vector< OneSudokuSquare > & values );
                    
    size_t id() const { return id_; }
    
    SudokuSliceType type() const { return type_; }

    OneSudokuSquare square( const size_t i ) const;

    void set_square( const size_t i, const OneSudokuSquare & value );

    size_t size() const { return values_.size(); }
    
    bool solved() const;

    SetOfNumbers collect_solved_numbers() const;
    
    // Only makes sense if size() == 6, otherwise it just returns 1 through 9
    SetOfNumbers collect_all_possible_numbers() const;

    std::vector< size_t > indices_of_unsolved_squares() const;

    std::vector< size_t > indices_of_solved_squares() const;

    void show() const;

private:
    size_t id_;
    SudokuSliceType type_;
    std::vector< OneSudokuSquare > values_;
};

#endif // ONESUDOKUSLICE_H

