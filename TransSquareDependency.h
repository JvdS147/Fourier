#ifndef TRANSSQUAREDEPENDENCY_H
#define TRANSSQUAREDEPENDENCY_H

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

#include "OneSudokuSlice.h"

/*

| D | D | D |   |   |   |   |   |   |
| D | D | D |   |   |   |   |   |   |
|   |   |   | C | C | C | C | C | C |

Values not possible in the squares marked D, are not allowed in the squares marked C.

So we take all the possible values in the squares marked D, then take the set 1 through 9 minus those values
and eliminate those values from the squares marked C.

*/
class TransSquareDependency
{
public:

    // Default constructor
//    TransSquareDependency();

    TransSquareDependency( const size_t id, const std::vector< OneSudokuSquare > & determining_values, const std::vector< OneSudokuSquare > & values_to_be_changed );

    size_t id() const { return id_; }

    OneSudokuSlice determining_values() const { return determining_values_; }
    OneSudokuSlice values_to_be_changed() const { return values_to_be_changed_; }
    
    void update_values_to_be_changed( const OneSudokuSlice & values_to_be_changed );

    size_t determining_values_size() const { return 6; }
    size_t values_to_be_changed_size() const { return 6; }
    
   // bool solved() const;

    void show() const;

private:
    size_t id_;
    OneSudokuSlice determining_values_;
    OneSudokuSlice values_to_be_changed_;
};

#endif // TRANSSQUAREDEPENDENCY_H

