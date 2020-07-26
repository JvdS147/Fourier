#ifndef BAGOFNUMBERS_H
#define BAGOFNUMBERS_H

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

#include "RandomNumberGenerator.h"

#include <cstddef> // For definition of size_t
#include <vector>

/*
 * This class represents a bag of numbers.
 *
 * Numbers are not necessarily unique.
 *
 * Helpful for randomising sequences.
 *
 */
class BagOfNumbers
{
public:

    // Default constructor, creates an empty bag,
    explicit BagOfNumbers( const int idum );

    // In keeping with C++ convention: zero-based.
    // Fills the bag with the numbers 0, 1, ..., nnumbers-1.
    explicit BagOfNumbers( const size_t nnumbers, const int idum  );
    
    // Does nothing if value currently not in the bag.
    // If multiple occurrences present, only removes one.
    void remove( const int value );

    // If value currently already in the bag, adds it again
    void add( const int value );

    // Returns one of the numbers at random and removes it from the bag.
    // Throws if the bag is empty.
    size_t draw();
    
    // Returns the number of numbers in the bag.
    size_t size() const { return values_.size(); }
    
    // Returns one of the numbers at random, it is NOT removed from the bag
    // Throws if the bag is empty.
//    size_t draw_with_replace() const;
    size_t draw_with_replace();

private:
    std::vector< int > values_;
    RandomNumberGenerator_integer RNG_int_;
    
    void remove_position( const size_t iPos );
};

#endif // BAGOFNUMBERS_H

