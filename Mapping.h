#ifndef MAPPING_H
#define MAPPING_H

/* *********************************************
Copyright (c) 2013-2025, Cornelis Jan (Jacco) van de Streek
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

#include <cstddef> // For definition of size_t
#include <vector>

/*
    Must be a one-to-one mapping that maps e.g. 0, 1, 2 to 2, 0, 1.
*/
class Mapping
{
public:

    explicit Mapping( const size_t nvalues = 0 );

    explicit Mapping( const std::vector< size_t > & values );

    size_t size() const { return mapping_.size(); }

    void swap( const size_t i, const size_t j );
    
    // Increments size by one. Can of course only push back a mapping to itself.
    void push_back();

    // Note that this is all const, i.e. it is read-only.
    size_t operator[]( const size_t i ) const;

    std::vector< size_t >::iterator begin() { return mapping_.begin(); }
    std::vector< size_t >::iterator end() { return mapping_.end(); }
    std::vector< size_t >::const_iterator begin() const { return mapping_.begin(); }
    std::vector< size_t >::const_iterator end() const { return mapping_.end(); }

    void invert();

private:
    std::vector< size_t > mapping_;
    
    void check_consistency() const;
};

#endif // MAPPING_H

