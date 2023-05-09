#ifndef REFCODELIST_H
#define REFCODELIST_H

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

class FileName;

#include "Mapping.h"
#include "Refcode.h"

#include <vector>

// @@ There is currently no sort() function.
class RefcodeList
{
public:

    RefcodeList();

    RefcodeList( const std::vector< std::string > & values );

    void push_back( const Refcode & refcode );

    void reserve( const size_t nvalues );
    size_t size() const { return refcodes_.size(); }

    bool contains( const Refcode & refcode ) const;

    // Returns size() if not found.
    size_t index( const Refcode & refcode ) const;

    // The index is zero-based.
    // We don't actually sort the lists, but create a sorted map.
    Refcode refcode( const size_t i ) const { return refcodes_[ sorted_map_[ i ] ];   }

    void set_refcode( const size_t i, const Refcode & refcode ) { refcodes_[ sorted_map_[ i ] ] = refcode; }

    // Converts each entry to its refcode family, then removes duplicates.
    void convert_to_unique_families();

    void read_gcd( const FileName & file_name );

    // For debugging.
    void show() const;

    // The extension should be .gcd .
    void save( const FileName & file_name ) const;

private:
    std::vector< Refcode > refcodes_;
    // We don't actually sort the list, but create a sorted map.
    Mapping sorted_map_;

};

#endif // REFCODELIST_H

