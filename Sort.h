#ifndef SORT_H
#define SORT_H

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

#include "Mapping.h"

#include <algorithm>
#include <cstddef> // For definition of size_t
#include <vector>

// Uses std::stable_sort().
// The original data are not sorted, instead a sorted map is returned. This has several advantages:
// 1. The map can be used to address multiple lists in parallel
// 2. No copying of large objects
// 3. The original data is not changed and can be const
// at the expense of some slight overhead in the form of a std::vector< size_t >
// and the requirement that all elements in the original data structure must now be addressed as:
// value = values[ sorted_map[i] ]; and values[ sorted_map[i] ] = value;

template <class T>
class Compare
{
public:

    Compare( const std::vector< T > & values, const bool reverse_order = false ): values_(values),reverse_order_(reverse_order) {}

    bool operator()( const size_t lhs, const size_t rhs ) { if ( reverse_order_ )
                                                                return ( values_[rhs] < values_[lhs] );
                                                            else
                                                                return ( values_[lhs] < values_[rhs] );
                                                          }

private:
    const std::vector< T > & values_;
    bool reverse_order_;
};

// ********************************************************************************

template <class T>
Mapping sort( const std::vector< T > & values, const bool reverse_order = false )
{
    // We don't actually sort the list, but create a sorted map
    // We use std::stable_sort() with a functor
    Mapping sorted_map( values.size() );
    std::stable_sort( sorted_map.begin(), sorted_map.end(), Compare<T>( values, reverse_order ) );
    return sorted_map;
}

// ********************************************************************************

#endif // SORT_H

