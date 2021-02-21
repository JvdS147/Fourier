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

#include "GenerateCombinations.h"

#include <algorithm>

// ********************************************************************************

GenerateCombinations::GenerateCombinations( const std::vector< size_t > & values, const size_t k )
{
    another_one_is_available_ = true;
    values_ = values;
    k_ = k;
    bitmask_ = std::vector< size_t >( k_, 1 ); // k leading 1's
    bitmask_.resize( values_.size(), 0 ); // N-k trailing 0's
}

// ********************************************************************************

bool GenerateCombinations::next_combination( std::vector< size_t > & result ) const
{
    if ( ! another_one_is_available_ )
        return false;
    result.clear();
    for ( size_t i( 0 ); i != bitmask_.size(); ++i )
    {
        if ( bitmask_[i] )
            result.push_back( values_[i] );
    }
    another_one_is_available_ = std::prev_permutation( bitmask_.begin(), bitmask_.end() );
    return true;
}

// ********************************************************************************

