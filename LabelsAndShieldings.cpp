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

#include "LabelsAndShieldings.h"
#include "TextFileWriter.h"
#include "Utilities.h"

#include <algorithm>

namespace {

class Compare
{
public:

    Compare( const std::vector< Element > & elements, const std::vector< double > & shieldings ): elements_(elements), shieldings_(shieldings)
    {
    }

    bool operator()( const size_t lhs, const size_t rhs )
    {
        if ( elements_[rhs] < elements_[lhs] )
            return true;
        if ( elements_[lhs] < elements_[rhs] )
            return false;
        return shieldings_[rhs] < shieldings_[lhs];
    }

private:
    const std::vector< Element > & elements_;
    const std::vector< double > & shieldings_;
};

} // namespace

// ********************************************************************************

void LabelsAndShieldings::push_back( const std::string & label, const double shielding )
{
    labels_.push_back( label );
    std::string element_string;
    if ( label.size() == 1 )
        element_string = label;
    else
    {
        if ( isalpha(label[1]) )
            element_string = label.substr( 0, 2 );
        else
            element_string = label.substr( 0, 1 );
    }
    elements_.push_back( Element( element_string ) );
    shieldings_.push_back( shielding );
    sorted_map_.push_back( size() - 1 );
}

// ********************************************************************************

void LabelsAndShieldings::sort() const
{
    // We don't actually sort the lists, but create a sorted map
    // We use std::sort() with a functor
    std::sort( sorted_map_.begin(), sorted_map_.end(), Compare( elements_, shieldings_ ) );
}

// ********************************************************************************

void LabelsAndShieldings::save( const FileName & output_file_name ) const
{
    sort();
    TextFileWriter text_file_writer( output_file_name );
    for ( size_t i( 0 ); i != size(); ++i )
        text_file_writer.write_line( label(i) + " " + double2string( shielding(i) ) );
}

// ********************************************************************************
