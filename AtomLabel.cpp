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

#include "AtomLabel.h"
#include "StringFunctions.h"
#include "Utilities.h"

#include <stdexcept>

// ********************************************************************************

AtomLabel::AtomLabel( const std::string & label ):
label_(label),
index_(0),
symmetry_operator_index_(0),
u_(0),
v_(0),
w_(0)
{
    if ( label.empty() )
        throw std::runtime_error( "AtomLabel::AtomLabel(): Error: label is empty." );
    std::vector< std::string > words = split( label, '_' );
    if ( words.size() == 5 )
    {
        symmetry_operator_index_ = string2integer( words[1] );
        u_ = string2integer( words[2] );
        v_ = string2integer( words[3] );
        w_ = string2integer( words[4] );
    }
    else if ( words.size() == 4 )
    {
        u_ = string2integer( words[1] );
        v_ = string2integer( words[2] );
        w_ = string2integer( words[3] );
    }
    else if ( words.size() != 1 )
        throw std::runtime_error( "AtomLabel::AtomLabel(): Error." );
    std::string label_2 = words[0];
    if ( ! is_letter( label_2[0] ) )
        throw std::runtime_error( "AtomLabel::AtomLabel(): Error." );
    element_ = label_2[0];
    label_2 = label_2.substr( 1 );
    if ( label_2.empty() )
        return;
    if ( is_letter( label_2[0] ) )
    {
        element_ += label_2[0];
        label_2 = label_2.substr( 1 );
        if ( label_2.empty() )
            return;
    }
    if ( ! is_digit( label_2[0] ) ) // "Unh" etc. exist, but encountering them would still be an error.
        throw std::runtime_error( "AtomLabel::AtomLabel(): Error." );
    size_t iPos( 1 );
    while ( ( iPos < label_2.length() ) && is_digit( label_2[iPos] ) )
    {
        ++iPos;
    }
    index_ = string2integer( label_2.substr( 0, iPos ) );
    subindex_ = label_2.substr( iPos );
}

// ********************************************************************************

