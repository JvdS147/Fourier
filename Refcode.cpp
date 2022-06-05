/* *********************************************
Copyright (c) 2013-2022, Cornelis Jan (Jacco) van de Streek
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

#include "Refcode.h"
#include "StringFunctions.h"

#include <stdexcept>

// ********************************************************************************

Refcode::Refcode(): refcode_("XXXXXX01")
{
}

// ********************************************************************************

// Throws if format incorrect
Refcode::Refcode( const std::string & refcode )
{
    if ( ! format_is_correct( refcode ) )
        throw std::runtime_error( "Refcode::Refcode(): incorrect format: " + refcode );
    refcode_ = to_upper( refcode );
}

// ********************************************************************************

// Strips digits and returns first six alphabetical characters
std::string Refcode::family() const
{
    return refcode_.substr( 0, 6 );
}

// ********************************************************************************

// Convenience function
// Capitalisation is ignored
// XXXXXX00 is allowed
bool Refcode::format_is_correct( std::string refcode )
{
    // Length must be six or eight
    if ( ( refcode.length() != 6 ) &&
         ( refcode.length() != 8 ) )
        return false;
    // First six characters must be alphabetical characters
    for ( size_t i( 0 ); i != 6; ++i )
    {
        if ( ! std::isalpha( refcode[i] ) )
            return false;
    }
    if ( refcode.length() == 6 )
        return true;
    // Last two must be digits
    if ( ! std::isdigit( refcode[6] ) )
        return false;
    if ( ! std::isdigit( refcode[7] ) )
        return false;
    return true;
}

// ********************************************************************************

