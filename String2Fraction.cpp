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

#include "String2Fraction.h"
#include "StringFunctions.h"
#include "Utilities.h"

#include <stdexcept>

// ********************************************************************************

// Expects "1/2", can cope with "0.5", "-1-1/-2"
Fraction string2Fraction( std::string fraction_str )
{
    fraction_str = strip( fraction_str );
    if ( fraction_str.length() == 0 )
        throw std::runtime_error( "string2Fraction(): empty string" );
    if ( fraction_str[0] == '+' )
        fraction_str = fraction_str.substr( 1 );
    if ( fraction_str.length() == 0 )
        throw std::runtime_error( "string2Fraction(): empty string" );
    size_t iPos = fraction_str.find( "." );
    if ( iPos != std::string::npos )
        return double2fraction( string2double( fraction_str ), Fraction( 1, 8 ) );
    iPos = fraction_str.find( "/" );
    if ( iPos == 0 )
        throw std::runtime_error( "string2Fraction(): incorrect format: |" + fraction_str + "|" );
    if ( iPos == fraction_str.length()-1 )
        throw std::runtime_error( "string2Fraction(): incorrect format: |" + fraction_str + "|" );
    if ( iPos == std::string::npos )
        return Fraction( string2integer( fraction_str ) );
    int denominator = string2integer( strip( fraction_str.substr( iPos+1 ) ) );
    fraction_str = strip( fraction_str.substr( 0, iPos ) );
    size_t nplus  = count_characters( fraction_str, '+' );
    size_t nminus = count_characters( fraction_str, '-' );
    if ( nplus > 1 )
        throw std::runtime_error( "string2Fraction(): incorrect format: |" + fraction_str + "|" );
    if ( ( nplus + nminus ) > 2 )
        throw std::runtime_error( "string2Fraction(): incorrect format: |" + fraction_str + "|" );
    // "1", "1+1", "-1+1", "1-1", "-1", "-1-1"
    if ( ( nplus + nminus ) == 0 )
        return Fraction( string2integer( fraction_str ), denominator );
    if ( nplus == 1 )
    {
        iPos = fraction_str.find( "+" );
        return Fraction( string2integer( strip( fraction_str.substr( 0, iPos ) ) ), string2integer( strip( fraction_str.substr( iPos ) ) ), denominator );
    }
    // "1-1", "-1", "-1-1" left
    iPos = fraction_str.find_last_of( "-" );
    if ( iPos == 0 )
        return Fraction( string2integer( fraction_str ), denominator );
    else
        return Fraction( string2integer( strip( fraction_str.substr( 0, iPos ) ) ), string2integer( strip( fraction_str.substr( iPos ) ) ), denominator );
}

// ********************************************************************************

