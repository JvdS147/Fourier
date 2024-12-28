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

#include "StringConversions.h"
#include "Angle.h"
#include "Matrix3D.h"
#include "MillerIndices.h"
#include "StringFunctions.h"
#include "Utilities.h"
#include "Vector3D.h"

#include <stdexcept>

// ********************************************************************************

Angle Angle_from_string( const std::string & input )
{
    return Angle( string2double( input ), Angle::DEGREES );
}

// ********************************************************************************

Matrix3D Matrix3D_from_string( std::string input )
{
    input = remove( input, ' ' );
    input = remove( input, '[' );
    input = remove( input, ']' );
    std::vector< std::string > words = split( input, ',' );
    if ( words.size() != 9 )
        throw std::runtime_error( "Matrix3D_from_string(): error: incorrect format. |" + input + "|" );
    return Matrix3D( string2double( words[0] ), string2double( words[1] ), string2double( words[2] ),
                     string2double( words[3] ), string2double( words[4] ), string2double( words[5] ),
                     string2double( words[6] ), string2double( words[7] ), string2double( words[8] ) );
}

// ********************************************************************************

MillerIndices MillerIndices_from_string( std::string input )
{
    input = remove( input, ' ' );
    input = remove( input, '(' );
    input = remove( input, ')' );
    std::vector< std::string > words = split( input, ',' );
    if ( words.size() != 3 )
        throw std::runtime_error( "MillerIndices_from_string(): error: incorrect format. |" + input + "|" );
    return MillerIndices( string2integer( words[0] ), string2integer( words[1] ), string2integer( words[2] ) );
}

// ********************************************************************************

Vector3D Vector3D_from_string( std::string input )
{
    input = remove( input, ' ' );
    input = remove( input, '[' );
    input = remove( input, ']' );
    std::vector< std::string > words = split( input, ',' );
    if ( words.size() != 3 )
        throw std::runtime_error( "Vector3D_from_string(): error: incorrect format. |" + input + "|" );
    return Vector3D( string2double( words[0] ), string2double( words[1] ), string2double( words[2] ) );
}

// ********************************************************************************

std::string to_string( const Angle input )
{
    return double2string( input.value_in_degrees() );
}

// ********************************************************************************

std::string to_string( const Matrix3D & input )
{
    return "[ [ " + double2string( input.value( 0, 0 ) ) + ", " + double2string( input.value( 0, 1 ) ) + ", " + double2string( input.value( 0, 2 ) ) + " ], " +
            " [ " + double2string( input.value( 1, 0 ) ) + ", " + double2string( input.value( 1, 1 ) ) + ", " + double2string( input.value( 1, 2 ) ) + " ], " +
            " [ " + double2string( input.value( 2, 0 ) ) + ", " + double2string( input.value( 2, 1 ) ) + ", " + double2string( input.value( 2, 2 ) ) + " ] ]";
}

// ********************************************************************************

std::string to_string( const Vector3D & input )
{
    return "[ " + double2string( input.value( 0 ) ) + ", " + double2string( input.value( 1 ) ) + ", " + double2string( input.value( 2 ) ) + " ]";
}

// ********************************************************************************

std::string to_string( const MillerIndices & input )
{
    return "( " + int2string( input.h() ) + ", " + int2string( input.k() ) + ", " + int2string( input.l() ) + " )";
}

// ********************************************************************************

