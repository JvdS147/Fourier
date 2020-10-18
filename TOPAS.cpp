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

#include "TOPAS.h"

#include "Angle.h"
#include "CrystalLattice.h"
#include "TextFileReader_2.h"
#include "Utilities.h"

#include <iostream> // For debugging
#include <stdexcept>

// ********************************************************************************

std::string insert_at_sign( std::string input )
{
//     CS_L( , 10000.00000_LIMIT_MIN_0.3)
//   Zero_Error( ,-0.01424)
    size_t iPos = input.find( "(" );
    if ( iPos == std::string::npos )
        throw std::runtime_error( "insert_at_sign(): opening parenthesis not found: >" + input + "<" );
    return input.substr( 0, iPos+1 ) + "@" + input.substr( iPos + 1 );
}

// ********************************************************************************

// Not allowed: "prm!a 1", "prm ! a 1", it must be "prm !a 1" or "prm a 1"
double read_keyword( const std::string & keyword, TextFileReader_2 & input_file )
{
    size_t iLine( 0 );
    // First try to find "!keyword", then "keyword"
    do
    {
        iLine = input_file.find_whole_word( "!" + keyword, iLine+1 );
        if ( iLine == std::string::npos )
            break;
        // Remove the bit that has been commented out
        std::string line = input_file.line( iLine );
        size_t iPos = line.find( "'" );
        line = line.substr( 0, iPos );
        iPos = line.find( "!" + keyword );
        if ( iPos == std::string::npos )
            continue;
        std::vector< std::string > words = split( line );
        if ( words.size() < 3 )
            continue;
        if ( ( words[0] == "prm" ) && ( words[1] == "!" + keyword ) )
            return TOPASstring2double( words[2] );
    }
    while ( true );
    iLine = 0;
    do
    {
        iLine = input_file.find_whole_word( keyword, iLine+1 );
        if ( iLine == std::string::npos )
            throw std::runtime_error( "read_keyword(): keyword not found." );
        // Remove the bit that has been commented out
        std::string line = input_file.line( iLine );
        size_t iPos = line.find( "'" );
        line = line.substr( 0, iPos );
        iPos = line.find( keyword );
        if ( iPos == std::string::npos )
            continue;
        std::vector< std::string > words = split( line );
        if ( words.size() < 3 )
            continue;
        if ( ( words[0] == "prm" ) && ( words[1] == keyword ) )
            return TOPASstring2double( words[2] );
    }
    while ( true );
}

// ********************************************************************************

CrystalLattice read_lattice_parameters( TextFileReader_2 & input_file )
{
    bool found( false );
    size_t iLine( 0 );
    std::vector< std::string > words;
    while ( ( ! found ) && ( iLine != input_file.size() ) )
    {
        words = split_2( input_file.line( iLine ) );
        if ( words.size() > 1 )
        {
            if ( words[0] == "a" )
            {
                if ( words[1] == "@" )
                {
                    if ( words.size() > 2 )
                    {
                        try
                        {
                            /* double dummy = */ TOPASstring2double( words[2] );
                            found = true;
                        }
                        catch ( std::exception & e ) {}
                    }
                }
                else
                {
                    try
                    {
                        /* double dummy = */ TOPASstring2double( words[1] );
                        found = true;
                    }
                    catch ( std::exception & e ) {}
                }
            }
        }
        if ( ! found )
            ++iLine;
    }
    // When we are here, either we are at the end of the file or we are at the line "a @ 8.743".
    if ( ! found )
        throw std::runtime_error( "read_lattice_parameters(): could not find lattice parameter a." );
    // Assume that the following lines are a, b, c, al, be and ga
    double a;
    double b;
    double c;
    double al;
    double be;
    double ga;
    found = false;
    words = split_2( input_file.line( iLine ) );
    if ( words.size() > 1 )
    {
        if ( words[0] == "a" )
        {
            if ( words[1] == "@" )
            {
                if ( words.size() > 2 )
                {
                    a = TOPASstring2double( words[2] );
                    found = true;
                }
            }
            else
            {
                a = TOPASstring2double( words[1] );
                found = true;
            }
        }
    }
    if ( ! found )
        throw std::runtime_error( "read_lattice_parameters(): could not find lattice parameter a." );
    ++iLine;
    found = false;
    words = split_2( input_file.line( iLine ) );
    if ( words.size() > 1 )
    {
        if ( words[0] == "b" )
        {
            if ( words[1] == "@" )
            {
                if ( words.size() > 2 )
                {
                    b = TOPASstring2double( words[2] );
                    found = true;
                }
            }
            else
            {
                b = TOPASstring2double( words[1] );
                found = true;
            }
        }
    }
    if ( ! found )
        throw std::runtime_error( "read_lattice_parameters(): could not find lattice parameter b." );
    ++iLine;
    found = false;
    words = split_2( input_file.line( iLine ) );
    if ( words.size() > 1 )
    {
        if ( words[0] == "c" )
        {
            if ( words[1] == "@" )
            {
                if ( words.size() > 2 )
                {
                    c = TOPASstring2double( words[2] );
                    found = true;
                }
            }
            else
            {
                c = TOPASstring2double( words[1] );
                found = true;
            }
        }
    }
    if ( ! found )
        throw std::runtime_error( "read_lattice_parameters(): could not find lattice parameter c." );
    ++iLine;
    found = false;
    words = split_2( input_file.line( iLine ) );
    if ( words.size() > 1 )
    {
        if ( words[0] == "al" )
        {
            if ( words[1] == "@" )
            {
                if ( words.size() > 2 )
                {
                    al = TOPASstring2double( words[2] );
                    found = true;
                }
            }
            else
            {
                al = TOPASstring2double( words[1] );
                found = true;
            }
        }
    }
    if ( ! found )
        throw std::runtime_error( "read_lattice_parameters(): could not find lattice parameter al." );
    ++iLine;
    found = false;
    words = split_2( input_file.line( iLine ) );
    if ( words.size() > 1 )
    {
        if ( words[0] == "be" )
        {
            if ( words[1] == "@" )
            {
                if ( words.size() > 2 )
                {
                    be = TOPASstring2double( words[2] );
                    found = true;
                }
            }
            else
            {
                be = TOPASstring2double( words[1] );
                found = true;
            }
        }
    }
    if ( ! found )
        throw std::runtime_error( "read_lattice_parameters(): could not find lattice parameter be." );
    ++iLine;
    found = false;
    words = split_2( input_file.line( iLine ) );
    if ( words.size() > 1 )
    {
        if ( words[0] == "ga" )
        {
            if ( words[1] == "@" )
            {
                if ( words.size() > 2 )
                {
                    ga = TOPASstring2double( words[2] );
                    found = true;
                }
            }
            else
            {
                ga = TOPASstring2double( words[1] );
                found = true;
            }
        }
    }
    if ( ! found )
        throw std::runtime_error( "read_lattice_parameters(): could not find lattice parameter ga." );
    return CrystalLattice( a, b, c,
                           Angle::from_degrees( al ),
                           Angle::from_degrees( be ),
                           Angle::from_degrees( ga ) );
}

// ********************************************************************************

// Understands "118.34201`_0.68292", "118.34201`" and "118.34201".
double TOPASstring2double( std::string input )
{
    size_t iPos = input.find_first_of( "_" );
    // Check that the string is formatted properly
    if ( iPos != std::string::npos )
    {
        double dummy = string2double_2( input.substr( iPos+1 ), true );
        input.erase( iPos );
        if ( input.size() == 0 )
            throw std::runtime_error( "TOPASstring2double(): no number before ESD :  >" + input + "<" );
        if ( input[input.length()-1] != '`' )
            throw std::runtime_error( "TOPASstring2double(): no '' after number with ESD :  >" + input + "<" );
    }
    if ( input.size() == 0 )
        throw std::runtime_error( "TOPASstring2double(): empty number : " + input );
    if ( input[input.length()-1] == '`' )
        input.erase( input.length()-1 );
    return string2double_2( input, true );
}

// ********************************************************************************

