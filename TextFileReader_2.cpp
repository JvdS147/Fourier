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

#include "TextFileReader_2.h"
#include "FileName.h"
#include "Utilities.h"

#include <stdexcept>
#include <iostream>

// ********************************************************************************

TextFileReader_2::TextFileReader_2( const FileName & file_name )
{
    read_file( file_name );
}

// ********************************************************************************

void TextFileReader_2::read_file( const FileName & file_name )
{
    std::ifstream input_file( file_name.full_name().c_str() );
    if ( ! input_file )
       throw std::runtime_error( std::string( "TextFileReader::read_file(): Could not open file " ) + file_name.full_name() );
    lines_.clear();
    std::string line;
    do
    {
        if ( ! getline( input_file, line ) )
        {
            input_file.close();
            return;
        }
        // remove \r
        line = remove( line, '\r');
        lines_.push_back( line );
    }
    while ( true );
}

// ********************************************************************************

// Starts search from line i
size_t TextFileReader_2::find( const std::string & word, const size_t i_start ) const
{
    if ( word == "" )
       return std::string::npos;
    for ( size_t i( i_start ); i != lines_.size(); ++i )
    {
        if ( lines_[i].find( word ) != std::string::npos )
            return i;
    }
    return std::string::npos;
}

// ********************************************************************************

// Starts search from line i
size_t TextFileReader_2::find_whole_word( const std::string & word, const size_t i_start ) const
{
    if ( word == "" )
       return std::string::npos;
    for ( size_t iLine( i_start ); iLine != lines_.size(); ++iLine )
    {
        // @@ I probably need a "find_whole_word" class (so that "whole word characters" can be configured).
        // If you get a match, and that turns out not to be a whole word, you still have to scan the rest
        // of the string to see if a "whole word" match occurs later in the string.
        size_t iPos = lines_[iLine].find( word );
        while ( iPos != std::string::npos )
        {
            bool start_ok( false );
            bool end_ok( false );
            if ( ( iPos == 0 ) || ( ! is_whole_word_character( lines_[iLine][iPos-1] ) ) )
                start_ok = true;
            if ( ( (iPos+word.length()) == lines_[iLine].length() ) || ( ! is_whole_word_character( lines_[iLine][iPos+word.length()] ) ) )
                end_ok = true;
            if ( start_ok && end_ok )
                return iLine;
            iPos = lines_[iLine].find( word, iPos+1 );
        }
    }
    return std::string::npos;
}

// ********************************************************************************
