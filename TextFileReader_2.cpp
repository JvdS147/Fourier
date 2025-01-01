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

#include "TextFileReader_2.h"
#include "FileName.h"
#include "StringFunctions.h"

#include <stdexcept>
#include <iostream>

// ********************************************************************************

TextFileReader_2::TextFileReader_2( const FileName & file_name )
{
    read_file( file_name );
}

// ********************************************************************************

TextFileReader_2::TextFileReader_2( const std::vector< std::string > & lines )
{
    lines_ = lines;
}

// ********************************************************************************

void TextFileReader_2::read_file( const FileName & file_name )
{
    std::ifstream input_file( file_name.full_name().c_str() );
    if ( ! input_file )
       throw std::runtime_error( std::string( "TextFileReader_2::read_file(): Could not open file " ) + file_name.full_name() );
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

std::string TextFileReader_2::line( const size_t i ) const
{
    if ( i < lines_.size() )
        return lines_[i];
    throw std::runtime_error( "TextFileReader_2::line( size_t ): out of bounds." );
}

// ********************************************************************************

// Deletes entire lines.
// comment_identifier must start at the start of the line.
void TextFileReader_2::purge_comment_lines( std::string comment_identifier, const bool case_sensitive )
{
    if ( ! case_sensitive )
        comment_identifier = to_upper( comment_identifier );
    size_t iPos( 0 );
    for ( size_t i( 0 ); i != lines_.size(); ++i )
    {
        bool is_comment( false );
        if ( case_sensitive )
            is_comment = ( lines_[i].substr( 0, comment_identifier.length() ) == comment_identifier );
        else
            is_comment = ( to_upper( lines_[i].substr( 0, comment_identifier.length() ) ) == comment_identifier );
        if ( ! is_comment )
        {
            if ( iPos != i )
                lines_[iPos] = lines_[i];
            ++iPos;
        }
    }
    lines_.resize( iPos );
}

// ********************************************************************************

// Empty means empty, a line with only whitespace is not empty
void TextFileReader_2::purge_empty_lines()
{
    size_t iPos( 0 );
    for ( size_t i( 0 ); i != lines_.size(); ++i )
    {
        if ( ! lines_[i].empty() )
        {
            if ( iPos != i )
                lines_[iPos] = lines_[i];
            ++iPos;
        }
    }
    lines_.resize( iPos );    
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

