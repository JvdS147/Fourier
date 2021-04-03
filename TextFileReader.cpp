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

#include "TextFileReader.h"
#include "FileName.h"
#include "Utilities.h"

#include <stdexcept>
#include <iostream>

// ********************************************************************************

TextFileReader::TextFileReader( const FileName & file_name ):
input_file_( file_name.full_name().c_str() ),
line_number_(0),
skip_empty_lines_(false),
allow_single_quotes_(false),
push_back_last_line_(false)
{
    if ( ! input_file_ )
       throw std::runtime_error( std::string( "TextFileReader::TextFileReader(): Could not open file " ) + file_name.full_name() );
}

// ********************************************************************************

bool TextFileReader::get_next_line( std::vector< std::string > & words )
{
    std::string line;
    bool return_code = get_next_line( line );
    if ( return_code )
    {
        if ( allow_single_quotes_ )
            words = split_2( line_ );
        else
            words = split( line_ );
    }
    return return_code;
}

// ********************************************************************************

bool TextFileReader::get_next_line( std::string & line )
{
    if ( push_back_last_line_ )
    {
        line = line_;
        push_back_last_line_ = false;
        return true;
    }
    do
    {
        if ( ! getline( input_file_, line_ ) )
            return false;
        // remove \r
        line_ = remove( line_, '\r');
        ++line_number_;
        // Skip comments
        bool is_comment( false );
        for ( size_t i( 0 ); i != comment_identifiers_.size(); ++i )
        {
            if ( comment_identifiers_[i].size() <= line_.size() )
            {
                 if ( line_.substr( 0, comment_identifiers_[i].size() ) == comment_identifiers_[i] )
                 {
                    is_comment = true;
                    break;
                 }
            }
        }
        if ( is_comment )
            continue;
        if ( skip_empty_lines_ && line_.empty() )
            continue;
        if ( skip_whitespace_only_lines_ && strip(line_).empty() )
            continue;            
        line = line_;
        return true;
    }
    while ( true );
}

// ********************************************************************************

