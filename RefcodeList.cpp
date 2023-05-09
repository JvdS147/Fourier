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

#include "RefcodeList.h"
//#include "Sort.h"
#include "TextFileReader.h"
#include "TextFileWriter.h"
//#include "Utilities.h"

//#include <algorithm>
#include <iostream>
#include <fstream>
#include <stdexcept>

// ********************************************************************************

RefcodeList::RefcodeList()
{
}

// ********************************************************************************

RefcodeList::RefcodeList( const std::vector< std::string > & values )
{
    for ( size_t i( 0 ); i != values.size(); ++i )
        push_back( Refcode( values[ i ] ) );
}

// ********************************************************************************

void RefcodeList::push_back( const Refcode & refcode )
{
    refcodes_.push_back( refcode );
    sorted_map_.push_back();
}

// ********************************************************************************

void RefcodeList::reserve( const size_t nvalues )
{
    refcodes_.reserve( nvalues );
}

// ********************************************************************************

bool RefcodeList::contains( const Refcode & refcode ) const
{
    for ( size_t i( 0 ); i != size(); ++i )
    {
        if ( refcodes_[ i ] == refcode )
            return true;
    }
    return false;
}

// ********************************************************************************

size_t RefcodeList::index( const Refcode & refcode ) const
{
    for ( size_t i( 0 ); i != size(); ++i )
    {
        if ( this->refcode( i ) == refcode )
            return i;
    }
    return size();
}

// ********************************************************************************

// Converts each entry to its refcode family, then removes duplicates.
void RefcodeList::convert_to_unique_families()
{
    RefcodeList new_list;
    for ( size_t i( 0 ); i != refcodes_.size(); ++i )
    {
        if ( ! new_list.contains( refcode( i ).family() ) )
            new_list.push_back( refcode( i ).family() );
    }
    *this = new_list;
}

// ********************************************************************************

void RefcodeList::read_gcd( const FileName & file_name )
{
    *this = RefcodeList();
    TextFileReader text_file_reader( file_name );
    std::vector< std::string > words;
    while ( text_file_reader.get_next_line( words ) )
    {
        if ( words.size() != 1 )
            throw std::runtime_error( "RefcodeList::read_gcd(): cannot interpret line \"" + text_file_reader.get_line() + "\"" );
        push_back( Refcode( words[0] ) );
    }
}

// ********************************************************************************

void RefcodeList::show() const
{
    for ( size_t i( 0 ); i != size(); ++i )
        std::cout << refcode( i ).value() << std::endl;
}

// ********************************************************************************

void RefcodeList::save( const FileName & file_name ) const
{
    TextFileWriter text_file_writer( file_name );
    for ( size_t i( 0 ); i != refcodes_.size(); ++i )
        text_file_writer.write_line( refcode( i ).value() );
}

// ********************************************************************************

