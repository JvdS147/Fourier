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

#include "ReflectionList.h"
#include "Sort.h"
#include "TextFileReader.h"
#include "TextFileWriter.h"
#include "Utilities.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <stdexcept>

// ********************************************************************************

ReflectionList::ReflectionList()
{
}

// ********************************************************************************

void ReflectionList::push_back( const MillerIndices & miller_indices, const double F_squared, const double d_spacing, const size_t multiplicity )
{
    miller_indices_.push_back( miller_indices );
    F_squared_.push_back( F_squared );
    d_spacings_.push_back( d_spacing );
    multiplicity_.push_back( multiplicity );
    sorted_map_.push_back();
    sort_by_d_spacing();
}

// ********************************************************************************

void ReflectionList::reserve( const size_t nvalues )
{
    miller_indices_.reserve( nvalues );
    F_squared_.reserve( nvalues );
    d_spacings_.reserve( nvalues );
    multiplicity_.reserve( nvalues );
}

// ********************************************************************************

size_t ReflectionList::index( const MillerIndices & miller_indices )
{
    for ( size_t i( 0 ); i != size(); ++i )
    {
        if ( this->miller_indices( i ) == miller_indices )
            return i;
    }
    return size();
}

// ********************************************************************************

void ReflectionList::sort_by_d_spacing()
{
    // We don't actually sort the lists, but create a sorted map
    sorted_map_ = sort( d_spacings_, true );
}

// ********************************************************************************

//  -6  -3  -8    1.71    4.74
//   6   3   8   -1.06    4.75
//  -6  -3  -9    4.34    4.63
//   6   3   9   -0.95    4.77
//    void push_back( const MillerIndices & miller_indices, const double F_squared, const double d_spacing, const size_t multiplicity );

void ReflectionList::read_hkl( const FileName & file_name )
{
    *this = ReflectionList();
    TextFileReader text_file_reader( file_name );
    std::vector< std::string > words;
    while ( text_file_reader.get_next_line( words ) )
    {
        if ( words.size() != 5 )
            throw std::runtime_error( "ReflectionList::read_hkl(): cannot interpret line \"" + text_file_reader.get_line() + "\"" );
        push_back( MillerIndices( string2integer( words[0] ), string2integer( words[1] ), string2integer( words[2] ) ), string2double( words[3] ), 0.0, 0 );
    }
}

// ********************************************************************************

void ReflectionList::show() const
{
    for ( size_t i( 0 ); i != size(); ++i )
    {
        std::cout << miller_indices(i) << std::endl;
        std::cout << F_squared(i) << std::endl;
        std::cout << d_spacing(i) << std::endl;
    }
}

// ********************************************************************************

void ReflectionList::save( const FileName & file_name ) const
{
    TextFileWriter text_file_writer( file_name );
    text_file_writer.write_line( " h   k   l  d-spacing   F^2  mult" );
    for ( size_t i( 0 ); i != d_spacings_.size(); ++i )
        text_file_writer.write_line( int2string( miller_indices(i).h(), 3, ' ' ) + " " +
                                     int2string( miller_indices(i).k(), 3, ' ' ) + " " +
                                     int2string( miller_indices(i).l(), 3, ' ' ) + " " +
                                     double2string( d_spacing(i), 5, 8 ) + " " +
                                     double2string( F_squared(i), 2, 12 ) + " " +
                                     size_t2string( multiplicity(i), 2, ' ' ) );
}

// ********************************************************************************

