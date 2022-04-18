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

#include "RESULTSFile.h"
#include "MathsFunctions.h"
#include "TextFileReader.h"
#include "Utilities.h"

#include <stdexcept>

//set1 structure_001_00 -28213.7328 done ---- -1 jobs/job0 ---- none
//set1 structure_001_01 -28213.70532 done ---- -1 jobs/job1 ---- none
//set1 structure_002_00 -28213.73293 done ---- -1 jobs/job2 ---- none
//set1 structure_002_01 -28213.70614 done ---- -1 jobs/job3 ---- none

// ********************************************************************************

RESULTSFile::RESULTSFile( const FileName & file_name, const size_t natoms )
{
    natoms_ = natoms;
    TextFileReader text_file_reader( file_name );
    std::vector< std::string > words;
    while ( text_file_reader.get_next_line( words ) )
    {
        if ( words.size() != 9 )
            throw std::runtime_error( "RESULTSFile::RESULTSFile( FileName ): words.size() != 9." );
        if ( words[3] != "done" )
            continue;
        names_.push_back( words[1] );
        raw_energies_.push_back( string2double( words[2] ) );
    }
    lowest_energy_ = empty() ? 0.0 : calculate_minimum( raw_energies_ );
}

// ********************************************************************************

std::string RESULTSFile::name( const size_t i ) const
{
    if ( i < size() )
        return names_[i];
    throw std::runtime_error( "RESULTSFile::name( size_t ): out of bounds." );
}

// ********************************************************************************

double RESULTSFile::energy( const size_t i ) const
{
    if ( i < size() )
        return static_cast< double >( natoms_ ) * ( raw_energies_[i] - lowest_energy_ );
    throw std::runtime_error( "RESULTSFile::energy( size_t ): out of bounds." );
}

// ********************************************************************************

std::vector< double > RESULTSFile::energies() const
{
    std::vector< double > result;
    result.reserve( raw_energies_.size() );
    for ( size_t i( 0 ); i != raw_energies_.size(); ++i )
        result.push_back( energy( i ) );
    return result;
}

// ********************************************************************************

double RESULTSFile::raw_energy( const size_t i ) const
{
    if ( i < size() )
        return raw_energies_[i];
    throw std::runtime_error( "RESULTSFile::raw_energy( size_t ): v." );
}

// ********************************************************************************

void RESULTSFile::reduce_to_base_name()
{
    for ( size_t i( 0 ); i != names_.size(); ++i )
    {
        if ( names_[i].length() < 4 )
            continue;
        if ( names_[i][ names_[i].length()-3 ] != '_' )
            continue;
        if ( ! is_digit( names_[i][ names_[i].length()-2 ] ) )
            continue;
        if ( ! is_digit( names_[i][ names_[i].length()-1 ] ) )
            continue;
        names_[i] = names_[i].substr( 0, names_[i].length() - 3 );
    }
}

// ********************************************************************************

