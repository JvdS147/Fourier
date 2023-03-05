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

#include "PowderMatchTable.h"
#include "MillerIndices.h"
#include "StringFunctions.h"
#include "TextFileReader.h"
#include "Utilities.h"

#include <stdexcept>

//  set_name |  structure_name  | directory  | status | comment |     energy      |       density      |     cell_volume    | reduced_cell_volume |  space_group  | number_of_dof |      penalty        | energy_with_penalty |        FOM         |          B             |    MarchDollase      |            h          |          k          |          l            |    CellDeformation   |        a           |        b           |        c           |      alpha       |        beta      |      gamma
//           |                  |            |        |         | [kcal/mol/atom] |       [g/cm3]      |         [A3]       |       [A3]          |               |               |   [kcal/mol/atom]   |   [kcal/mol/atom]   |                    |                        |                      |                       |                     |                       |                      |                    |                    |                    |                  |                  |
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//    set1   | structure_000001 | jobs/job0  |  done  |  none   |       0         | 1.2529500320242304 | 1795.8846350943263 | 1795.8846350943263  | P 2_1 2_1 2_1 |       60      | 0.57516559613365126 | 0.57516559613365126 | 0.4316837101620814 | -0.037706868362922898  |  1.5513171064752962  |  -0.36060802990789309 |  0.9327174538765467 |          0            | 0.013636944428995774 | 8.2593307806227418 | 13.57628501823058  | 15.470539514891119 |        90        |        90        |         90

// ********************************************************************************

PowderMatchTable::PowderMatchTable( const FileName & file_name )
{
    TextFileReader text_file_reader( file_name );
    std::string line;
    // Read the three header lines
    if ( ! text_file_reader.get_next_line( line ) )
        throw std::runtime_error( "PowderMatchTable::PowderMatchTable( FileName ): header line missing 1." );
    if ( ! text_file_reader.get_next_line( line ) )
        throw std::runtime_error( "PowderMatchTable::PowderMatchTable( FileName ): header line missing 2." );
    if ( ! text_file_reader.get_next_line( line ) )
        throw std::runtime_error( "PowderMatchTable::PowderMatchTable( FileName ): header line missing 3." );
    Splitter splitter( "|" );
    std::vector< std::string > words;
    while ( text_file_reader.get_next_line( line ) )
    {
        words = splitter.split( line );
        if ( words.size() != 26 )
            throw std::runtime_error( "PowderMatchTable::PowderMatchTable( FileName ): words.size() != 26." );
        strip( words );
        names_.push_back( words[1] );
        densities_.push_back( string2double( words[6] ) );
        space_group_names_.push_back( words[9] );
        figures_of_merit_.push_back( string2double( words[13] ) );
        Biso_values_.push_back( string2double( words[14] ) );
        MarchDollase_values_.push_back( string2double( words[15] ) );
        PO_directions_.push_back( MillerIndices( string2integer( words[16] ), string2integer( words[17] ), string2integer( words[18] ) ) );
        cell_deformations_.push_back( string2double( words[19] ) );
        double a = string2double( words[20] );
        double b = string2double( words[21] );
        double c = string2double( words[22] );
        Angle alpha = Angle::from_degrees( string2double( words[23] ) );
        Angle beta  = Angle::from_degrees( string2double( words[24] ) );
        Angle gamma = Angle::from_degrees( string2double( words[25] ) );
        crystal_lattices_.push_back( CrystalLattice( a, b, c, alpha, beta, gamma ) );
    }
}

// ********************************************************************************

