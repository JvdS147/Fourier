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

#include "ReadXYZ.h"
#include "Atom.h"
#include "Element.h"
#include "TextFileReader.h"
#include "Utilities.h"
#include "Vector3D.h"

//#include <fstream>
#include <stdexcept>
#include <iostream>

//#include <iostream> // For debugging

/*
160
WUBDOM
C      -1.031928   -0.865325    1.394524
C      -2.766437    0.404252    0.425293
C      -0.901004    1.591030    1.463277
C       0.838377    0.204561    2.425207
C       0.933200    2.644682    2.515255
O      -0.667888    3.912635    1.491873
C      -2.196361   -0.787397    0.723425
...
*/

// The coordinates of the atoms are Cartesian...
void read_xyz( const FileName & file_name, std::vector< Atom > & atoms )
{
    atoms.clear();
    TextFileReader text_file_reader( file_name );
    text_file_reader.set_skip_empty_lines( false );
    std::vector< std::string > words;
    std::string line;

    // Number of atoms
    if ( ! text_file_reader.get_next_line( words ) )
        throw std::runtime_error( "read_xyz(): file is empty." );
    if ( words.size() != 1 )
        throw std::runtime_error( "read_xyz(): first line must contain number of atoms." );
    int number_of_atoms = string2integer( words[ 0 ] );
    if ( number_of_atoms < 0 )
        throw std::runtime_error( "read_xyz(): file is empty." );
    // Commment
    if ( ! text_file_reader.get_next_line( line ) )
        throw std::runtime_error( "read_xyz(): file is empty." );
    size_t atom_number( 1 );
    while ( text_file_reader.get_next_line( words ) )
    {
        // Empty lines do not count, but do generate a warning
        if ( words.empty() )
        {
            std::cout << "read_xyz(): WARNING: empty line found." << std::endl;
        }
        else if ( atom_number <= number_of_atoms )
        {
            if ( words.size() != 4 )
                throw std::runtime_error( "read_xyz(): atom line does not contain four items." );
            Atom atom( Element( words[ 0 ] ), Vector3D( string2double( words[ 1 ] ),
                                                        string2double( words[ 2 ] ),
                                                        string2double( words[ 3 ] ) ),
                                                       words[ 0 ]+size_t2string( atom_number ) );
            atoms.push_back( atom );
        }
        ++atom_number;
    }
    if ( atom_number != ( number_of_atoms + 1 ) )
        std::cout << "read_xyz(): WARNING: actual number of atoms (" + size_t2string( atom_number-1 ) + ") not equal to number specified (" + size_t2string( number_of_atoms ) + ")." << std::endl;
}
