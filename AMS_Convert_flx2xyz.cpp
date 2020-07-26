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

#include "AMS_Convert_flx2xyz.h"
#include "FileName.h"
#include "TextFileReader_2.h"
#include "TextFileWriter.h"
#include "Utilities.h"

#include <stdexcept>

// 3
//      C    1        O1_0        -5.837466094        0.4745566298        -2.042695814    7           O2                   0    2    3
//      N    2        C1_0        -5.222830376        0.3042018964       -0.8392345753    3           C3                   0    1    4    5
//      O    3       H36_0        -6.437558802       -0.2925495275        -2.141493922    0            H                   0    1
//
//Flexible torsions
//    0
//
//Fixed torsions
//    0
//
//Flexible inversions
//    0
//
//Fixed inversions
//    0
//
//Rings
//    0

// ********************************************************************************

void convert_flx2xyz( const FileName & input_file_name )
{
    TextFileReader_2 input_file( input_file_name );
    TextFileWriter output_file( replace_extension( input_file_name, "xyz" ) );
    // Read number of atoms
    size_t iLine( 0 );
    std::string line = input_file.line( iLine );
    std::vector< std::string > words = split( line );
    if ( words.size() != 1 )
        throw std::runtime_error( "convert_flx2xyz(): unexpected file format." );
    // Technically an error here if number of atoms is negative
    size_t natoms = string2integer( words[0] );
    if ( natoms == 0 )
        throw std::runtime_error( "convert_flx2xyz(): unexpected file format." );
    output_file.write_line( words[0] );
    output_file.write_line();
    for ( size_t i( 0 ); i != natoms; ++i )
    {
        ++iLine;
        words = split( input_file.line( iLine ) );
        output_file.write_line( pad( words[0], 2 ) + " " + pad_plus( words[3], 15 ) + " " + pad_plus( words[4], 15 ) + " " + pad_plus( words[5], 15 ) );
    }
}

// ********************************************************************************

