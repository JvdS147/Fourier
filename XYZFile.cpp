/* *********************************************
Copyright (c) 2013-2026, Cornelis Jan (Jacco) van de Streek
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

#include "XYZFile.h"
#include "CrystalStructure.h"
#include "FileName.h"
#include "StringFunctions.h"
#include "TextFileWriter.h"
#include "Utilities.h"

#include <vector>

// ********************************************************************************

void save_as_xyz( const CrystalStructure & crystal_structure, const FileName & file_name )
{
    TextFileWriter text_file_writer( file_name );
    text_file_writer.write_line( size_t2string( crystal_structure.natoms() ) );
    text_file_writer.write_line( "Cell=\"" + double2string( crystal_structure.crystal_lattice().a() ) + " " +
                                             double2string( crystal_structure.crystal_lattice().b() ) + " " +
                                             double2string( crystal_structure.crystal_lattice().c() ) + " " +
                                             double2string( crystal_structure.crystal_lattice().alpha().value_in_degrees() ) + " " +
                                             double2string( crystal_structure.crystal_lattice().beta().value_in_degrees() ) + " " +
                                             double2string( crystal_structure.crystal_lattice().gamma().value_in_degrees() ) + "\" Name=\"" + crystal_structure.name() + "\"" );
    for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
    {
        Vector3D position = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom(i).position() );
        text_file_writer.write_line( pad( crystal_structure.atom(i).element().symbol(), 2 ) + " " + double2string( position.x() ) + " " +
                                                                                                    double2string( position.y() ) + " " + 
                                                                                                    double2string( position.z() ) );
    }
}

// ********************************************************************************

