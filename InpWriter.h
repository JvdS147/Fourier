#ifndef INPWRITER_H
#define INPWRITER_H

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

class CrystalLattice;
class FileName;
class TextFileWriter;

#include <string>

// This is the bit that is constant and global.
void write_preamble( TextFileWriter & text_file_writer );

// From .cif file
// aal = additional atom label, a string that is appended to each atom label
void inp_writer( const FileName & input_cif_file_name, const FileName & input_xye_file_name, const std::string & aal = "" );

// Each atoms is tethered to its initial postion with a little spring.
void inp_writer_distance_restraints( const FileName & input_cif_file_name, const FileName & input_xye_file_name, const std::string & aal = "" );

// Only call this function if the refinement is parametric.
void write_unit_cell_variables( TextFileWriter & text_file_writer, const CrystalLattice & crystal_lattice, const bool variable_temperature = false );

void write_unit_cell( TextFileWriter & text_file_writer, const CrystalLattice & crystal_lattice, const bool parametric_refinement = false, const bool variable_temperature = false );

#endif // INPWRITER_H

