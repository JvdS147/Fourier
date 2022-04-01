#ifndef GENERATEPOWDERCIF_H
#define GENERATEPOWDERCIF_H

/* *********************************************
Copyright (c) 2013-2022, Cornelis Jan (Jacco) van de Streek
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

#include "ChemicalFormula.h"
#include "CrystalLattice.h"
#include "CrystalStructure.h"
#include "PowderPattern.h"
#include "PowderPatternCalculator.h"
#include "SpaceGroup.h"
#include "TextFileReader_2.h"
#include "TextFileWriter.h"
#include "Wavelength.h"

#include <string>

// TODO: relabel
class GeneratePowderCIF
{
public:

    // @@ Calling this constructor clears all output files, that is not really desired behaviour...
    GeneratePowderCIF( const std::string & directory, const std::string & base_name );

    void generate();
    
    // The new labels of the non-hydrogen atoms are checked. They must be either C1, C2, C3, N1, N2, N3, etc. or C1, C2, N3, C4, N5, etc.
    // Hydrogen atoms are ignored.
    // @@ Could add check on type of characters that are allowed, e.g. "'" would be problematic in a cif
    // @@@ Very weird, currently can only be called *after* generate().
    void check_new_labels() const;

    enum ZoomPolicy { ALWAYS_ZOOM, NEVER_ZOOM, ZOOM_OVER_40 };

    void generate_R_input_file( const ZoomPolicy zoom_policy ); // Not const because technically the output_file_R_ object is changed. Could make it mutable and this function const.

private:
    std::string directory_;
    std::string base_name_;
    std::string file_name_pro_;
    std::string file_name_tic_;
    TextFileReader_2 file_cif_;
    TextFileReader_2 file_inp_;
    TextFileReader_2 file_Hmi_;
    TextFileReader_2 file_ext_;
    TextFileReader_2 file_bond_lengths_;
    TextFileReader_2 file_valence_angles_;
    TextFileWriter output_file_;
    TextFileWriter output_file_xrpd_;
    TextFileWriter output_file_R_;
    bool replace_hydrogen_atoms_;
    bool relabel_;
    std::map< std::string, std::string > labels_;
    size_t padding_length_;
    size_t longest_label_size_;
    size_t longest_x_;
    size_t longest_y_;
    size_t longest_z_;
    ChemicalFormula chemical_formula_;
    Wavelength wavelength_;
    CrystalLattice crystal_lattice_;
    SpaceGroup space_group_;
    CrystalStructure crystal_structure_;
    PowderPattern powder_pattern_;

// @@@ Should any of these now be const (quite a lot has changed since they were added)

void insert_keyword( const std::string & keyword );
void insert_keyword_and_value( const std::string & keyword, const std::string & value );
void write_bond_part();
void write_angle_part();
std::string relabel( const std::string & old_label ) const;

};

#endif // GENERATEPOWDERCIF_H

