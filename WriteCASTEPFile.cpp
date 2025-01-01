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

#include "WriteCASTEPFile.h"
#include "FileName.h"
#include "StringFunctions.h"
#include "TextFileWriter.h"
#include "Utilities.h"

#include <stdexcept>

// ********************************************************************************

WriteCASTEPFile::WriteCASTEPFile() :
job_type_(UNIT_CELL_FIXED),
cut_off_energy_(520.0)
{
}

// ********************************************************************************

WriteCASTEPFile::WriteCASTEPFile( const CrystalStructure & crystal_structure,
                                  const std::string & directory,
                                  const std::string & base_name ) :
job_type_(UNIT_CELL_FIXED),
base_name_(base_name),
crystal_structure_(crystal_structure),
cut_off_energy_(520.0)
{
    directory_ = append_backslash( directory );
    crystal_structure_.apply_space_group_symmetry();
}

// ********************************************************************************

void WriteCASTEPFile::set_directory( const std::string & directory )
{
    directory_ = append_backslash( directory );
}

// ********************************************************************************

void WriteCASTEPFile::write() const
{
    if ( base_name_.empty() )
        throw std::runtime_error( "WriteCASTEPFile::write(): base_name not defined." );
    TextFileWriter text_file_writer( FileName( directory_, base_name_, "cell" ) );

// MEXZOG
//a 5.10141(12) b 5.53079(13) c 9.0323(2)
//a 79.0266(10) b 75.7579(9) g 78.8261(9)
//%BLOCK LATTICE_CART
//       4.885292166468590       0.763649905894596       1.255047113684490
//       0.000000000000000       5.429663167226080       1.052803834799450
//       0.000000000000000       0.000000000000000       9.032299999999999
//%ENDBLOCK LATTICE_CART

    Matrix3D lattice = crystal_structure_.crystal_lattice().for_CASTEP();
    text_file_writer.write_line( "%BLOCK LATTICE_CART" );
    text_file_writer.write_line( "       " + double2string( lattice.value( 0, 0 ), 5, 9 ) +
                                       " " + double2string( lattice.value( 0, 1 ), 5, 9 ) +
                                       " " + double2string( lattice.value( 0, 2 ), 5, 9 ) );
    text_file_writer.write_line( "       " + double2string( lattice.value( 1, 0 ), 5, 9 ) +
                                       " " + double2string( lattice.value( 1, 1 ), 5, 9 ) +
                                       " " + double2string( lattice.value( 1, 2 ), 5, 9 ) );
    text_file_writer.write_line( "       " + double2string( lattice.value( 2, 0 ), 5, 9 ) +
                                       " " + double2string( lattice.value( 2, 1 ), 5, 9 ) +
                                       " " + double2string( lattice.value( 2, 2 ), 5, 9 ) );
    text_file_writer.write_line( "%ENDBLOCK LATTICE_CART" );
    std::set< Element > elements = crystal_structure_.elements();
    text_file_writer.write_line();
    text_file_writer.write_line( "%BLOCK POSITIONS_FRAC" );
    for ( std::set< Element >::const_iterator it( elements.begin() ); it != elements.end(); ++it )
    {
        for ( size_t i( 0 ); i != crystal_structure_.natoms(); ++i )
        {
            if ( *it == crystal_structure_.atom( i ).element() )
                text_file_writer.write_line( "  " + crystal_structure_.atom( i ).element().symbol() + " " +
                                             double2string_pad_plus( crystal_structure_.atom( i ).position().x(), 5 ) + " " +
                                             double2string_pad_plus( crystal_structure_.atom( i ).position().y(), 5 ) + " " +
                                             double2string_pad_plus( crystal_structure_.atom( i ).position().z(), 5 ) + " " );
        }
    }
    text_file_writer.write_line( "%ENDBLOCK POSITIONS_FRAC" );
    text_file_writer.write_line();
    text_file_writer.write_line( "KPOINTS_MP_SPACING 0.07" );
    text_file_writer.write_line();
    text_file_writer.write_line( "%BLOCK SYMMETRY_OPS" );
    for ( size_t i( 0 ); i != crystal_structure_.space_group().nsymmetry_operators(); ++i )
    {
        Matrix3D rotation = crystal_structure_.space_group().symmetry_operator( i ).rotation();
        Vector3D translation = crystal_structure_.space_group().symmetry_operator( i ).translation();
        text_file_writer.write_line( " " + double2string_pad_plus( rotation.value( 0, 0 ), 5 ) +
                                     " " + double2string_pad_plus( rotation.value( 0, 1 ), 5 ) +
                                     " " + double2string_pad_plus( rotation.value( 0, 2 ), 5 ) );
        text_file_writer.write_line( " " + double2string_pad_plus( rotation.value( 1, 0 ), 5 ) +
                                     " " + double2string_pad_plus( rotation.value( 1, 1 ), 5 ) +
                                     " " + double2string_pad_plus( rotation.value( 1, 2 ), 5 ) );
        text_file_writer.write_line( " " + double2string_pad_plus( rotation.value( 2, 0 ), 5 ) +
                                     " " + double2string_pad_plus( rotation.value( 2, 1 ), 5 ) +
                                     " " + double2string_pad_plus( rotation.value( 2, 2 ), 5 ) );
        text_file_writer.write_line( " " + double2string_pad_plus( translation.value( 0 ), 5 ) +
                                     " " + double2string_pad_plus( translation.value( 1 ), 5 ) +
                                     " " + double2string_pad_plus( translation.value( 2 ), 5 ) );
    }
    text_file_writer.write_line( "%ENDBLOCK SYMMETRY_OPS" );
    text_file_writer.write_line();
    if ( job_type_ == UNIT_CELL_FREE )
        text_file_writer.write_line( "FIX_ALL_CELL : false" );
    else
        text_file_writer.write_line( "FIX_ALL_CELL : true" );
    text_file_writer.write_line();
    text_file_writer.write_line( "FIX_COM : false" );
    text_file_writer.write_line();
    if ( job_type_ == H_ATOMS_ONLY )
    {
        text_file_writer.write_line( "%BLOCK IONIC_CONSTRAINTS" );
        size_t counter_1( 0 );
        for ( std::set< Element >::const_iterator it( elements.begin() ); it != elements.end(); ++it )
        {
            if ( it->is_H_or_D() )
                continue;
            size_t counter_2( 0 );
            for ( size_t i( 0 ); i != crystal_structure_.natoms(); ++i )
            {
                if ( *it == crystal_structure_.atom( i ).element() )
                {
                    ++counter_2;
                    ++counter_1;
                    text_file_writer.write_line( "  " + size_t2string( counter_1, 4, ' ' ) + " " +
                                                 crystal_structure_.atom( i ).element().symbol() + " " +
                                                 size_t2string( counter_2, 4, ' ' ) + " " +
                                                 "1.00000   0.00000   0.00000" );
                    ++counter_1;
                    text_file_writer.write_line( "  " + size_t2string( counter_1, 4, ' ' ) + " " +
                                                 crystal_structure_.atom( i ).element().symbol() + " " +
                                                 size_t2string( counter_2, 4, ' ' ) + " " +
                                                 "0.00000   1.00000   0.00000" );
                    ++counter_1;
                    text_file_writer.write_line( "  " + size_t2string( counter_1, 4, ' ' ) + " " +
                                                 crystal_structure_.atom( i ).element().symbol() + " " +
                                                 size_t2string( counter_2, 4, ' ' ) + " " +
                                                 "0.00000   0.00000   1.00000" );
                }
            }
        }
        text_file_writer.write_line( "%ENDBLOCK IONIC_CONSTRAINTS" );
    }
//    text_file_writer.write_line();
//    text_file_writer.write_line( "%BLOCK EXTERNAL_PRESSURE" );
//    text_file_writer.write_line( "    0.0000000000    0.0000000000    0.0000000000" );
//    text_file_writer.write_line( "                    0.0000000000    0.0000000000" );
//    text_file_writer.write_line( "                                    0.0000000000" );
//    text_file_writer.write_line( "%ENDBLOCK EXTERNAL_PRESSURE" );
//    text_file_writer.write_line();
//    text_file_writer.write_line( "%BLOCK SPECIES_MASS" );
//    for ( std::set< Element >::const_iterator it( elements.begin() ); it != elements.end(); ++it )
//        text_file_writer.write_line( " " + it->symbol() + " " + double2string( it->atomic_weight() ) );
//    text_file_writer.write_line( "%ENDBLOCK SPECIES_MASS" );
    text_file_writer.write_line();
    text_file_writer.write_line( "%BLOCK SPECIES_POT" );
    for ( std::set< Element >::const_iterator it( elements.begin() ); it != elements.end(); ++it )
        text_file_writer.write_line( " " + it->symbol() + " " + it->symbol() + "_00PBE.usp" );
    text_file_writer.write_line( "%ENDBLOCK SPECIES_POT" );
    text_file_writer.write_line();
    text_file_writer.write_line( "%BLOCK SPECIES_LCAO_STATES" );
    for ( std::set< Element >::const_iterator it( elements.begin() ); it != elements.end(); ++it )
    {
        switch ( it->atomic_number() )
        {
            case  1 : text_file_writer.write_line( " " + it->symbol() + " 1" ); break;
            case 20 : text_file_writer.write_line( " " + it->symbol() + " 4" ); break; // Ca
            default : text_file_writer.write_line( " " + it->symbol() + " 2" ); // C, N, O, S, Cl
        }
    }
    text_file_writer.write_line( "%ENDBLOCK SPECIES_LCAO_STATES" );

    TextFileWriter text_file_writer_2( FileName( directory_, base_name_, "param" ) );
    if ( job_type_ != SS_NMR )
        text_file_writer_2.write_line( "task : GeometryOptimization" );
    text_file_writer_2.write_line( "xc_functional : PBE" );
    text_file_writer_2.write_line( "sedc_apply : true" );
    text_file_writer_2.write_line( "sedc_scheme : G06" );
    text_file_writer_2.write_line( "spin_polarized : false" );
    text_file_writer_2.write_line( "opt_strategy : Default" );
    text_file_writer_2.write_line( "page_wvfns :        0" );
    text_file_writer_2.write_line( "cut_off_energy :      " + double2string( cut_off_energy_, 5 ) );
    text_file_writer_2.write_line( "grid_scale :        1.500000000000000" );
    text_file_writer_2.write_line( "fine_grid_scale :        1.500000000000000" );
    if ( job_type_ == UNIT_CELL_FREE )
    {
        text_file_writer_2.write_line( "finite_basis_corr :        2" );
        text_file_writer_2.write_line( "finite_basis_npoints :        3" );
    }
    else
        text_file_writer_2.write_line( "finite_basis_corr :        0" );
    text_file_writer_2.write_line( "elec_energy_tol :   1.000000000000000e-006" );
    text_file_writer_2.write_line( "max_scf_cycles :      200" );
    text_file_writer_2.write_line( "fix_occupancy : true" );
    text_file_writer_2.write_line( "metals_method : dm" );
    text_file_writer_2.write_line( "mixing_scheme : Pulay" );
    text_file_writer_2.write_line( "mix_charge_amp :        0.500000000000000" );
    text_file_writer_2.write_line( "mix_charge_gmax :        1.500000000000000" );
    text_file_writer_2.write_line( "mix_history_length :       20" );
    text_file_writer_2.write_line( "nextra_bands : 0" );
    if ( job_type_ != SS_NMR )
    {
        text_file_writer_2.write_line( "geom_energy_tol :   1.000000000000000e-005" );
        text_file_writer_2.write_line( "geom_force_tol :        0.030000000000000" );
        text_file_writer_2.write_line( "geom_stress_tol :        0.050000000000000" );
        text_file_writer_2.write_line( "geom_disp_tol :   1.000000000000000e-003" );
        text_file_writer_2.write_line( "geom_max_iter :     1000" );
        text_file_writer_2.write_line( "geom_method : BFGS" );
        text_file_writer_2.write_line( "fixed_npw : false" );
    }
    if ( job_type_ == UNIT_CELL_FREE )
        text_file_writer_2.write_line( "geom_modulus_est :       250.00000000000000  GPa" );
    text_file_writer_2.write_line( "calculate_ELF : false" );
    if ( job_type_ == UNIT_CELL_FREE )
        text_file_writer_2.write_line( "calculate_stress : true" );
    else
        text_file_writer_2.write_line( "calculate_stress : false" );
    text_file_writer_2.write_line( "popn_calculate : false" );
    text_file_writer_2.write_line( "calculate_hirshfeld : false" );
    text_file_writer_2.write_line( "calculate_densdiff : false" );
    text_file_writer_2.write_line( "pdos_calculate_weights : false" );
    text_file_writer_2.write_line( "num_dump_cycles : 0" );
    if ( job_type_ == SS_NMR )
    {
        text_file_writer_2.write_line( "task : MagRes" );
        text_file_writer_2.write_line( "magres_method : crystal" );
        text_file_writer_2.write_line( "magres_task : NMR" );
        text_file_writer_2.write_line( "magres_max_cg_steps :      250" );
        text_file_writer_2.write_line( "bs_max_iter : 250" );
        text_file_writer_2.write_line( "bs_max_cg_steps : 5" );
        text_file_writer_2.write_line( "bs_eigenvalue_tol : 1.0e-9" );
        text_file_writer_2.write_line( "bs_write_eigenvalues : false" );
    }
}

// ********************************************************************************

