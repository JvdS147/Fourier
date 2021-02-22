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

#include "AnalyseTrajectory.h"
#include "3DCalculations.h"
#include "AnisotropicDisplacementParameters.h"
#include "CrystalStructure.h"
#include "DoubleWithESD.h"
#include "Element.h"
#include "ReadCif.h"
#include "TextFileWriter.h"
#include "Utilities.h"

// ********************************************************************************

AnalyseTrajectory::AnalyseTrajectory( const FileList file_list,
                                      const size_t u,
                                      const size_t v,
                                      const size_t w,
                                      const SpaceGroup & space_group,
                                      const Matrix3D & transformation ) :
file_list_(file_list),
u_(u),
v_(v),
w_(w),
space_group_(space_group),
transformation_(transformation),
drift_correction_(USE_FIRST_FRAME),
write_lean_(false),
write_average_(true),
write_average_noH_(false),
write_average_ESDs_(false),
write_sum_(false)
{
    analyse();
}

// ********************************************************************************

void AnalyseTrajectory::analyse()
{
    file_list_.set_prepend_file_name_with_basedirectory( true );
    std::string directory;
    std::string file_name;
    std::string extension;
    std::vector< Element > elements;
    std::vector< RunningAverageAndESD< Vector3D > > average_positions; // Fractional coordinates
    size_t natoms;
    std::vector< std::vector< Vector3D > > fractional_positions_trajectory;
    // Read the first cif file and initialise everything
    {
    CrystalStructure crystal_structure;
    std::cout << "Now reading cif... " + file_list_.value( 0 ).full_name() << std::endl;
    read_cif( file_list_.value( 0 ), crystal_structure );
    if ( write_lean_ )
        crystal_structure.save_cif( append_to_file_name( file_list_.value( 0 ), "_lean" ) );
    CrystalLattice crystal_lattice = crystal_structure.crystal_lattice();
    crystal_structure.set_space_group( space_group_ );
    average_a_.add_value( crystal_lattice.a() / u_ );
    average_b_.add_value( crystal_lattice.b() / v_ );
    average_c_.add_value( crystal_lattice.c() / w_ );
    average_alpha_.add_value( crystal_lattice.alpha() );
    average_beta_.add_value( crystal_lattice.beta() );
    average_gamma_.add_value( crystal_lattice.gamma() );
    average_volume_.add_value( crystal_lattice.volume() / ( u_ * v_ * w_ ) );
    std::vector< std::vector< Vector3D > > fractional_positions_frame;
    Vector3D actual_centre;
    // Returns a std::vector of atomic coordinates for each atom in the asymmetric unit
    if ( ( drift_correction_ == NONE ) ||
         ( drift_correction_ == USE_FIRST_FRAME ) )
    {
        crystal_structure.collapse_supercell( u_, v_, w_, 0, drift_correction_vector_, transformation_, actual_centre, fractional_positions_frame );
        drift_correction_vector_ = actual_centre;
    }
    else
        crystal_structure.collapse_supercell( u_, v_, w_, drift_correction_, drift_correction_vector_, transformation_, actual_centre, fractional_positions_frame );
    centres_of_mass_.push_back( actual_centre );
    crystal_structure.transform( transformation_ );
    natoms = fractional_positions_frame.size();
    elements.reserve( natoms );
    average_positions.reserve( natoms );
    for ( size_t i( 0 ); i != natoms; ++i )
    {
        elements.push_back( crystal_structure.atom( i ).element() );
        RunningAverageAndESD< Vector3D > average_position;
        std::vector< Vector3D > temp_fractional_positions;
        for ( size_t j( 0 ); j != fractional_positions_frame[i].size(); ++j )
        {
            average_position.add_value( fractional_positions_frame[i][j] );
            temp_fractional_positions.push_back( fractional_positions_frame[i][j] );
        }
        average_positions.push_back( average_position );
        fractional_positions_trajectory.push_back( temp_fractional_positions );
    }
    }
    // Read the remaining cif files
    for ( size_t i( 1 ); i != file_list_.size(); ++i )
    {
        std::cout << "Now reading cif... " + file_list_.value( i ).full_name() << std::endl;
        CrystalStructure crystal_structure;
        read_cif( file_list_.value( i ), crystal_structure );
        if ( write_lean_ )
            crystal_structure.save_cif( append_to_file_name( file_list_.value( i ), "_lean" ) );
        CrystalLattice crystal_lattice = crystal_structure.crystal_lattice();
        crystal_structure.set_space_group( space_group_ );
        average_a_.add_value( crystal_lattice.a() / u_ );
        average_b_.add_value( crystal_lattice.b() / v_ );
        average_c_.add_value( crystal_lattice.c() / w_ );
        average_alpha_.add_value( crystal_lattice.alpha() );
        average_beta_.add_value( crystal_lattice.beta() );
        average_gamma_.add_value( crystal_lattice.gamma() );
        average_volume_.add_value( crystal_lattice.volume() / ( u_ * v_ * w_ ) );
        std::vector< std::vector< Vector3D > > fractional_positions_frame;
        // Returns a std::vector of atomic coordinates for each atom in the asymmetric unit
        Vector3D actual_centre;
        crystal_structure.collapse_supercell( u_, v_, w_, drift_correction_, drift_correction_vector_, transformation_, actual_centre, fractional_positions_frame );
        centres_of_mass_.push_back( actual_centre );
        crystal_structure.transform( transformation_ );
        if ( fractional_positions_frame.size() != natoms )
            throw std::runtime_error( "AnalyseTrajectory::analyse(): The number of atoms in the cif files is not the same, the average cif could not be generated." );
        for ( size_t i( 0 ); i != natoms; ++i )
        {
            for ( size_t j( 0 ); j != fractional_positions_frame[i].size(); ++j )
            {
                average_positions[i].add_value( fractional_positions_frame[i][j] );
                fractional_positions_trajectory[i].push_back( fractional_positions_frame[i][j] );
            }
        }
    }
    CrystalLattice crystal_lattice_average( average_a_.average(),
                                            average_b_.average(),
                                            average_c_.average(),
                                            average_alpha_.average(),
                                            average_beta_.average(),
                                            average_gamma_.average() );
    if ( write_average_ )
    {
        TextFileWriter text_file_writer( FileName( file_list_.base_directory(), "average_adps", "cif" ) );
        text_file_writer.write_line( "data_average" );
        text_file_writer.write_line( "_symmetry_space_group_name_H-M  '" + space_group_.name() + "'" );
    //    text_file_writer.write_line( "_symmetry_Int_Tables_number     1" );
    //    text_file_writer.write_line( "_symmetry_cell_setting          triclinic" );
        text_file_writer.write_line( "_cell_length_a    " + double2string( average_a_.average(), 6 ) );
        text_file_writer.write_line( "_cell_length_b    " + double2string( average_b_.average(), 6 ) );
        text_file_writer.write_line( "_cell_length_c    " + double2string( average_c_.average(), 6 ) );
        text_file_writer.write_line( "_cell_angle_alpha " + double2string( average_alpha_.average().value_in_degrees(), 6 ) );
        text_file_writer.write_line( "_cell_angle_beta  " + double2string( average_beta_.average().value_in_degrees() , 6 ) );
        text_file_writer.write_line( "_cell_angle_gamma " + double2string( average_gamma_.average().value_in_degrees(), 6 ) );
        text_file_writer.write_line( "_cell_volume      " + double2string( average_volume_.average(), 6 ) );
        text_file_writer.write_line( "loop_" );
        text_file_writer.write_line( "_symmetry_equiv_pos_site_id" );
        text_file_writer.write_line( "_symmetry_equiv_pos_as_xyz" );
        for ( size_t i( 0 ); i != space_group_.nsymmetry_operators(); ++i )
            text_file_writer.write_line( size_t2string( i+1 ) + " " + space_group_.symmetry_operator( i ).to_string() );
        text_file_writer.write_line( "loop_" );
        text_file_writer.write_line( "_atom_site_label" ); // Needed for Materials Studio
        text_file_writer.write_line( "_atom_site_type_symbol" );
        text_file_writer.write_line( "_atom_site_fract_x" );
        text_file_writer.write_line( "_atom_site_fract_y" );
        text_file_writer.write_line( "_atom_site_fract_z" );
        size_t len( 2 );
        size_t current_size = 99;
        while ( natoms >= current_size )
        {
            ++len;
            current_size = 10 * current_size + 9;
        }
        for ( size_t i( 0 ); i != natoms; ++i )
        {
            // This is just too weird, need std::vector< DoubleWithESD > for this.
            text_file_writer.write_line( elements[i].symbol() + size_t2string( i + 1, len, '0' ) + " " +
                                         elements[i].symbol() + " " +
                                         double2string( adjust_for_translations( average_positions[ i ].average().x() ), 6 ) + " " +
                                         double2string( adjust_for_translations( average_positions[ i ].average().y() ), 6 ) + " " +
                                         double2string( adjust_for_translations( average_positions[ i ].average().z() ), 6 ) );
        }
        text_file_writer.write_line( "loop_" );
        text_file_writer.write_line( "_atom_site_aniso_label" );
        text_file_writer.write_line( "_atom_site_aniso_U_11" );
        text_file_writer.write_line( "_atom_site_aniso_U_22" );
        text_file_writer.write_line( "_atom_site_aniso_U_33" );
        text_file_writer.write_line( "_atom_site_aniso_U_12" );
        text_file_writer.write_line( "_atom_site_aniso_U_13" );
        text_file_writer.write_line( "_atom_site_aniso_U_23" );
        for ( size_t i( 0 ); i != natoms; ++i )
        {
            std::vector< Vector3D > cartesian_positions;
            for ( size_t j( 0 ); j != fractional_positions_trajectory[i].size(); ++j )
                cartesian_positions.push_back( crystal_lattice_average.fractional_to_orthogonal_matrix() * fractional_positions_trajectory[i][j] );
            AnisotropicDisplacementParameters adps( cartesian_positions );
            SymmetricMatrix3D Ucif = adps.U_cif( crystal_lattice_average );
            // This can go wrong if the ADP is very small, it may be printed like "1E-14" which will not be recognised in the cif.
            text_file_writer.write_line( elements[i].symbol() + size_t2string( i + 1, len, '0' ) + " " +
                                         double2string( Ucif.value( 0, 0 ), 6 ) + " " +
                                         double2string( Ucif.value( 1, 1 ), 6 ) + " " +
                                         double2string( Ucif.value( 2, 2 ), 6 ) + " " +
                                         double2string( Ucif.value( 0, 1 ), 6 ) + " " +
                                         double2string( Ucif.value( 0, 2 ), 6 ) + " " +
                                         double2string( Ucif.value( 1, 2 ), 6 ) );
        }
        text_file_writer.write_line();
        text_file_writer.write_line( "#END" );
    }
    if ( write_average_noH_ )
    {
        TextFileWriter text_file_writer( FileName( file_list_.base_directory(), "average_noH", "cif" ) );
        text_file_writer.write_line( "data_average" );
        text_file_writer.write_line( "_symmetry_space_group_name_H-M  '" + space_group_.name() + "'" );
    //    text_file_writer.write_line( "_symmetry_Int_Tables_number     1" );
    //    text_file_writer.write_line( "_symmetry_cell_setting          triclinic" );
        text_file_writer.write_line( "_cell_length_a    " + double2string( average_a_.average(), 6 ) );
        text_file_writer.write_line( "_cell_length_b    " + double2string( average_b_.average(), 6 ) );
        text_file_writer.write_line( "_cell_length_c    " + double2string( average_c_.average(), 6 ) );
        text_file_writer.write_line( "_cell_angle_alpha " + double2string( average_alpha_.average().value_in_degrees(), 6 ) );
        text_file_writer.write_line( "_cell_angle_beta  " + double2string( average_beta_.average().value_in_degrees(), 6 ) );
        text_file_writer.write_line( "_cell_angle_gamma " + double2string( average_gamma_.average().value_in_degrees(), 6 ) );
        text_file_writer.write_line( "_cell_volume      " + double2string( average_volume_.average(), 6 ) );
        text_file_writer.write_line( "loop_" );
        text_file_writer.write_line( "_symmetry_equiv_pos_site_id" );
        text_file_writer.write_line( "_symmetry_equiv_pos_as_xyz" );
        for ( size_t i( 0 ); i != space_group_.nsymmetry_operators(); ++i )
            text_file_writer.write_line( size_t2string( i+1 ) + " " + space_group_.symmetry_operator( i ).to_string() );
        text_file_writer.write_line( "loop_" );
        text_file_writer.write_line( "_atom_site_label" ); // Needed for Materials Studio
        text_file_writer.write_line( "_atom_site_type_symbol" );
        text_file_writer.write_line( "_atom_site_fract_x" );
        text_file_writer.write_line( "_atom_site_fract_y" );
        text_file_writer.write_line( "_atom_site_fract_z" );
        size_t len( 2 );
        size_t current_size = 99;
        while ( natoms >= current_size )
        {
            ++len;
            current_size = 10 * current_size + 9;
        }
        Element hydrogen( "H" );
        for ( size_t i( 0 ); i != natoms; ++i )
        {
            if ( elements[i] == hydrogen )
                continue;
            // This is just too weird, need std::vector< DoubleWithESD > for this.
            text_file_writer.write_line( elements[i].symbol() + size_t2string( i + 1, len, '0' ) + " " +
                                         elements[i].symbol() + " " +
                                         double2string( adjust_for_translations( average_positions[ i ].average().x() ), 6 ) + " " +
                                         double2string( adjust_for_translations( average_positions[ i ].average().y() ), 6 ) + " " +
                                         double2string( adjust_for_translations( average_positions[ i ].average().z() ), 6 ) );
        }
        text_file_writer.write_line( "loop_" );
        text_file_writer.write_line( "_atom_site_aniso_label" );
        text_file_writer.write_line( "_atom_site_aniso_U_11" );
        text_file_writer.write_line( "_atom_site_aniso_U_22" );
        text_file_writer.write_line( "_atom_site_aniso_U_33" );
        text_file_writer.write_line( "_atom_site_aniso_U_12" );
        text_file_writer.write_line( "_atom_site_aniso_U_13" );
        text_file_writer.write_line( "_atom_site_aniso_U_23" );
        for ( size_t i( 0 ); i != natoms; ++i )
        {
            if ( elements[i] == hydrogen )
                continue;
            std::vector< Vector3D > cartesian_positions;
            for ( size_t j( 0 ); j != fractional_positions_trajectory[i].size(); ++j )
                cartesian_positions.push_back( crystal_lattice_average.fractional_to_orthogonal_matrix() * fractional_positions_trajectory[i][j] );
            AnisotropicDisplacementParameters adps( cartesian_positions );
            SymmetricMatrix3D Ucif = adps.U_cif( crystal_lattice_average );
            // This can go wrong if the ADP is very small, it may be printed like "1E-14" which will not be recognised in the cif.
            text_file_writer.write_line( elements[i].symbol() + size_t2string( i + 1, len, '0' ) + " " +
                                         double2string( Ucif.value( 0, 0 ), 6 ) + " " +
                                         double2string( Ucif.value( 1, 1 ), 6 ) + " " +
                                         double2string( Ucif.value( 2, 2 ), 6 ) + " " +
                                         double2string( Ucif.value( 0, 1 ), 6 ) + " " +
                                         double2string( Ucif.value( 0, 2 ), 6 ) + " " +
                                         double2string( Ucif.value( 1, 2 ), 6 ) );
        }
        text_file_writer.write_line();
        text_file_writer.write_line( "#END" );
    }
    if ( write_average_ESDs_ )
    {
        TextFileWriter text_file_writer( FileName( file_list_.base_directory(), "average_ESDs_adps", "cif" ) );
        text_file_writer.write_line( "data_avg_ESDs" );
        text_file_writer.write_line( "_symmetry_space_group_name_H-M  '" + space_group_.name() + "'" );
    //    text_file_writer.write_line( "_symmetry_Int_Tables_number     1" );
    //    text_file_writer.write_line( "_symmetry_cell_setting          triclinic" );
        text_file_writer.write_line( "_cell_length_a    " + crystallographic_style( average_a_.average(), average_a_.estimated_standard_deviation() ) );
        text_file_writer.write_line( "_cell_length_b    " + crystallographic_style( average_b_.average(), average_b_.estimated_standard_deviation() ) );
        text_file_writer.write_line( "_cell_length_c    " + crystallographic_style( average_c_.average(), average_c_.estimated_standard_deviation() ) );
        text_file_writer.write_line( "_cell_angle_alpha " + crystallographic_style( average_alpha_.average().value_in_degrees(), average_alpha_.estimated_standard_deviation().value_in_degrees() ) );
        text_file_writer.write_line( "_cell_angle_beta  " + crystallographic_style( average_beta_.average().value_in_degrees() , average_beta_.estimated_standard_deviation().value_in_degrees()  ) );
        text_file_writer.write_line( "_cell_angle_gamma " + crystallographic_style( average_gamma_.average().value_in_degrees(), average_gamma_.estimated_standard_deviation().value_in_degrees() ) );
        text_file_writer.write_line( "_cell_volume      " + crystallographic_style( average_volume_.average(), average_volume_.estimated_standard_deviation() ) );
        text_file_writer.write_line( "loop_" );
        text_file_writer.write_line( "_symmetry_equiv_pos_site_id" );
        text_file_writer.write_line( "_symmetry_equiv_pos_as_xyz" );
        for ( size_t i( 0 ); i != space_group_.nsymmetry_operators(); ++i )
            text_file_writer.write_line( size_t2string( i+1 ) + " " + space_group_.symmetry_operator( i ).to_string() );
        text_file_writer.write_line( "loop_" );
        text_file_writer.write_line( "_atom_site_label" ); // Needed for Materials Studio
        text_file_writer.write_line( "_atom_site_type_symbol" );
        text_file_writer.write_line( "_atom_site_fract_x" );
        text_file_writer.write_line( "_atom_site_fract_y" );
        text_file_writer.write_line( "_atom_site_fract_z" );
        size_t len( 2 );
        size_t current_size = 99;
        while ( natoms >= current_size )
        {
            ++len;
            current_size = 10 * current_size + 9;
        }
        for ( size_t i( 0 ); i != natoms; ++i )
        {
            // This is just too weird, need std::vector< DoubleWithESD > for this.
            text_file_writer.write_line( elements[i].symbol() + size_t2string( i + 1, len, '0' ) + " " +
                                         elements[i].symbol() + " " +
                                         crystallographic_style( adjust_for_translations( average_positions[ i ].average().x() ), average_positions[ i ].estimated_standard_deviation().x() ) + " " +
                                         crystallographic_style( adjust_for_translations( average_positions[ i ].average().y() ), average_positions[ i ].estimated_standard_deviation().y() ) + " " +
                                         crystallographic_style( adjust_for_translations( average_positions[ i ].average().z() ), average_positions[ i ].estimated_standard_deviation().z() ) );
        }
        text_file_writer.write_line( "loop_" );
        text_file_writer.write_line( "_atom_site_aniso_label" );
        text_file_writer.write_line( "_atom_site_aniso_U_11" );
        text_file_writer.write_line( "_atom_site_aniso_U_22" );
        text_file_writer.write_line( "_atom_site_aniso_U_33" );
        text_file_writer.write_line( "_atom_site_aniso_U_12" );
        text_file_writer.write_line( "_atom_site_aniso_U_13" );
        text_file_writer.write_line( "_atom_site_aniso_U_23" );
        for ( size_t i( 0 ); i != natoms; ++i )
        {
            std::vector< Vector3D > cartesian_positions;
            for ( size_t j( 0 ); j != fractional_positions_trajectory[i].size(); ++j )
                cartesian_positions.push_back( crystal_lattice_average.fractional_to_orthogonal_matrix() * fractional_positions_trajectory[i][j] );
            AnisotropicDisplacementParameters adps( cartesian_positions );
            SymmetricMatrix3D Ucif = adps.U_cif( crystal_lattice_average );
            // This can go wrong if the ADP is very small, it may be printed like "1E-14" which will not be recognised in the cif.
            text_file_writer.write_line( elements[i].symbol() + size_t2string( i + 1, len, '0' ) + " " +
                                         double2string( Ucif.value( 0, 0 ), 6 ) + " " +
                                         double2string( Ucif.value( 1, 1 ), 6 ) + " " +
                                         double2string( Ucif.value( 2, 2 ), 6 ) + " " +
                                         double2string( Ucif.value( 0, 1 ), 6 ) + " " +
                                         double2string( Ucif.value( 0, 2 ), 6 ) + " " +
                                         double2string( Ucif.value( 1, 2 ), 6 ) );
        }
        text_file_writer.write_line();
        text_file_writer.write_line( "#END" );
    }

    if ( write_sum_ )
    {
        TextFileWriter text_file_writer( FileName( file_list_.base_directory(), "average_sum", "cif" ) );
        text_file_writer.write_line( "data_sum" );
        text_file_writer.write_line( "_symmetry_space_group_name_H-M  '" + space_group_.name() + "'" );
    //    text_file_writer.write_line( "_symmetry_Int_Tables_number     1" );
    //    text_file_writer.write_line( "_symmetry_cell_setting          triclinic" );
        text_file_writer.write_line( "_cell_length_a    " + crystallographic_style( average_a_.average(), average_a_.estimated_standard_deviation() ) );
        text_file_writer.write_line( "_cell_length_b    " + crystallographic_style( average_b_.average(), average_b_.estimated_standard_deviation() ) );
        text_file_writer.write_line( "_cell_length_c    " + crystallographic_style( average_c_.average(), average_c_.estimated_standard_deviation() ) );
        text_file_writer.write_line( "_cell_angle_alpha " + crystallographic_style( average_alpha_.average().value_in_degrees(), average_alpha_.estimated_standard_deviation().value_in_degrees() ) );
        text_file_writer.write_line( "_cell_angle_beta  " + crystallographic_style( average_beta_.average().value_in_degrees() , average_beta_.estimated_standard_deviation().value_in_degrees()  ) );
        text_file_writer.write_line( "_cell_angle_gamma " + crystallographic_style( average_gamma_.average().value_in_degrees(), average_gamma_.estimated_standard_deviation().value_in_degrees() ) );
        text_file_writer.write_line( "_cell_volume      " + crystallographic_style( average_volume_.average(), average_volume_.estimated_standard_deviation() ) );
        text_file_writer.write_line( "loop_" );
        text_file_writer.write_line( "_symmetry_equiv_pos_site_id" );
        text_file_writer.write_line( "_symmetry_equiv_pos_as_xyz" );
        for ( size_t i( 0 ); i != space_group_.nsymmetry_operators(); ++i )
            text_file_writer.write_line( size_t2string( i+1 ) + " " + space_group_.symmetry_operator( i ).to_string() );
        text_file_writer.write_line( "loop_" );
        text_file_writer.write_line( "_atom_site_label" ); // Needed for Materials Studio
        text_file_writer.write_line( "_atom_site_type_symbol" );
        text_file_writer.write_line( "_atom_site_fract_x" );
        text_file_writer.write_line( "_atom_site_fract_y" );
        text_file_writer.write_line( "_atom_site_fract_z" );
        size_t len( 2 );
        size_t current_size = 99;
        while ( natoms >= current_size )
        {
            ++len;
            current_size = 10 * current_size + 9;
        }
        for ( size_t i( 0 ); i != natoms; ++i )
        {
            for ( size_t j( 0 ); j != fractional_positions_trajectory[i].size(); ++j )
            text_file_writer.write_line( elements[i].symbol() + size_t2string( i + 1, len, '0' ) + " " +
                                         elements[i].symbol() + " " +
                                         double2string_pad_plus( fractional_positions_trajectory[i][j].x(), 5, ' ' ) + " " +
                                         double2string_pad_plus( fractional_positions_trajectory[i][j].y(), 5, ' ' ) + " " +
                                         double2string_pad_plus( fractional_positions_trajectory[i][j].z(), 5, ' ' ) );
        }
        text_file_writer.write_line();
        text_file_writer.write_line( "#END" );
    }
}

// ********************************************************************************

CrystalLattice AnalyseTrajectory::average_crystal_lattice() const
{
    return CrystalLattice( average_a_.average(),
                           average_b_.average(),
                           average_c_.average(),
                           average_alpha_.average(),
                           average_beta_.average(),
                           average_gamma_.average() );
}

// ********************************************************************************

void AnalyseTrajectory::save_centres_of_mass() const
{
    TextFileWriter text_file_writer( FileName( file_list_.base_directory(), "centres_of_mass", "txt" ) );
    CrystalLattice average_crystal_lattice( this->average_crystal_lattice() );
    for ( size_t i( 0 ); i != centres_of_mass_.size(); ++i )
        text_file_writer.write_line( double2string( centres_of_mass_[i].x() ) + " " +
                                     double2string( centres_of_mass_[i].y() ) + " " +
                                     double2string( centres_of_mass_[i].z() ) + " " +
                                     double2string( (average_crystal_lattice.fractional_to_orthogonal( centres_of_mass_[i] ) - average_crystal_lattice.fractional_to_orthogonal( centres_of_mass_[0] )).length() ) );
}

// ********************************************************************************

