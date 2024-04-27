/* *********************************************
Copyright (c) 2013-2024, Cornelis Jan (Jacco) van de Streek
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

#include "InpWriter.h"
#include "CopyTextFile.h"
#include "CrystalStructure.h"
#include "FileName.h"
#include "PowderPattern.h"
#include "ReadCif.h"
#include "StringFunctions.h"
#include "TextFileReader_2.h"
#include "TextFileWriter.h"
//#include "TOPAS.h"
#include "Utilities.h"

#include <iostream>
#include <stdexcept>
#include <vector>

// ********************************************************************************

void inp_writer( const FileName & input_cif_file_name, const FileName & input_xye_file_name, const std::string & aal  )
{
    bool generate_restraints( false );
    CrystalStructure crystal_structure;
    std::cout << "Now reading cif... " + input_cif_file_name.full_name() << std::endl;
    read_cif( input_cif_file_name, crystal_structure );
    // The powder diffraction file must contain a third column with estimated standard deviations or TOPAS cannot read the file.
    // So create that file now.
    if ( FileName( input_cif_file_name.directory(), input_cif_file_name.file_name(), "xye" ).exists() )
        throw std::runtime_error( "inp_writer(): can't create " + FileName( input_cif_file_name.directory(), input_cif_file_name.file_name(), "xye" ).full_name() + " because it already exists." );
    PowderPattern powder_pattern( input_xye_file_name );
    powder_pattern.recalculate_estimated_standard_deviations();
    powder_pattern.save_xye( FileName( input_cif_file_name.directory(), input_cif_file_name.file_name(), "xye" ), false );
    TextFileWriter text_file_writer( replace_extension( input_cif_file_name, "inp" ) );
    bool file_00000001_Bonds_tsv_exists = FileName( input_cif_file_name.directory(), "00000001-Bonds", "tsv" ).exists();
    bool file_00000001_AllAngles_tsv_exists = FileName( input_cif_file_name.directory(), "00000001-AllAngles", "tsv" ).exists();
    std::vector< std::string > bond_labels_1;
    std::vector< std::string > bond_labels_2;
    std::vector< double > bond_target_values;
    std::vector< std::string > angle_labels_1;
    std::vector< std::string > angle_labels_2;
    std::vector< std::string > angle_labels_3;
    std::vector< double > angle_target_values;
    if ( file_00000001_Bonds_tsv_exists )
    {
        std::cout << "Bond restraints file found, bond restraints will be written out." << std::endl;
        TextFileReader_2 file_bond_restraints( FileName( input_cif_file_name.directory(), "00000001-Bonds", "tsv" ) );
        std::vector< std::string > words;
        for ( size_t i( 1 ); i != file_bond_restraints.size(); ++i )
        {
            words = split( file_bond_restraints.line( i ) );
            if ( words.size() != 4 )
                throw std::runtime_error( "inp_writer(): unexpected format in bond restraints file: >" + file_bond_restraints.line( i ) + "<" );
            bond_labels_1.push_back( words[1]+aal );
            bond_labels_2.push_back( words[2]+aal );
            bond_target_values.push_back( string2double( words[3] ) );
        }
    }
    else if ( generate_restraints )
    {
        std::vector< std::string > temp_labels_1;
        std::vector< std::string > temp_labels_2;
        crystal_structure.list_all_bonds( temp_labels_1, temp_labels_2, bond_target_values );
        for ( size_t i( 0 ); i != bond_target_values.size(); ++i )
        {
            bond_labels_1.push_back( temp_labels_1[i]+aal );
            bond_labels_2.push_back( temp_labels_2[i]+aal );
        }
    }
    else
        std::cout << "No bond restraints file found, no bond restraints will be written out." << std::endl;
    if ( file_00000001_AllAngles_tsv_exists )
    {
        std::cout << "Angle restraints file found, angle restraints will be written out." << std::endl;
        TextFileReader_2 file_angle_restraints( FileName( input_cif_file_name.directory(), "00000001-AllAngles", "tsv" ) );
        std::vector< std::string > words;
        for ( size_t i( 1 ); i != file_angle_restraints.size(); ++i )
        {
            words = split( file_angle_restraints.line( i ) );
            if ( words.size() != 5 )
                throw std::runtime_error( ".inp writer: unexpected format in angle restraints file: >" + file_angle_restraints.line( i ) + "<" );
            angle_labels_1.push_back( words[1]+aal );
            angle_labels_2.push_back( words[2]+aal );
            angle_labels_3.push_back( words[3]+aal );
            angle_target_values.push_back( string2double( words[4] ) );
        }
    }
    else if ( generate_restraints )
    {
        std::vector< std::string > temp_labels_1;
        std::vector< std::string > temp_labels_2;
        std::vector< std::string > temp_labels_3;
        crystal_structure.list_all_angles( temp_labels_1, temp_labels_2, temp_labels_3, angle_target_values );
        for ( size_t i( 0 ); i != angle_target_values.size(); ++i )
        {
            angle_labels_1.push_back( temp_labels_1[i]+aal );
            angle_labels_2.push_back( temp_labels_2[i]+aal );
            angle_labels_3.push_back( temp_labels_3[i]+aal );
        }
    }
    else
        std::cout << "No angle restraints file found, no angle restraints will be written out." << std::endl;
    // We must find the unique combinations, e.g.:
    // angle r1-r2-r3 and bond r3-r2 only require two combinations: r1-r2 and r3-r2
    // Create list of unique bonds
//    std::vector< std::string > unique_labels_1;
//    std::vector< std::string > unique_labels_2;
//    for ( size_t i( 0 ); i != angle_labels_1.size(); ++i )
//    {
//        bool found( false );
//        // Find bond 1
//        for ( size_t j( 0 ); j != unique_labels_1.size(); ++j )
//        {
//            if ( ( unique_labels_1[j] == angle_labels_1[i] ) && ( unique_labels_2[j] == angle_labels_2[i] ) )
//            {
//                 found = true;
//                 break;
//            }
//        }
//        if ( ! found )
//        {
//            unique_labels_1.push_back( angle_labels_1[i] );
//            unique_labels_2.push_back( angle_labels_2[i] );
//        }
//        found = false;
//        // Find bond 2
//        for ( size_t j( 0 ); j != unique_labels_1.size(); ++j )
//        {
//            if ( ( unique_labels_1[j] == angle_labels_3[i] ) && ( unique_labels_2[j] == angle_labels_2[i] ) )
//            {
//                 found = true;
//                 break;
//            }
//        }
//        if ( ! found )
//        {
//            unique_labels_1.push_back( angle_labels_3[i] );
//            unique_labels_2.push_back( angle_labels_2[i] );
//        }
//    }
    text_file_writer.write_line( "penalties_weighting_K1 5" );
    text_file_writer.write_line( "'do_errors_include_penalties" );
    text_file_writer.write_line( "prm JvdS_shift = Get(refine_ls_shift_on_su_max); : 0.0" );
    text_file_writer.write_line( "prm JvdS_numpar = Get(number_independent_parameters); :  0.0" );
    text_file_writer.write_line( "r_exp  0.0" );
    text_file_writer.write_line( "r_exp_dash  0.0" );
    text_file_writer.write_line( "r_wp  0.0" );
    text_file_writer.write_line( "r_wp_dash  0.0" );
    text_file_writer.write_line( "r_p  0.0" );
    text_file_writer.write_line( "r_p_dash  0.0" );
    text_file_writer.write_line( "gof  0.0" );
    text_file_writer.write_line( "'continue_after_convergence" );
    text_file_writer.write_line( "xdd " + FileName( input_cif_file_name.directory(), input_cif_file_name.file_name(), "xye" ).full_name() + " xye_format" );
    text_file_writer.write_line( "  bkg @" );
    for ( size_t i( 0 ); i != 20; ++i )
        text_file_writer.write_line( "    0.0" );
    text_file_writer.write_line( "  start_X       " + double2string( powder_pattern.two_theta_start().value_in_degrees() ) );
    text_file_writer.write_line( "  finish_X      " + double2string( powder_pattern.two_theta_end().value_in_degrees() ) );
    text_file_writer.write_line( "  x_calculation_step " + double2string( powder_pattern.average_two_theta_step().value_in_degrees() ) );
    text_file_writer.write_line( "'  Specimen_Displacement(@ , 0.0 )" );
    text_file_writer.write_line( "'  Absorption(@ , 0,0 )" );
    text_file_writer.write_line( "  Zero_Error(@ , 0.02 )" );
    text_file_writer.write_line( "'Synchrotron use: LP_Factor( 90.0 )" );
    text_file_writer.write_line( "'Neutrons use: LP_Factor( 90.0 )" );
    text_file_writer.write_line( "'No monochromator use: LP_Factor( 0.0 )" );
    text_file_writer.write_line( "'Ge Monochromator, Cu radiation, use LP_Factor( 27.3 )" );
    text_file_writer.write_line( "'Graphite Monochromator, Cu radiation, use LP_Factor( 26.4 )" );
    text_file_writer.write_line( "'Quartz Monochromator, Cu radiation, use LP_Factor( 26.6 )" );
    text_file_writer.write_line( "  LP_Factor( 26.5 )" );
    text_file_writer.write_line( "'  Variable_Divergence(@ , 30.0 )" );
    text_file_writer.write_line( "  axial_conv" );
    text_file_writer.write_line( "    filament_length @ 4.97890" );
    text_file_writer.write_line( "    sample_length @ 2.43658" );
    text_file_writer.write_line( "    receiving_slit_length @ 5.25714" );
    text_file_writer.write_line( "    axial_n_beta 50" );
    text_file_writer.write_line( "  lam" );
    text_file_writer.write_line( "    ymin_on_ymax 0.001" );
    text_file_writer.write_line( "    la 1 lo  1.540560" );
    text_file_writer.write_line( "  str" );
    text_file_writer.write_line( "    r_bragg  0.0" );
    text_file_writer.write_line( "    CS_G(@ , 107.03272`)" );
    text_file_writer.write_line( "    CS_L(@ , 9999.99881`)" );
    text_file_writer.write_line( "    Strain_G(@ , 0.49554`)" );
    text_file_writer.write_line( "    Strain_L(@ , 0.03347`)" );
    text_file_writer.write_line( "    prm  sh_scale_l" + aal + " " + double2string( powder_pattern.cumulative_intensity() * 4.2E-10 ) );
    text_file_writer.write_line( "    spherical_harmonics_hkl sh_l" + aal );
    text_file_writer.write_line( "      sh_order 6" );
    text_file_writer.write_line( "    lor_fwhm = Abs( sh_scale_l" + aal + " * sh_l" + aal + " );" );
    text_file_writer.write_line( "    prm  sh_scale_g" + aal + "  0.01948" );
    text_file_writer.write_line( "    spherical_harmonics_hkl sh_g" + aal );
    text_file_writer.write_line( "      sh_order 6" );
    text_file_writer.write_line( "    gauss_fwhm = Abs( sh_scale_g" + aal + " * sh_g" + aal + " );" );
    if ( crystal_structure.space_group().crystal_system() == "triclinic" )
    {
        text_file_writer.write_line( "    a  @ " + double2string( crystal_structure.crystal_lattice().a() ) );
        text_file_writer.write_line( "    b  @ " + double2string( crystal_structure.crystal_lattice().b() ) );
        text_file_writer.write_line( "    c  @ " + double2string( crystal_structure.crystal_lattice().c() ) );
        text_file_writer.write_line( "    al @ " + double2string( crystal_structure.crystal_lattice().alpha().value_in_degrees() ) );
        text_file_writer.write_line( "    be @ " + double2string( crystal_structure.crystal_lattice().beta().value_in_degrees() ) );
        text_file_writer.write_line( "    ga @ " + double2string( crystal_structure.crystal_lattice().gamma().value_in_degrees() ) );
    }
    else if ( crystal_structure.space_group().crystal_system() == "monoclinic" )
    {
        text_file_writer.write_line( "    a  @ " + double2string( crystal_structure.crystal_lattice().a() ) );
        text_file_writer.write_line( "    b  @ " + double2string( crystal_structure.crystal_lattice().b() ) );
        text_file_writer.write_line( "    c  @ " + double2string( crystal_structure.crystal_lattice().c() ) );
        text_file_writer.write_line( "    al   90.0" );
        text_file_writer.write_line( "    be @ " + double2string( crystal_structure.crystal_lattice().beta().value_in_degrees() ) );
        text_file_writer.write_line( "    ga   90.0" );
    }
    else if ( crystal_structure.space_group().crystal_system() == "orthorhombic" )
    {
        text_file_writer.write_line( "    a  @ " + double2string( crystal_structure.crystal_lattice().a() ) );
        text_file_writer.write_line( "    b  @ " + double2string( crystal_structure.crystal_lattice().b() ) );
        text_file_writer.write_line( "    c  @ " + double2string( crystal_structure.crystal_lattice().c() ) );
        text_file_writer.write_line( "    al 90.0" );
        text_file_writer.write_line( "    be 90.0" );
        text_file_writer.write_line( "    ga 90.0" );
    }
    else if ( crystal_structure.space_group().crystal_system() == "tetragonal" )
    {
        text_file_writer.write_line( "    prm uc_a " + double2string( crystal_structure.crystal_lattice().a() ) );
        text_file_writer.write_line( "    a = uc_a;" );
        text_file_writer.write_line( "    b = uc_a;" );
        text_file_writer.write_line( "    c  @ " + double2string( crystal_structure.crystal_lattice().c() ) );
        text_file_writer.write_line( "    al 90.0" );
        text_file_writer.write_line( "    be 90.0" );
        text_file_writer.write_line( "    ga 90.0" );
    }
    else if ( crystal_structure.space_group().crystal_system() == "cubic" )
    {
        text_file_writer.write_line( "    prm uc_a " + double2string( crystal_structure.crystal_lattice().a() ) );
        text_file_writer.write_line( "    a = uc_a;" );
        text_file_writer.write_line( "    b = uc_a;" );
        text_file_writer.write_line( "    c = uc_a;" );
        text_file_writer.write_line( "    al 90.0" );
        text_file_writer.write_line( "    be 90.0" );
        text_file_writer.write_line( "    ga 90.0" );
    }
    else if ( crystal_structure.space_group().crystal_system() == "hexagonal" )
    {
        text_file_writer.write_line( "    prm uc_a " + double2string( crystal_structure.crystal_lattice().a() ) );
        text_file_writer.write_line( "    a = uc_a;" );
        text_file_writer.write_line( "    b = uc_a;" );
        text_file_writer.write_line( "    c  @ " + double2string( crystal_structure.crystal_lattice().c() ) );
        text_file_writer.write_line( "    al 90.0" );
        text_file_writer.write_line( "    be 90.0" );
        text_file_writer.write_line( "    ga 120.0" );
    }
    else if ( crystal_structure.space_group().crystal_system() == "trigonal" )
    {
        // Two options: trigonal or rhombohedral
        if ( nearly_equal( crystal_structure.crystal_lattice().alpha(), Angle::angle_90_degrees() ) && nearly_equal( crystal_structure.crystal_lattice().beta(), Angle::angle_120_degrees() ) )
        {
            text_file_writer.write_line( "    prm uc_a " + double2string( crystal_structure.crystal_lattice().a() ) );
            text_file_writer.write_line( "    a = uc_a;" );
            text_file_writer.write_line( "    b = uc_a;" );
            text_file_writer.write_line( "    c  @ " + double2string( crystal_structure.crystal_lattice().c() ) );
            text_file_writer.write_line( "    al 90.0" );
            text_file_writer.write_line( "    be 90.0" );
            text_file_writer.write_line( "    ga 120.0" );
        }
        else
        {
            text_file_writer.write_line( "    prm uc_a " + double2string( crystal_structure.crystal_lattice().a() ) );
            text_file_writer.write_line( "    a = uc_a;" );
            text_file_writer.write_line( "    b = uc_a;" );
            text_file_writer.write_line( "    c = uc_a;" );
            text_file_writer.write_line( "    prm uc_alpha " + double2string( crystal_structure.crystal_lattice().alpha().value_in_degrees() ) );
            text_file_writer.write_line( "    al = uc_alpha;" );
            text_file_writer.write_line( "    be = uc_alpha;" );
            text_file_writer.write_line( "    ga = uc_alpha;" );
        }
    }
    else
        throw std::runtime_error( "inp_writer(): we should never be here." );
    text_file_writer.write_line( "    MVW( 0.0, 0.0, 0.0 )");
    text_file_writer.write_line( "    space_group \"" + remove( remove( crystal_structure.space_group().name(), '_' ), ' ') + "\"");
    text_file_writer.write_line( "    scale @  0.0001" );
    text_file_writer.write_line( "'    PO(@ , 1.0, , 1 0 0 )" );
    text_file_writer.write_line( "'    PO_Spherical_Harmonics( sh, 6 )" );
    text_file_writer.write_line( "    macro ref_flag" + aal + " { @ }" );
    text_file_writer.write_line( "    prm    bnonh" + aal + " 3.0" );
    text_file_writer.write_line( "    prm bh" + aal + " = 1.2 * bnonh" + aal + ";" );
    for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
    {
        text_file_writer.write( "    site " + crystal_structure.atom( i ).label() + aal + " x ref_flag" + aal + " " + double2string_pad_plus( crystal_structure.atom( i ).position().x(), 5, ' ' ) +
                                                                                          " y ref_flag" + aal + " " + double2string_pad_plus( crystal_structure.atom( i ).position().y(), 5, ' ' ) +
                                                                                          " z ref_flag" + aal + " " + double2string_pad_plus( crystal_structure.atom( i ).position().z(), 5, ' ' ) +
                                     " occ " + pad( crystal_structure.atom( i ).element().symbol(), 2, ' ' ) + " 1 beq = " );
        if ( crystal_structure.atom( i ).element().is_H_or_D() )
            text_file_writer.write_line( "bh" + aal + ";" );
        else
            text_file_writer.write_line( "bnonh" + aal + ";" );
    }
    if ( file_00000001_Bonds_tsv_exists )
    {
        text_file_writer.write_line( "    prm !bond_width   0" );
        text_file_writer.write_line( "    prm !bond_weight  10000" );
        for ( size_t i( 0 ); i != bond_labels_1.size(); ++i )
        {
            if ( element_from_atom_label( bond_labels_1[i] ).is_H_or_D() || element_from_atom_label( bond_labels_2[i] ).is_H_or_D() )
                text_file_writer.write_line( "    Distance_Restrain( " + bond_labels_1[i] + " " + bond_labels_2[i] + ", 0.95, 0.0, bond_width, bond_weight )" );
            else
                text_file_writer.write_line( "    Distance_Restrain( " + bond_labels_1[i] + " " + bond_labels_2[i] + ", " + double2string( bond_target_values[i] ) + ", 0.0, bond_width, bond_weight )" );
        }
    }
    if ( file_00000001_AllAngles_tsv_exists )
    {
        text_file_writer.write_line( "    prm !angle_width  0" );
        text_file_writer.write_line( "    prm !angle_weight 1" );
        for ( size_t i( 0 ); i != angle_labels_1.size(); ++i )
        {
           text_file_writer.write_line( "    Angle_Restrain( " + angle_labels_1[i] + " " + angle_labels_2[i] + " " + angle_labels_3[i] + ", " + double2string( angle_target_values[i] ) + ", 0.0, angle_width, angle_weight )" );
        }
    }
    text_file_writer.write_line( "    prm !flatten_width 0" );
    text_file_writer.write_line( "    prm !flatten_weight    100000" );
    text_file_writer.write_line( "    'Flatten( C19 N11 N35 C47 H54 C50 H58 C45 H49 C34 C46, , 0.0, flatten_width, flatten_weight )" );
    text_file_writer.write_line( "    Out_CIF_STR( " + FileName( input_cif_file_name.directory(), input_cif_file_name.file_name() + "_RR", "cif" ).full_name() + " )" );
    text_file_writer.write_line( "  xdd_out " + FileName( input_cif_file_name.directory(), input_cif_file_name.file_name() + "_profile", "txt" ).full_name() + " load out_record out_fmt out_eqn" );
    text_file_writer.write_line( "  {" );
    text_file_writer.write_line( "      \" %11.5f \" = X;" );
    text_file_writer.write_line( "      \" %11.5f \" = Yobs;" );
    text_file_writer.write_line( "      \" %11.5f \" = Ycalc;" );
    text_file_writer.write_line( "      \" %11.5f\\n\" = SigmaYobs;" );
    text_file_writer.write_line( "  }" );
    text_file_writer.write_line( "  phase_out " + FileName( input_cif_file_name.directory(), input_cif_file_name.file_name() + "_tickmarks", "txt" ).full_name() + " load out_record out_fmt out_eqn" );
    text_file_writer.write_line( "  {" );
    text_file_writer.write_line( "      \" %11.5f -200\\n\" = 2.0 * Rad * Th;" );
    text_file_writer.write_line( "  }" );
    text_file_writer.write_line( "  ' Structure factors should be in *.fcf, intensities should be in *.hkl." );
    text_file_writer.write_line( "  out " + FileName( input_cif_file_name.directory(), input_cif_file_name.file_name(), "fcf" ).full_name() );
    text_file_writer.write_line( "  Out_String(\"\\ndata_\")" );
    text_file_writer.write_line( "  Out(Get(a), \"\\n_cell_length_a %V\")" );
    text_file_writer.write_line( "  Out(Get(b), \"\\n_cell_length_b %V\")" );
    text_file_writer.write_line( "  Out(Get(c), \"\\n_cell_length_c %V\")" );
    text_file_writer.write_line( "  Out(Get(al), \"\\n_cell_angle_alpha %V\")" );
    text_file_writer.write_line( "  Out(Get(be), \"\\n_cell_angle_beta  %V\")" );
    text_file_writer.write_line( "  Out(Get(ga), \"\\n_cell_angle_gamma %V\")" );
    text_file_writer.write_line( "  Out_String(\"\\n_shelx_F_squared_multiplier 1\")" );
    text_file_writer.write_line( "  Out_String(\"\\nloop_\\n_symmetry_equiv_pos_as_xyz\")" );
    text_file_writer.write_line( "  Out(Get(sp_xyzs_txt), \"%s\")" );
    text_file_writer.write_line( "  Out_String(\"\\nloop_\")" );
    text_file_writer.write_line( "  Out_String(\"\\n_refln_index_h\")" );
    text_file_writer.write_line( "  Out_String(\"\\n_refln_index_k\")" );
    text_file_writer.write_line( "  Out_String(\"\\n_refln_index_l\")" );
    text_file_writer.write_line( "  Out_String(\"\\n_refln_F_squared_calc\")" );
    text_file_writer.write_line( "  Out_String(\"\\n_refln_F_squared_meas\")" );
    text_file_writer.write_line( "  Out_String(\"\\n_refln_F_squared_sigma\")" );
    text_file_writer.write_line( "  Out_String(\"\\n_refln_observed_status\")" );
    text_file_writer.write_line( "  phase_out " + FileName( input_cif_file_name.directory(), input_cif_file_name.file_name(), "fcf" ).full_name() + " append" );
    text_file_writer.write_line( "    load out_record out_fmt out_eqn" );
    text_file_writer.write_line( "    {" );
    text_file_writer.write_line( "      \"\\n%4.0f\" = H;" );
    text_file_writer.write_line( "      \" %4.0f\" = K;" );
    text_file_writer.write_line( "      \" %4.0f\" = L;" );
    text_file_writer.write_line( "      \" %12.2f\" = I_no_scale_pks / ( Get( scale ) * M ); ' This is the calculated F^2, same as F2_Merged" );
    text_file_writer.write_line( "      \" %12.2f\" = Iobs_no_scale_pks / ( Get( scale ) * M ); ' M is the multiplicity" );
    text_file_writer.write_line( "      \" %9.3f o\" = Iobs_no_scale_pks_err / ( Get( scale ) * M );" );
    text_file_writer.write_line( "    }" );
    copy_text_file( replace_extension( input_cif_file_name, "inp" ), replace_extension( input_cif_file_name, "org" ) );
    std::cout << "DO NOT FORGET TO RESET THE SPACE GROUP." << std::endl;
}

// ********************************************************************************

