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

#include "GeneratePowderCIF.h"
#include "3DCalculations.h"
#include "BasicMathsFunctions.h"
#include "DoubleWithESD.h"
#include "Element.h"
#include "FileName.h"
#include "PhysicalConstants.h"
#include "PowderPatternCalculator.h"
#include "String2Fraction.h"
#include "StringFunctions.h"
#include "SymmetryOperator.h"
#include "Utilities.h"

#include <iostream>
#include <set>
#include <stdexcept>

// ********************************************************************************

GeneratePowderCIF::GeneratePowderCIF( const std::string & directory, const std::string & base_name ) :
directory_(directory),
base_name_(base_name),
output_file_( FileName( directory_, base_name_ + "_final", "cif" ) ),
output_file_xrpd_( FileName( directory_, base_name_, "pcif" ) ),
output_file_R_( FileName( directory_, base_name_ + "_plot", "R" ) ),
replace_hydrogen_atoms_(false),
relabel_(false),
padding_length_(31),
longest_label_size_(0),
longest_x_(0),
longest_y_(0),
longest_z_(0)
{
    file_name_pro_ = base_name_ + "_profile.txt";
    file_name_tic_ = base_name_ + "_tickmarks.txt";
    if ( FileName( directory_, base_name_, "xye" ).exists() )
        powder_pattern_ = PowderPattern( FileName( directory_, base_name_, "xye" ) );
    else if ( FileName( directory_, "temp", "xye" ).exists() )
        powder_pattern_ = PowderPattern( FileName( directory_, "temp", "xye" ) );
    else
        throw std::runtime_error( "GeneratePowderCIF::GeneratePowderCIF(): .xye file not found." );
}

// ********************************************************************************

void GeneratePowderCIF::generate()
{
    // TODO @@ remove duplicates from bond list and angle list
    // Needless to say, atom labels must be consistent
    file_cif_.read_file( FileName( directory_, base_name_, "cif" ) );
    file_inp_.read_file( FileName( directory_, base_name_, "inp" ) );
    replace_hydrogen_atoms_ = FileName( directory_, base_name_ + "_Hmi", "cif" ).exists();
    relabel_ = FileName( directory_, base_name_ + "_relabel", "txt" ).exists();
    TextFileReader_2 file_pro( FileName( directory_, base_name_ + "_profile", "txt" ) );
    file_ext_.read_file( FileName( directory_, base_name_ + "_additional_info", "txt" ) );
    file_bond_lengths_.read_file( FileName( directory_, base_name_ + "-Bonds", "tsv" ) );
    file_valence_angles_.read_file( FileName( directory_, base_name_ + "-AllAngles", "tsv" ) );
    std::vector< std::string > words;
    size_t iLine;
    Splitter splitter( "," );
    splitter.set_merge_delimiters( false );
    iLine = file_cif_.find( "_atom_site_U_iso_or_equiv" );
    for ( size_t i( iLine + 1 ); i != file_cif_.size(); ++i )
    {
//C1 C   2 0.2611(2) -0.22243(12) 0.1923(12) 1 4.3(3)     0.05426
        words = split( file_cif_.line( i ) );
        if ( words.size() == 0 )
            continue;
         chemical_formula_.add_element( Element( words[1] ) );
       if ( words[3].length() > longest_x_ )
            longest_x_ = words[3].length();
        if ( words[4].length() > longest_y_ )
            longest_y_ = words[4].length();
        if ( words[5].length() > longest_z_ )
            longest_z_ = words[5].length();
        if ( relabel_ )
            continue;
        if ( words[0].length() > longest_label_size_ )
            longest_label_size_ = words[0].length();
    }
    // Bit of a fudge in case the longest does NOT have a minus sign (which then gets padded as well).
    ++longest_x_;
    ++longest_y_;
    ++longest_z_;
    std::vector< std::string > old_label;
    std::vector< std::string > new_label;
    if ( relabel_ )
    {
        std::cout << "Relabel file found, atoms will be relabelled." << std::endl;
        TextFileReader_2 file_relabel( FileName( directory_, base_name_ + "_relabel", "txt" ) );
        for ( size_t i( 0 ); i != file_relabel.size(); ++i )
        {
            words = split( file_relabel.line( i ) );
            if ( words.size() != 2 )
                throw std::runtime_error( "GeneratePowderCIF::generate(): unexpected format in relabelling file: >" + file_relabel.line( i ) + "<" );
            if ( labels_.find( words[0] ) != labels_.end() )
                throw std::runtime_error( "GeneratePowderCIF::generate(): relabelling file contains duplicate label: >" + words[0] + "<" );
            if ( element_from_atom_label( words[0] ) != element_from_atom_label( words[1] ) )
                throw std::runtime_error( "GeneratePowderCIF::generate(): relabelling would change an element: >" + words[0] + " != " + words[1] + "<" );
            labels_[ words[0] ] = words[1];
            if ( words[1].length() > longest_label_size_ )
                longest_label_size_ = words[1].length();
        }
    }
    else
        std::cout << "No relabel file found, atoms will not be relabelled." << std::endl;
    if ( replace_hydrogen_atoms_ )
        std::cout << "File with new H atom positions found, H atoms will be replaced." << std::endl;
    else
        std::cout << "No file with new H atom positions found, H atoms are assumed to have been refined." << std::endl;
//  lam
//    ymin_on_ymax 0.001
//    la 0.6666667 lo 1.54060
//    la 0.3333333 lo 1.54439
    iLine = file_inp_.find( "ymin_on_ymax" );
    words = split( file_inp_.line( iLine+1 ) );
    double wavelength_1 = string2double( words[3] );
    std::string ratio_wavelength_1 = words[1];
    words = split( file_inp_.line( iLine+2 ) );
    if ( ( words.size() == 4 ) && ( words[0] == "la" ) )
    {
        double wavelength_2 = string2double( words[3] );
        wavelength_ = Wavelength( wavelength_1, wavelength_2, string2double( words[1] ) / string2double( ratio_wavelength_1 ) );
    }
    else
        wavelength_ = Wavelength( wavelength_1 );
    iLine = file_ext_.find( "Z_prime" );
    std::string line = file_ext_.line( iLine );
    size_t iPos = line.find( "Z_prime" );
    Fraction Z_prime = string2Fraction( line.substr( iPos + 8 ) );
    Fraction Z_prime_inverse = Z_prime;
    Z_prime_inverse.reciprocal();
    chemical_formula_.multiply( Z_prime_inverse );
    output_file_.write_line( "data_global" );
    output_file_.write_line();
    output_file_.write_line( "#==============================================================================" );
    output_file_.write_line();
    output_file_.write_line( "# 1. SUBMISSION DETAILS" );
    output_file_.write_line();
    insert_keyword_and_value( "_audit_creation_method", "'manual editing'" );
    insert_keyword_and_value( "_publ_contact_author_name", "'Jacco van de Streek'" );
    output_file_.write_line( "_publ_contact_author_address" );
    output_file_.write_line( ";" );
    output_file_.write_line( "Institute for Inorganic and Analytical Chemistry" );
    output_file_.write_line( "Goethe-University Frankfurt" );
    output_file_.write_line( "Max-von-Laue-Str. 7" );
    output_file_.write_line( "Frankfurt D-60438" );
    output_file_.write_line( "Germany" );
    output_file_.write_line( ";" );
    insert_keyword_and_value( "_publ_contact_author_email", "'jaccovandestreek@yahoo.co.uk'" );
    insert_keyword_and_value( "_publ_contact_author_fax", "?" );
    insert_keyword_and_value( "_publ_contact_author_phone", "'+49-69-79829178'" );
    output_file_.write_line();
    output_file_.write_line( "#==============================================================================" );
    output_file_.write_line();
    output_file_.write_line( "# 2. TITLE AND AUTHOR LIST" );
    output_file_.write_line();
    output_file_.write_line( "_publ_section_title" );
    output_file_.write_line( ";" );
    output_file_.write_line( "@@@@" );
    output_file_.write_line( ";" );
    output_file_.write_line( "loop_" );
    output_file_.write_line( "    _publ_author_name" );
    output_file_.write_line( "    _publ_author_footnote" );
    output_file_.write_line( "    _publ_author_address" );
    output_file_.write_line( "    'van de Streek, Jacco' ." );
    output_file_.write_line( ";" );
    output_file_.write_line( "Institute for Inorganic and Analytical Chemistry" );
    output_file_.write_line( "Goethe-University Frankfurt" );
    output_file_.write_line( "Max-von-Laue-Str. 7" );
    output_file_.write_line( "Frankfurt D-60438" );
    output_file_.write_line( "Germany" );
    output_file_.write_line( ";" );
    output_file_.write_line( "data_I" );
    insert_keyword_and_value( "_database_code_depnum_ccdc_archive", "'CCDC @@@@'" );
    output_file_.write_line();
    output_file_.write_line( "#==============================================================================" );
    output_file_.write_line();
    output_file_.write_line( "# 3. CHEMICAL DATA" );
    output_file_.write_line();
    output_file_.write_line( "_chemical_name_systematic" );
    output_file_.write_line( ";" );
    output_file_.write_line( "@@@@" );
    output_file_.write_line( ";" );
    insert_keyword( "_chemical_name_common" );
    insert_keyword_and_value( "_chemical_formula_moiety", "'" + chemical_formula_.to_string( true ) + "'" );
    insert_keyword_and_value( "_chemical_formula_sum", "'" + chemical_formula_.to_string( true ) + "'" );
    insert_keyword_and_value( "_chemical_formula_weight", double2string( chemical_formula_.molecular_weight() ) );
    output_file_.write_line();
    output_file_.write_line( "#==============================================================================" );
    output_file_.write_line();
    output_file_.write_line( "# 4. POWDER SPECIMEN AND CRYSTAL DATA" );
    output_file_.write_line();
//_space_group P-1
//loop_
//_symmetry_equiv_pos_as_xyz
//	 'x, y, z '
//	 '-x, -y, -z '
    iLine = file_cif_.find( "_cell_length_a" );
    words = split( file_cif_.line( iLine ) );
    DoubleWithESD a( words[1] );
    iLine = file_cif_.find( "_cell_length_b" );
    words = split( file_cif_.line( iLine ) );
    DoubleWithESD b( words[1] );
    iLine = file_cif_.find( "_cell_length_c" );
    words = split( file_cif_.line( iLine ) );
    DoubleWithESD c( words[1] );
    iLine = file_cif_.find( "_cell_angle_alpha" );
    words = split( file_cif_.line( iLine ) );
    DoubleWithESD alpha( words[1] );
    iLine = file_cif_.find( "_cell_angle_beta" );
    words = split( file_cif_.line( iLine ) );
    DoubleWithESD beta( words[1] );
    iLine = file_cif_.find( "_cell_angle_gamma" );
    words = split( file_cif_.line( iLine ) );
    DoubleWithESD gamma( words[1] );
    crystal_lattice_ = CrystalLattice( a.value(),
                                       b.value(),
                                       c.value(),
                                       Angle::from_degrees( alpha.value() ),
                                       Angle::from_degrees( beta.value() ),
                                       Angle::from_degrees( gamma.value() ) );
    insert_keyword_and_value( "_space_group_crystal_system", LatticeSystem2string( crystal_lattice_.lattice_system() ) );
    insert_keyword( "_space_group_name_H-M_alt" );
    insert_keyword( "_space_group_name_Hall" );
    output_file_.write_line( "loop_" );
    output_file_.write_line( "    _space_group_symop_operation_xyz" );
    iLine = file_cif_.find( "_symmetry_equiv_pos_as_xyz" );
    ++iLine;
    std::vector< SymmetryOperator > symmetry_operators;
    while ( file_cif_.line( iLine ) != "loop_" )
    {
        SymmetryOperator symmetry_operator( extract_delimited_text( file_cif_.line( iLine ), "'", "'" ) );
        output_file_.write_line( "    '" + symmetry_operator.to_string() + "'" );
        symmetry_operators.push_back( symmetry_operator );
        ++iLine;
    }
    Fraction Z = symmetry_operators.size() * Z_prime;
    if ( ! Z.is_integer() )
        std::cout << "Warning: Z is not an integer, Z = " << Z.to_string() << std::endl;
    space_group_ = SpaceGroup( symmetry_operators, "P21/c");
    crystal_structure_.set_crystal_lattice( crystal_lattice_ );
    crystal_structure_.set_space_group( space_group_ );
    PowderPatternCalculator powder_pattern_calculator( crystal_structure_ );
    powder_pattern_calculator.set_wavelength( wavelength_ );
    iLine = file_cif_.find( "_cell_length_a" );
    words = split( file_cif_.line( iLine ) );
    insert_keyword_and_value( "_cell_length_a", words[1] );
    iLine = file_cif_.find( "_cell_length_b" );
    words = split( file_cif_.line( iLine ) );
    insert_keyword_and_value( "_cell_length_b", words[1] );
    iLine = file_cif_.find( "_cell_length_c" );
    words = split( file_cif_.line( iLine ) );
    insert_keyword_and_value( "_cell_length_c", words[1] );
    iLine = file_cif_.find( "_cell_angle_alpha" );
    words = split( file_cif_.line( iLine ) );
    insert_keyword_and_value( "_cell_angle_alpha", words[1] );
    iLine = file_cif_.find( "_cell_angle_beta" );
    words = split( file_cif_.line( iLine ) );
    insert_keyword_and_value( "_cell_angle_beta", words[1] );
    iLine = file_cif_.find( "_cell_angle_gamma" );
    words = split( file_cif_.line( iLine ) );
    insert_keyword_and_value( "_cell_angle_gamma", words[1] );
    iLine = file_cif_.find( "_cell_volume" );
    words = split( file_cif_.line( iLine ) );
    insert_keyword_and_value( "_cell_volume", words[1] );
    insert_keyword_and_value( "_cell_formula_units_Z", Z.to_string() );
    insert_keyword( "_cell_measurement_temperature" );
    insert_keyword_and_value( "_exptl_crystal_density_diffrn", double2string( ( ( chemical_formula_.molecular_weight() * Z ) / crystal_lattice_.volume() ) / ( Avogadros_constant / 1.0E24 ), 3 ) );
    insert_keyword_and_value( "_exptl_crystal_F_000", double2string( chemical_formula_.nelectrons() * Z.to_double() ) );
    insert_keyword( "_pd_char_colour" );
    insert_keyword_and_value( "_pd_char_particle_morphology", "powder" );
    insert_keyword( "_exptl_absorpt_coefficient_mu" );
    output_file_.write_line();
    output_file_.write_line( "#==============================================================================" );
    output_file_.write_line();
    output_file_.write_line( "# 5. EXPERIMENTAL DATA" );
    output_file_.write_line();
    insert_keyword_and_value( "_exptl_special_details", "?" );
    insert_keyword( "_diffrn_ambient_temperature" );
    iLine = file_ext_.find( "_diffrn_radiation_probe" );
    if ( iLine == std::string::npos )
        insert_keyword_and_value( "_diffrn_radiation_probe", "x-ray" );
    else
        output_file_.write_line( file_ext_.line( iLine ) );
    insert_keyword_and_value( "_diffrn_radiation_type", wavelength_.cif_style() );
    insert_keyword_and_value( "_diffrn_radiation_wavelength", double2string( wavelength_.average_wavelength() ) );
    insert_keyword( "_diffrn_radiation_monochromator" );
    insert_keyword( "_diffrn_measurement_device_type" );
    insert_keyword_and_value( "_pd_meas_number_of_points", size_t2string( powder_pattern_.size() ) );
    insert_keyword_and_value( "_pd_meas_2theta_range_min", double2string( powder_pattern_.two_theta_start().value_in_degrees() ) );
    insert_keyword_and_value( "_pd_meas_2theta_range_max", double2string( powder_pattern_.two_theta_end().value_in_degrees() ) );
    insert_keyword_and_value( "_pd_meas_2theta_range_inc", double2string( powder_pattern_.average_two_theta_step().value_in_degrees() ) );
    insert_keyword( "_pd_spec_mount_mode" );
    insert_keyword( "_pd_spec_mounting" );
    insert_keyword( "_pd_spec_shape" );
    insert_keyword( "_pd_spec_size_axial" );
    insert_keyword( "_pd_spec_size_equat" );
    insert_keyword( "_pd_spec_size_thick" );
    output_file_.write_line();
    output_file_.write_line( "#==============================================================================" );
    output_file_.write_line();
    output_file_.write_line( "# 6. REFINEMENT DATA" );
    output_file_.write_line();
    powder_pattern_calculator.set_two_theta_start( powder_pattern_.two_theta_start() );
    powder_pattern_calculator.set_two_theta_end( powder_pattern_.two_theta_end() );
    powder_pattern_calculator.set_two_theta_step( powder_pattern_.average_two_theta_step() );
    powder_pattern_calculator.calculate_reflection_list( true );
    ReflectionList reflection_list = powder_pattern_calculator.reflection_list();
    insert_keyword_and_value( "_cell_measurement_reflns_used", size_t2string( reflection_list.size() ) );
    insert_keyword_and_value( "_cell_measurement_theta_min", double2string( powder_pattern_.two_theta_start().value_in_degrees() ) );
    insert_keyword_and_value( "_cell_measurement_theta_max", double2string( powder_pattern_.two_theta_end().value_in_degrees() ) );
    insert_keyword_and_value( "_diffrn_reflns_number", size_t2string( reflection_list.size() ) );
    insert_keyword_and_value( "_pd_meas_scan_method", "cont" );
    insert_keyword_and_value( "_refine_ls_matrix_type", "fullcycle" );
    // Write gof^2
    iLine = file_inp_.find( "gof" );
    words = split( file_inp_.line( iLine ) );
    if ( words.size() == 2 )
        insert_keyword_and_value( "_refine_ls_goodness_of_fit_all", words[1] );
    // Write r_bragg
    iLine = file_inp_.find( "r_bragg" );
    if ( iLine != std::string::npos )
    {
        words = split( file_inp_.line( iLine ) );
        if ( words.size() == 2 )
            insert_keyword_and_value( "_refine_ls_R_I_factor", words[1] );
    }
    insert_keyword( "_refine_ls_number_parameters" );
    // Count all the distance, angle and planarity restraints
    // Ignore the ones that have been commented out
    size_t nrestraints( 0 );
    iLine = file_inp_.find( "Distance_Restrain" );
    while ( iLine != std::string::npos )
    {
        std::string line_str = strip( file_inp_.line( iLine ) );
        if ( line_str[0] != '\'' )
            ++nrestraints;
        iLine = file_inp_.find( "Distance_Restrain", iLine+1 );
    }
    iLine = file_inp_.find( "Angle_Restrain" );
    while ( iLine != std::string::npos )
    {
        std::string line_str = strip( file_inp_.line( iLine ) );
        if ( line_str[0] != '\'' )
            ++nrestraints;
        iLine = file_inp_.find( "Angle_Restrain", iLine+1 );
    }
    iLine = file_inp_.find( "Flatten" );
    while ( iLine != std::string::npos )
    {
        std::string line_str = strip( file_inp_.line( iLine ) );
        if ( line_str[0] != '\'' )
            ++nrestraints;
        iLine = file_inp_.find( "Flatten", iLine+1 );
    }
    insert_keyword_and_value( "_refine_ls_number_restraints", size_t2string( nrestraints ) );
    insert_keyword_and_value( "_refine_ls_number_constraints", "0" );
    if ( replace_hydrogen_atoms_ )
       insert_keyword_and_value( "_refine_ls_hydrogen_treatment", "noref" );
    else
       insert_keyword_and_value( "_refine_ls_hydrogen_treatment", "refxyz" );
    insert_keyword_and_value( "_refine_ls_weighting_scheme", "sigma" );
    insert_keyword_and_value( "_refine_ls_weighting_details", "'w=1/\\s[Y~obs~]^2^'" );
    insert_keyword_and_value( "_refine_ls_shift/su_max", "0.001" );
    insert_keyword( "_refine_diff_density_min" );
    insert_keyword( "_refine_diff_density_max" );
    insert_keyword_and_value( "_refine_ls_extinction_method", "none" );
//    output_file_.write_line( "loop_" );
//    output_file_.write_line( "    _atom_type_symbol" );
//    output_file_.write_line( "    _atom_type_description" );
//    output_file_.write_line( "    _atom_type_scat_dispersion_real" );
//    output_file_.write_line( "    _atom_type_scat_dispersion_imag" );
//    output_file_.write_line( "    _atom_type_scat_source" );
//    output_file_.write_line( "    S S 0.3331 0.5567 'International Tables for Crystallography, Vol. C'" );
//    output_file_.write_line( "    O O 0.0492 0.0322 'International Tables for Crystallography, Vol. C'" );
//    output_file_.write_line( "    N N 0.0311 0.0180 'International Tables for Crystallography, Vol. C'" );
//    output_file_.write_line( "    C C 0.0181 0.0091 'International Tables for Crystallography, Vol. C'" );
//    output_file_.write_line( "    H H 0.0000 0.0000 'International Tables for Crystallography, Vol. C'" );
    insert_keyword_and_value( "_computing_data_collection", "'<i>WINX^POW^</i> (Stoe & Cie, 2004)'" );
    insert_keyword_and_value( "_computing_cell_refinement", "'<i>TOPAS-Academic</i> (Coelho, 2012)'" );
    insert_keyword_and_value( "_computing_data_reduction", "'<i>DASH</i> 3.3 (David <i>et al.</i>, 2006)'" );
    insert_keyword_and_value( "_computing_structure_solution", "'<i>DASH</i> 3.3 (David <i>et al.</i>, 2006)'" );
    insert_keyword_and_value( "_computing_structure_refinement", "'<i>TOPAS-Academic</i> (Coelho, 2012)'" );
    insert_keyword_and_value( "_computing_molecular_graphics", "'<i>Mercury</i> (Macrae <i>et al.</i>, 2008)'" );
    insert_keyword_and_value( "_pd_proc_2theta_range_min", double2string( powder_pattern_.two_theta_start().value_in_degrees() ) );
    insert_keyword_and_value( "_pd_proc_2theta_range_max", double2string( powder_pattern_.two_theta_end().value_in_degrees() ) );
    insert_keyword_and_value( "_pd_proc_2theta_range_inc", double2string( powder_pattern_.average_two_theta_step().value_in_degrees() ) );
    insert_keyword_and_value( "_pd_calc_method", "'Rietveld Refinement'" );
    insert_keyword_and_value( "_pd_proc_info_data_reduction", "?" );
    insert_keyword_and_value( "_pd_proc_info_excluded_regions", "none" );
    size_t iLine1 = file_inp_.find( "bkg" ) + 1;
    size_t iLine2 = file_inp_.find( "start_X" );
    insert_keyword_and_value( "_pd_proc_ls_background_function", "'Chebyshev function with " + size_t2string( iLine2 - iLine1 ) + " terms'" );
    insert_keyword_and_value( "_pd_proc_ls_profile_function", "'fundamental parameters'" );
//    PO(@, 0.90656`_0.00378, , 0 1 0 )
    iLine = file_inp_.find( "PO(" );
    if ( iLine == std::string::npos )
        insert_keyword_and_value( "_pd_proc_ls_pref_orient_corr", "none" );
    else
    {
        std::string PO_line = strip( file_inp_.line( iLine ) );
        if ( PO_line[0] == '\'')
            insert_keyword_and_value( "_pd_proc_ls_pref_orient_corr", "none" );
        else
        {
            output_file_.write_line( "_pd_proc_ls_pref_orient_corr" );
            output_file_.write_line( ";" );
            output_file_.write_line( "March-Dollase" );
            PO_line = extract_delimited_text( PO_line, "(", ")" );
            words = splitter.split( PO_line );
            output_file_.write_line( "Direction: [ " + strip( words[3] ) + " ]" );
            words[1] = remove( words[1], '`' );
            words[1] = strip( words[1] );
            size_t iPos = words[1].find( "_" );
            if ( iPos == std::string::npos )
                output_file_.write_line( "Refined parameter: " + words[1] );
            else
            {
                DoubleWithESD dwe( string2double( words[1].substr( 0, iPos ) ), string2double( words[1].substr( iPos + 1 ) ) );
                output_file_.write_line( "Refined parameter: " + dwe.crystallographic_style() );
            }
            output_file_.write_line( ";" );
        }
    }
    // Rp
    iLine = file_inp_.find( "r_p_dash" );
    words = split( file_inp_.line( iLine ) );
    if ( words.size() == 2 )
        insert_keyword_and_value( "_pd_proc_ls_prof_R_factor", words[1] + " # background-subtracted" );
    // Rwp
    iLine = file_inp_.find( "r_wp_dash" );
    words = split( file_inp_.line( iLine ) );
    if ( words.size() == 2 )
        insert_keyword_and_value( "_pd_proc_ls_prof_wR_factor", words[1] + " # background-subtracted" );
    // Rexp (always weighted)
    iLine = file_inp_.find( "r_exp_dash" );
    words = split( file_inp_.line( iLine ) );
    if ( words.size() == 2 )
        insert_keyword_and_value( "_pd_proc_ls_prof_wR_expected", words[1] + " # background-subtracted" );
    insert_keyword_and_value( "_pd_proc_ls_special_details", "?" );
    insert_keyword_and_value( "_pd_proc_wavelength", double2string( wavelength_.average_wavelength() ) );
    output_file_.write_line();
    output_file_.write_line( "#==============================================================================" );
    output_file_.write_line();
    output_file_.write_line( "# 7. ATOMIC COORDINATES AND DISPLACEMENT PARAMETERS" );
    output_file_.write_line();

// Atom lines. In .cif from TOPAS:

//loop_
//_atom_site_label
//_atom_site_type_symbol
//_atom_site_symmetry_multiplicity
//_atom_site_fract_x
//_atom_site_fract_y
//_atom_site_fract_z
//_atom_site_occupancy
//_atom_site_B_iso_or_equiv
//_atom_site_U_iso_or_equiv
//
//C1 C   2 0.2611(2) -0.22243(12) 0.1923(12) 1 4.3(3)     0.05426
//C2 C   2 0.4083(2) -0.20534(12) 0.1533(12) 1 4.3(3)     0.05426
//C3 C   2 0.2133(2) -0.14785(13) 0.1479(13) 1 4.3(3)     0.05426
//C4 C   2 0.1882(2) -0.29827(14) 0.2756(14) 1 4.3(3)     0.05426
//C5 C   2 0.4535(2) -0.12043(13) 0.0638(13) 1 4.3(3)     0.05426

    output_file_.write_line( "loop_" );
    output_file_.write_line( "    _atom_site_type_symbol" );
    output_file_.write_line( "    _atom_site_label" );
    output_file_.write_line( "    _atom_site_fract_x" );
    output_file_.write_line( "    _atom_site_fract_y" );
    output_file_.write_line( "    _atom_site_fract_z" );
    output_file_.write_line( "    _atom_site_U_iso_or_equiv" );
    output_file_.write_line( "    _atom_site_occupancy" );
    output_file_.write_line( "    _atom_site_symmetry_multiplicity" );
    double conversion_factor = 1.0 / ( 8.0 * CONSTANT_PI * CONSTANT_PI );
    if ( replace_hydrogen_atoms_ )
        file_Hmi_.read_file( FileName( directory_, base_name_ + "_Hmi", "cif" ) );
    iLine = file_cif_.find( "_atom_site_U_iso_or_equiv" );
    for ( size_t i( iLine + 1 ); i != file_cif_.size(); ++i )
    {
        words = split( file_cif_.line( i ) );
        if ( words.size() == 0 )
            continue;
        DoubleWithESD dwe_1( words[7] );
        DoubleWithESD dwe_2( conversion_factor * dwe_1.value(), conversion_factor * dwe_1.estimated_standard_deviation() );
        if ( replace_hydrogen_atoms_ && element_from_atom_label( words[0] ).is_H_or_D() )
        {
            size_t iLine_Hmi = file_Hmi_.find_whole_word( words[0] );
            // In Mercury, GRACE and Materials Studio, the x,y,z coordinates are columns 2,3,4 (zero-based)
            if ( iLine_Hmi == std::string::npos )
                throw std::runtime_error( "GeneratePowderCIF::generate(): H atom not found: >" + words[0] + "<" );
            std::vector< std::string > words_2 = split( file_Hmi_.line( iLine_Hmi ) );
            output_file_.write_line( "    " + pad( words[1], 2 ) +
                                        " " + pad( relabel( words[0] ), longest_label_size_ ) +
                                        " " + pad_plus( words_2[2], longest_x_ ) +
                                        " " + pad_plus( words_2[3], longest_y_ ) +
                                        " " + pad_plus( words_2[4], longest_z_ ) +
                                        " " + pad( dwe_2.crystallographic_style(), 12 ) +
                                        " " + words[6] +
                                        " " + words[2] );
        }
        else
        {
            output_file_.write_line( "    " + pad( words[1], 2 ) +
                                        " " + pad( relabel( words[0] ), longest_label_size_ ) +
                                        " " + pad_plus( words[3], longest_x_ ) +
                                        " " + pad_plus( words[4], longest_y_ ) +
                                        " " + pad_plus( words[5], longest_z_ ) +
                                        " " + pad( dwe_2.crystallographic_style(), 12 ) +
                                        " " + words[6] +
                                        " " + words[2] );
        }
    }
    if ( replace_hydrogen_atoms_ )
    {
        output_file_.write_line( "_geom_special_details" );
        output_file_.write_line( ";" );
        output_file_.write_line( "The geometry corresponds to the unit cell and positions for the non-H atoms as" );
        output_file_.write_line( "determined from the Rietveld refinement, with CASTEP optimised positions for" );
        output_file_.write_line( "the H atoms." );
        output_file_.write_line( ";" );
    }
    else
        insert_keyword_and_value( "_geom_special_details", "?" );
    write_bond_part();
    write_angle_part();
    insert_keyword_and_value( "_iucr_refine_instructions_details", "?" );
    
// Generate powder cif file
    
//data_I
//loop_
//_pd_meas_2theta_scan
//_pd_meas_intensity_total
//_pd_calc_intensity_total
//_pd_proc_ls_weight
//     3.00000   1685.00000   1663.69555      0.00059
//     3.02000   1675.00000   1659.66414      0.00060
//     3.04000   1634.00000   1655.59337      0.00061
//     3.06000   1696.00000   1651.48754      0.00059
//     3.08000   1653.00000   1647.35073      0.00060
    output_file_xrpd_.write_line( "data_I" );
    output_file_xrpd_.write_line( "loop_" );
    output_file_xrpd_.write_line( "_pd_meas_2theta_scan" );
    output_file_xrpd_.write_line( "_pd_meas_intensity_total" );
    output_file_xrpd_.write_line( "_pd_calc_intensity_total" );
    output_file_xrpd_.write_line( "_pd_proc_ls_weight" );
    for ( size_t i( 0 ); i != file_pro.size(); ++i )
    {
        words = split( file_pro.line( i ) );
        output_file_xrpd_.write_line( words[0] + " " + words[1] + " " + words[2] + " " + double2string( 1.0 / square( string2double( words[3] ) ) ) );
    }

// Generate R input file
    generate_R_input_file( ZOOM_OVER_40 );
}

// ********************************************************************************

void GeneratePowderCIF::generate_R_input_file( const ZoomPolicy zoom_policy )
{
    output_file_R_.write_line( std::string( "setwd( \"" ) + escape_slashes( directory_ ) + "\" ) " );
    output_file_R_.write_line();
    output_file_R_.write_line( "# 2theta, experimental pattern, calculated pattern" );
    output_file_R_.write_line( "profile = read.table( \"" + file_name_pro_ + "\" )" );
    output_file_R_.write_line();
    output_file_R_.write_line( "# Tickmarks");
    output_file_R_.write_line( "tickmarks = read.table( \"" + file_name_tic_ + "\" )" );
    output_file_R_.write_line();
    output_file_R_.write_line( "ndata = length( profile[ , 1 ] )" );
    output_file_R_.write_line();
    output_file_R_.write_line( "tiff( \"refinement_fit.tif\", width = 140, height = 100, units = \"mm\", res = 600, compression = \"lzw\", pointsize = 12 )" );
    output_file_R_.write_line( "par( mar = c( 4, 5, 1, 1 ) )" );
    output_file_R_.write_line();
    output_file_R_.write_line( "y_min = -400" );
    output_file_R_.write_line( "y_max = max( profile[,2], profile[,3] )" );
    output_file_R_.write_line();
    output_file_R_.write_line( "diff_offset = -100" );
    output_file_R_.write_line();
    output_file_R_.write_line( "diff_w = profile[,2]-profile[,3]" );
    output_file_R_.write_line();
    output_file_R_.write_line( "# Make empty plot" );
    output_file_R_.write_line( "plot( 0, 0, col = \"white\", xlim = c( 0.0, " + double2string( powder_pattern_.two_theta_end().value_in_degrees() ) + " ), ylim = c( y_min, y_max )," );
    output_file_R_.write_line( " xlab = expression( 2*italic(theta) / \"\\u00B0\" ), las = 1, xaxs = \"i\", mgp = c( 2, 0.75, 0 ), ylab = \"\" )" );
    output_file_R_.write_line();
    output_file_R_.write_line( "title( ylab = expression( italic(I) / counts ), line = 3 )" );
    output_file_R_.write_line();
    bool include_zoom_at_end( false );
    if ( zoom_policy == ALWAYS_ZOOM )
        include_zoom_at_end = true;
    if ( zoom_policy == ZOOM_OVER_40 )
        include_zoom_at_end = powder_pattern_.two_theta_end().value_in_degrees() >= 40.0;
    if ( include_zoom_at_end )
    {
        output_file_R_.write_line( "# High-angle limit for zoom" );
        output_file_R_.write_line( "ndata_zoom = 1570" );
        output_file_R_.write_line();
        output_file_R_.write_line( "# Experimental data as red crosses" );
        output_file_R_.write_line( "points( profile[1:ndata_zoom,1], profile[1:ndata_zoom,2], pch = \"+\", cex = 0.4, col = \"red\" )" );
        output_file_R_.write_line();
        output_file_R_.write_line( "# Calculated pattern as blue line on top" );
        output_file_R_.write_line( "lines( profile[1:ndata_zoom,1], profile[1:ndata_zoom,3], lwd = 0.75, col = \"blue\" )" );
        output_file_R_.write_line();
        output_file_R_.write_line( "# Difference curve as black line at the bottom" );
        output_file_R_.write_line( "lines( profile[1:ndata_zoom,1], diff_w[1:ndata_zoom] + diff_offset, lwd = 0.75, col = \"black\" )" );
        output_file_R_.write_line();
        output_file_R_.write_line( "# Multiply high-angle region by zoom_factor" );
        output_file_R_.write_line( "zoom_factor = 6" );
        output_file_R_.write_line();
        output_file_R_.write_line( "# The dividing line for the zoomed-in high angle part" );
        output_file_R_.write_line( "lines( c( profile[ndata_zoom,1], profile[ndata_zoom,1] ), c( diff_offset, y_max ), lwd = 1.5, col = \"red\" )" );
        output_file_R_.write_line();
        output_file_R_.write_line( "# Experimental data as red crosses" );
        output_file_R_.write_line( "points( profile[ndata_zoom:ndata,1], zoom_factor*profile[ndata_zoom:ndata,2], pch = \"+\", cex = 0.4, col = \"red\" )" );
        output_file_R_.write_line();
        output_file_R_.write_line( "# Calculated pattern as blue line on top" );
        output_file_R_.write_line( "lines( profile[ndata_zoom:ndata,1], zoom_factor*profile[ndata_zoom:ndata,3], lwd = 0.75, col = \"blue\" )" );
        output_file_R_.write_line();
        output_file_R_.write_line( "# Difference curve as black line at the bottom" );
        output_file_R_.write_line( "lines( profile[ndata_zoom:ndata,1], zoom_factor*diff_w[ndata_zoom:ndata] + diff_offset, lwd = 0.75, col = \"black\" )" );
    }
    else
    {
        output_file_R_.write_line( "# Experimental data as red crosses" );
        output_file_R_.write_line( "points( profile[1:ndata,1], profile[1:ndata,2], pch = \"+\", cex = 0.4, col = \"red\" )" );
        output_file_R_.write_line();
        output_file_R_.write_line( "# Calculated pattern as blue line on top" );
        output_file_R_.write_line( "lines( profile[1:ndata,1], profile[1:ndata,3], lwd = 0.75, col = \"blue\" )" );
        output_file_R_.write_line();
        output_file_R_.write_line( "# Difference curve as black line at the bottom" );
        output_file_R_.write_line( "lines( profile[1:ndata,1], diff_w[1:ndata] + diff_offset, lwd = 0.75, col = \"black\" )" );
    }
    output_file_R_.write_line();
    output_file_R_.write_line( "# Reflection marks" );
    output_file_R_.write_line();
    output_file_R_.write_line( "# Reflection marks drawn from offset to offset+length" );
    output_file_R_.write_line( "mark_offset = -450" );
    output_file_R_.write_line( "mark_length = 100" );
    output_file_R_.write_line();
    output_file_R_.write_line( "ntick = length( tickmarks[,1] )" );
    output_file_R_.write_line( "arrows( tickmarks[,1], rep( mark_offset, ntick ), y1 = rep( mark_offset+mark_length, ntick ), length = 0, col = \"blue\", lwd = 0.75 )" );
    output_file_R_.write_line();
    output_file_R_.write_line( "dev.off()" );
}

// ********************************************************************************

// The new labels of the non-hydrogen atoms are checked. They must be either C1, C2, C3, N1, N2, N3, etc. or C1, C2, N3, C4, N5, etc.
void GeneratePowderCIF::check_new_labels() const
{
    if ( ! relabel_ )
        return;
    std::map< Element, std::set< size_t > > labels;
    std::set< size_t > all_numbers;
    bool all_numbers_unique( true );
    for ( std::map< std::string, std::string >::const_iterator it( labels_.begin() ); it != labels_.end(); ++it )
    {
        std::string atom_label = it->second;
        // Only keep what is a digit
        std::string digits_only;
        for ( size_t iPos( 0 ); iPos != atom_label.length(); ++iPos )
        {
            if ( isdigit( atom_label[ iPos ] ) )
                digits_only += atom_label[iPos];
        }
        size_t number = string2integer( digits_only );
        // Check if number already exists for this element
        Element element = element_from_atom_label( atom_label );
        if ( element.is_H_or_D() )
            continue;
        std::set< size_t >::const_iterator it1 = labels[ element ].find( number );
        if ( it1 == labels[ element ].end() )
            labels[ element ].insert( number );
        else
            throw std::runtime_error( "GeneratePowderCIF::check_new_labels(): number not unique per element: " + size_t2string( number ) );
        if ( ! all_numbers_unique )
            continue;
        std::set< size_t >::const_iterator it2 = all_numbers.find( number );
        if ( it2 == all_numbers.end() )
            all_numbers.insert( number );
        else
            all_numbers_unique = false;
    }
    // First check if label numbering starts at 1 and consecutive per element
    bool consecutive_per_element( true );
    for ( std::map< Element, std::set< size_t > >::const_iterator it1( labels.begin() ); it1 != labels.end(); ++it1 )
    {
        if ( it1->first.is_H_or_D() )
            continue;
        size_t iCurrent( 1 );
        for ( std::set< size_t >::const_iterator it2( it1->second.begin() ); it2 != it1->second.end(); ++it2 )
        {
            if ( *it2 != iCurrent )
            {
                consecutive_per_element = false;
                break;
            }
            ++iCurrent;
        }
        if ( ! consecutive_per_element )
            break;
    }
    if ( consecutive_per_element )
        return;
    // Check if label numbering starts at 1 and consecutive for all elements combined
    if ( ! all_numbers_unique )
        throw std::runtime_error( "GeneratePowderCIF::check_new_labels(): per element check failed and not all numbers are unique." );
    size_t iCurrent( 1 );
    for ( std::set< size_t >::const_iterator it2( all_numbers.begin() ); it2 != all_numbers.end(); ++it2 )
    {
        if ( *it2 != iCurrent )
            throw std::runtime_error( "GeneratePowderCIF::check_new_labels(): per element check failed and not all numbers are consecutive." );
        ++iCurrent;
    }
}

// ********************************************************************************

void GeneratePowderCIF::insert_keyword( const std::string & keyword )
{
    size_t iLine = file_ext_.find( keyword );
    if ( iLine == std::string::npos )
        insert_keyword_and_value( keyword, "?" );
    else
    {
//        std::vector< std::string > words = split( file_ext_.line( iLine ) );
        // The following does not work because split() removes the enclosing quotes so "_cif_keyword 'two words'" becomes "_cif_keyword two words"
//        if ( words.size() == 2 )
//            output_file_.write_line( pad( keyword, padding_length_ ) + " " + words[1] );
//        else
            output_file_.write_line( file_ext_.line( iLine ) );
    }
}

// ********************************************************************************

void GeneratePowderCIF::insert_keyword_and_value( const std::string & keyword, const std::string & value )
{
    output_file_.write_line( pad( keyword, padding_length_ ) + " " + value );
}

// ********************************************************************************

// @@ If Z' < 1, Mercury completes the molecule and adds *all* bonds to the list, so there will be duplicates
void GeneratePowderCIF::write_bond_part()
{
    output_file_.write_line( "loop_" );
    output_file_.write_line( "    _geom_bond_atom_site_label_1" );
    output_file_.write_line( "    _geom_bond_atom_site_label_2" );
    output_file_.write_line( "    _geom_bond_distance" );
    output_file_.write_line( "    _geom_bond_publ_flag" );
    std::string headers = file_bond_lengths_.line( 0 );
    std::vector< std::string > words = split( headers );
//Number	Atom1	Atom2	Type	Polymeric	Cyclicity	Length	SybylType
//1	C1	C2	Unknown	no	cyclic	1.400(3)	un
    size_t atom1_index( 0 ); // Stupid initialisation to silence compiler warnings
    size_t atom2_index( 0 ); // Stupid initialisation to silence compiler warnings
    size_t length_index( 0 ); // Stupid initialisation to silence compiler warnings
    // Very simple algorithm, not robust, find first occurence and that's it.
    for ( size_t i( 0 ); i != words.size(); ++i )
    {
        if ( words[i] == "Atom1" )
            atom1_index = i;
        if ( words[i] == "Atom2" )
            atom2_index = i;
        if ( words[i] == "Length" )
            length_index = i;
    }
    size_t iLine = file_Hmi_.find( "_atom_site_label" );
    for ( size_t i( 1 ); i != file_bond_lengths_.size(); ++i )
    {
        words = split( file_bond_lengths_.line( i ) );
        // If it is an X-H bond, replace it
        if ( replace_hydrogen_atoms_ && ( element_from_atom_label( words[atom1_index] ).is_H_or_D() || element_from_atom_label( words[atom2_index] ).is_H_or_D() ) )
        {
            // In Mercury, GRACE and Materials Studio, the x,y,z coordinates are columns 2,3,4 (zero-based)
            // Calculate the bond length
            Vector3D lhs;
            Vector3D rhs;
            bool atom_1_found( false );
            bool atom_2_found( false );
            for ( size_t j( iLine + 1 ); j != file_Hmi_.size(); ++j )
            {
                std::vector< std::string > words_2 = split( file_Hmi_.line( j ) );
                if ( words_2.size() == 0 )
                    continue;
                if ( ( ! atom_1_found ) && ( words_2[0] == words[atom1_index] ) )
                {
                    lhs = Vector3D( string2double( words_2[2] ), string2double( words_2[3] ), string2double( words_2[4] ) );
                    atom_1_found = true;
                }
                else if ( ( ! atom_2_found ) && ( words_2[0] == words[atom2_index] ) )
                {
                    rhs = Vector3D( string2double( words_2[2] ), string2double( words_2[3] ), string2double( words_2[4] ) );
                    atom_2_found = true;
                }
                if ( atom_1_found && atom_2_found )
                    break;
            }
            if ( ! ( atom_1_found && atom_2_found ) )
                std::cout << "Atom label mismatch. Hint: Mercury may change atom labels when adding hydrogen atoms."<< std::endl;
            if ( ! atom_1_found )
                throw std::runtime_error( "GeneratePowderCIF::write_bond_part(): atom label not found >" + words[atom1_index] + "<" );
            if ( ! atom_2_found )
                throw std::runtime_error( "GeneratePowderCIF::write_bond_part(): atom label not found >" + words[atom2_index] + "<" );
            double bond_length;
            Vector3D dummy;
            crystal_structure_.shortest_distance( lhs, rhs, bond_length, dummy );
            output_file_.write_line( "    " + pad( relabel( words[atom1_index] ), longest_label_size_ ) + " " + pad( relabel( words[atom2_index] ), longest_label_size_ ) + " " + double2string( bond_length, 5 ) + " no" );
        }
        else
            output_file_.write_line( "    " + pad( relabel( words[atom1_index] ), longest_label_size_ ) + " " + pad( relabel( words[atom2_index] ), longest_label_size_ ) + " " + words[length_index] + " no" );
    }
}

// ********************************************************************************

// @@ If Z' < 1, Mercury completes the molecule and adds *all* angles to the list, so there will be duplicates
void GeneratePowderCIF::write_angle_part()
{
    output_file_.write_line( "loop_" );
    output_file_.write_line( "    _geom_angle_atom_site_label_1" );
    output_file_.write_line( "    _geom_angle_atom_site_label_2" );
    output_file_.write_line( "    _geom_angle_atom_site_label_3" );
    output_file_.write_line( "    _geom_angle" );
    output_file_.write_line( "    _geom_angle_publ_flag" );
    std::string headers = file_valence_angles_.line( 0 );
    std::vector< std::string > words = split( headers );
//Number	Atom1	Atom2	Atom3	Angle
//1	C2	C1	C3	108.1(2)
    size_t atom1_index( 0 ); // Stupid initialisation to silence compiler warnings
    size_t atom2_index( 0 ); // Stupid initialisation to silence compiler warnings
    size_t atom3_index( 0 ); // Stupid initialisation to silence compiler warnings
    size_t angle_index( 0 ); // Stupid initialisation to silence compiler warnings
    // Very simple algorithm, not robust, find first occurence and that's it.
    for ( size_t i( 0 ); i != words.size(); ++i )
    {
        if ( words[i] == "Atom1" )
            atom1_index = i;
        if ( words[i] == "Atom2" )
            atom2_index = i;
        if ( words[i] == "Atom3" )
            atom3_index = i;
        if ( words[i] == "Angle" )
            angle_index = i;
    }
    size_t iLine = file_Hmi_.find( "_atom_site_label" );
    for ( size_t i( 1 ); i != file_valence_angles_.size(); ++i )
    {
        words = split( file_valence_angles_.line( i ) );
        // If it is an X-Y-H angle, replace it
        if ( replace_hydrogen_atoms_ && ( element_from_atom_label( words[atom1_index] ).is_H_or_D() || element_from_atom_label( words[atom3_index] ).is_H_or_D() ) )
        {
            // In Mercury, GRACE and Materials Studio, the x,y,z coordinates are columns 2,3,4 (zero-based)
            // Calculate the valence angle
            Vector3D atom_1;
            Vector3D atom_2;
            Vector3D atom_3;
            bool atom_1_found( false );
            bool atom_2_found( false );
            bool atom_3_found( false );
            for ( size_t j( iLine + 1 ); j != file_Hmi_.size(); ++j )
            {
                std::vector< std::string > words_2 = split( file_Hmi_.line( j ) );
                if ( words_2.size() == 0 )
                    continue;
                if ( ( ! atom_1_found ) && ( words_2[0] == words[atom1_index] ) )
                {
                    atom_1 = Vector3D( string2double( words_2[2] ), string2double( words_2[3] ), string2double( words_2[4] ) );
                    atom_1_found = true;
                }
                else if ( ( ! atom_2_found ) && ( words_2[0] == words[atom2_index] ) )
                {
                    atom_2 = Vector3D( string2double( words_2[2] ), string2double( words_2[3] ), string2double( words_2[4] ) );
                    atom_2_found = true;
                }
                else if ( ( ! atom_3_found ) && ( words_2[0] == words[atom3_index] ) )
                {
                    atom_3 = Vector3D( string2double( words_2[2] ), string2double( words_2[3] ), string2double( words_2[4] ) );
                    atom_3_found = true;
                }
                if ( atom_1_found && atom_2_found && atom_3_found )
                    break;
            }
            if ( ! ( atom_1_found && atom_2_found && atom_3_found ) )
                std::cout << "Atom label mismatch. Hint: Mercury may change atom labels when adding hydrogen atoms."<< std::endl;
            if ( ! atom_1_found )
                throw std::runtime_error( "GeneratePowderCIF::write_bond_part(): atom label not found >" + words[atom1_index] + "<" );
            if ( ! atom_2_found )
                throw std::runtime_error( "GeneratePowderCIF::write_bond_part(): atom label not found >" + words[atom2_index] + "<" );
            if ( ! atom_3_found )
                throw std::runtime_error( "GeneratePowderCIF::write_bond_part(): atom label not found >" + words[atom3_index] + "<" );
            double dummy;
            Vector3D difference_vector_1;
            crystal_structure_.shortest_distance( atom_2, atom_1, dummy, difference_vector_1 );
            Vector3D difference_vector_2;
            crystal_structure_.shortest_distance( atom_2, atom_3, dummy, difference_vector_2 );
            // Convert fractional coordinates to Cartesian coordinates
            difference_vector_1 = crystal_lattice_.fractional_to_orthogonal( difference_vector_1 );
            difference_vector_2 = crystal_lattice_.fractional_to_orthogonal( difference_vector_2 );
            Angle valence_angle = angle( difference_vector_1, difference_vector_2 );
            output_file_.write_line( "    " + pad( relabel( words[atom1_index] ), longest_label_size_ ) + " " +
                                              pad( relabel( words[atom2_index] ), longest_label_size_ ) + " " +
                                              pad( relabel( words[atom3_index] ), longest_label_size_ ) + " " + double2string( valence_angle.value_in_degrees() ) + " no" );
        }
        else
            output_file_.write_line( "    " + pad( relabel( words[atom1_index] ), longest_label_size_ ) + " " +
                                              pad( relabel( words[atom2_index] ), longest_label_size_ ) + " " +
                                              pad( relabel( words[atom3_index] ), longest_label_size_ ) + " " + words[angle_index] + " no" );
    }
}

// ********************************************************************************

std::string GeneratePowderCIF::relabel( const std::string & old_label ) const
{
    if ( ! relabel_ )
        return old_label;
    std::map< std::string, std::string >::const_iterator it = labels_.find( old_label );
    if ( it == labels_.end() )
        throw std::runtime_error( "GeneratePowderCIF::relabel(): old label not found: >" + old_label + "<" );
    return it->second;
}

// ********************************************************************************

