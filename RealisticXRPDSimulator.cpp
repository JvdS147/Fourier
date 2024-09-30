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

#include "RealisticXRPDSimulator.h"
#include "CrystalStructure.h"
#include "InpWriter.h"
#include "FileName.h"
#include "PowderPatternCalculator.h"
#include "StringFunctions.h"
#include "TextFileWriter.h"
#include "Utilities.h"

#include <iostream>
#include <stdexcept>

// ********************************************************************************

RealisticXRPDSimulator::RealisticXRPDSimulator( const CrystalStructure & crystal_structure, const RealisticXRPDSimulatorSettings & settings ):
settings_(settings),
scale_factor_(1.0),
crystal_structure_(crystal_structure),
pattern_has_been_calculated_(false)
{
    if ( ! crystal_structure.space_group_symmetry_has_been_applied() )
        std::cout << "RealisticXRPDSimulator::RealisticXRPDSimulator( CrystalStructure ): Warning: space-group symmetry has not been applied for input crystal structure."<< std::endl;
}

// ********************************************************************************

PowderPattern RealisticXRPDSimulator::calculate()
{
        std::cout << "Now calculating powder pattern... " << std::endl;
        PowderPatternCalculator powder_pattern_calculator( crystal_structure_ );
        powder_pattern_calculator.set_wavelength( settings_.wavelength() );
        powder_pattern_calculator.set_two_theta_start( settings_.two_theta_start() );
        powder_pattern_calculator.set_two_theta_end( settings_.two_theta_end() );
        powder_pattern_calculator.set_two_theta_step( settings_.two_theta_step() );
        powder_pattern_calculator.set_FWHM( settings_.FWHM() );
        if ( settings_.include_zero_point_error() )
            powder_pattern_calculator.set_zero_point_error( settings_.zero_point_error() );
        if ( settings_.include_preferred_orientation() )
            powder_pattern_calculator.set_preferred_orientation( settings_.preferred_orientation_direction(), settings_.r() );
        if ( settings_.include_finger_cox_jephcoat() )
            powder_pattern_calculator.set_finger_cox_jephcoat( FingerCoxJephcoat( settings_.A(), settings_.B() ) );
        powder_pattern_calculator.calculate( Bragg_diffraction_ );
        Bragg_diffraction_.normalise_total_signal( settings_.Bragg_total_signal_normalisation() );
        powder_pattern_ = Bragg_diffraction_;
        if ( settings_.include_background() )
        {
            PowderPatternCalculator background_powder_pattern_calculator( crystal_structure_ );
            background_powder_pattern_calculator.set_wavelength( settings_.wavelength() );
            background_powder_pattern_calculator.set_two_theta_start( settings_.two_theta_start() );
            background_powder_pattern_calculator.set_two_theta_end( settings_.two_theta_end() );
            background_powder_pattern_calculator.set_two_theta_step( settings_.two_theta_step() );
            background_powder_pattern_calculator.set_FWHM( 5.0 );
            if ( settings_.include_zero_point_error() )
                background_powder_pattern_calculator.set_zero_point_error( settings_.zero_point_error() );
            // We never include PO for the amorphous background.
            if ( settings_.include_finger_cox_jephcoat() )
                background_powder_pattern_calculator.set_finger_cox_jephcoat( FingerCoxJephcoat( settings_.A(), settings_.B() ) );
            background_powder_pattern_calculator.calculate( background_ );
            background_.normalise_total_signal( settings_.background_total_signal_normalisation() );
            powder_pattern_ += background_;
        }
        scale_factor_ = powder_pattern_.normalise_highest_peak( settings_.highest_peak() );
        Bragg_diffraction_.scale( scale_factor_ );
        Bragg_diffraction_.make_counts_integer();
        Bragg_diffraction_.recalculate_estimated_standard_deviations();
        powder_pattern_ = Bragg_diffraction_;
        if ( settings_.include_background() )
        {
            background_.scale( scale_factor_ );
            background_.add_constant_background( 20.0 );
            background_.make_counts_integer();
            background_.recalculate_estimated_standard_deviations();
            powder_pattern_ += background_;
        }
        if ( settings_.include_noise() )
        {
            if ( settings_.include_noise_for_zero_background() )
                noise_ = calculate_Poisson_noise_including_zero( powder_pattern_, settings_.noise_for_zero_background_threshold() );
            else
                noise_ = calculate_Poisson_noise( powder_pattern_ );
            powder_pattern_ += noise_;
        }
        // They are *estimated* standard deviations, so they should be calculated *after* the noise has been introduced.
        powder_pattern_.recalculate_estimated_standard_deviations();
        pattern_has_been_calculated_ = true;
        return powder_pattern_;
}

// ********************************************************************************

PowderPattern RealisticXRPDSimulator::Bragg_diffraction() const
{
    if ( ! pattern_has_been_calculated_ )
        throw std::runtime_error( "RealisticXRPDSimulator::Bragg_diffraction(): Error: the pattern has not been calculated yet." );
    return Bragg_diffraction_;
}

// ********************************************************************************

PowderPattern RealisticXRPDSimulator::background() const
{
    if ( ! pattern_has_been_calculated_ )
        throw std::runtime_error( "RealisticXRPDSimulator::background(): Error: the pattern has not been calculated yet." );
    return background_;
}

// ********************************************************************************

PowderPattern RealisticXRPDSimulator::noise() const
{
    if ( ! pattern_has_been_calculated_ )
        throw std::runtime_error( "RealisticXRPDSimulator::noise(): Error: the pattern has not been calculated yet." );
    return noise_;
}

// ********************************************************************************

PowderPattern RealisticXRPDSimulator::powder_pattern() const
{
    if ( ! pattern_has_been_calculated_ )
        throw std::runtime_error( "RealisticXRPDSimulator::powder_pattern(): Error: the pattern has not been calculated yet." );
    return powder_pattern_;
}

// ********************************************************************************

void RealisticXRPDSimulator::write_TOPAS_inp_file( const FileName & file_name ) const
{
    const std::string aal;
    TextFileWriter text_file_writer( replace_extension( file_name, "inp" ) );
    std::cout << "No restraints will be written out." << std::endl;
    write_preamble( text_file_writer );
    text_file_writer.write_line( "xdd " + FileName( file_name.directory(), file_name.name(), "xye" ).full_name() + " xye_format" );
    text_file_writer.write_line( "  bkg @" );
    for ( size_t i( 0 ); i != 20; ++i )
        text_file_writer.write_line( "    0.0" );
    text_file_writer.write_line( "  start_X       " + double2string( settings_.two_theta_start().value_in_degrees() ) );
    text_file_writer.write_line( "  finish_X      " + double2string( settings_.two_theta_end().value_in_degrees() ) );
    text_file_writer.write_line( "  x_calculation_step " + double2string( settings_.two_theta_step().value_in_degrees() ) );
    text_file_writer.write_line( "'  Specimen_Displacement(@ , 0.0 )" );
    text_file_writer.write_line( "'  Absorption(@ , 0,0 )" );
    text_file_writer.write_line( "  Zero_Error(@ , " + double2string( settings_.zero_point_error().value_in_degrees() ) + " )" );
    text_file_writer.write_line( "'Synchrotron use: LP_Factor( 90.0 )" );
    text_file_writer.write_line( "'Neutrons use: LP_Factor( 90.0 )" );
    text_file_writer.write_line( "'No monochromator use: LP_Factor( 0.0 )" );
    text_file_writer.write_line( "'Ge Monochromator, Cu radiation, use LP_Factor( 27.3 )" );
    text_file_writer.write_line( "'Graphite Monochromator, Cu radiation, use LP_Factor( 26.4 )" );
    text_file_writer.write_line( "'Quartz Monochromator, Cu radiation, use LP_Factor( 26.6 )" );
    text_file_writer.write_line( "  LP_Factor( 26.5 )" );
    text_file_writer.write_line( "'  Variable_Divergence(@ , 30.0 )" );
    if ( settings_.include_finger_cox_jephcoat() )
    {
        text_file_writer.write_line( "  axial_conv" );
        text_file_writer.write_line( "    filament_length @ 4.97890" );
        text_file_writer.write_line( "    sample_length @ 2.43658" );
        text_file_writer.write_line( "    receiving_slit_length @ 5.25714" );
        text_file_writer.write_line( "    axial_n_beta 50" );
    }
    text_file_writer.write_line( "  lam" );
    text_file_writer.write_line( "    ymin_on_ymax 0.001" );
    text_file_writer.write_line( "    la 1 lo " + double2string( settings_.wavelength().wavelength_1() ) );
    text_file_writer.write_line( "  str" );
    text_file_writer.write_line( "    r_bragg  0.0" );
    text_file_writer.write_line( "    CS_G(@ , 107.03272`)" );
    text_file_writer.write_line( "    CS_L(@ , 9999.99881`)" );
    text_file_writer.write_line( "    Strain_G(@ , 0.49554`)" );
    text_file_writer.write_line( "    Strain_L(@ , 0.03347`)" );
//    text_file_writer.write_line( "    prm sh_scale_l" + aal + " 1.0" );
//    text_file_writer.write_line( "    spherical_harmonics_hkl sh_l" + aal );
//    text_file_writer.write_line( "      sh_order 6" );
//    text_file_writer.write_line( "    lor_fwhm = Abs( sh_scale_l" + aal + " * sh_l" + aal + " );" );
//    text_file_writer.write_line( "    prm sh_scale_g" + aal + " 1.0" );
//    text_file_writer.write_line( "    spherical_harmonics_hkl sh_g" + aal );
//    text_file_writer.write_line( "      sh_order 6" );
//    text_file_writer.write_line( "    gauss_fwhm = Abs( sh_scale_g" + aal + " * sh_g" + aal + " );" );
    write_unit_cell( text_file_writer, crystal_structure_.crystal_lattice() );
    text_file_writer.write_line( "    MVW( 0.0, 0.0, 0.0 )");
    text_file_writer.write_line( "    space_group \"" + remove( remove( crystal_structure_.space_group().name(), '_' ), ' ') + "\"");
    text_file_writer.write_line( "    scale @  0.0001" );
    if ( settings_.include_preferred_orientation() )
        text_file_writer.write_line( "    PO(@ , " + double2string( settings_.r() ) + ", , " + settings_.preferred_orientation_direction().to_string() + " )" );
    else
        text_file_writer.write_line( "'    PO(@ , 1.0, , 1 0 0 )" );
    text_file_writer.write_line( "'    PO_Spherical_Harmonics( sh, 6 )" );
    text_file_writer.write_line( "    macro ref_flag" + aal + " { @ }" );
    text_file_writer.write_line( "    prm bnonh" + aal + " 3.0" );
    text_file_writer.write_line( "    prm bh" + aal + " = 1.2 * bnonh" + aal + ";" );
    for ( size_t i( 0 ); i != crystal_structure_.natoms(); ++i )
    {
        text_file_writer.write( "    site " + crystal_structure_.atom( i ).label() + aal + " x ref_flag" + aal + " " + double2string_pad_plus( crystal_structure_.atom( i ).position().x(), 5, ' ' ) +
                                                                                           " y ref_flag" + aal + " " + double2string_pad_plus( crystal_structure_.atom( i ).position().y(), 5, ' ' ) +
                                                                                           " z ref_flag" + aal + " " + double2string_pad_plus( crystal_structure_.atom( i ).position().z(), 5, ' ' ) +
                                     " occ " + pad( crystal_structure_.atom( i ).element().symbol(), 2, ' ' ) + " 1 beq = " );
        if ( crystal_structure_.atom( i ).element().is_H_or_D() )
            text_file_writer.write_line( "bh" + aal + ";" );
        else
            text_file_writer.write_line( "bnonh" + aal + ";" );
    }
    text_file_writer.write_line( "    Out_CIF_STR( " + FileName( file_name.directory(), file_name.name() + "_RR", "cif" ).full_name() + " )" );
    text_file_writer.write_line( "  xdd_out " + FileName( file_name.directory(), file_name.name() + "_profile", "txt" ).full_name() + " load out_record out_fmt out_eqn" );
    text_file_writer.write_line( "  {" );
    text_file_writer.write_line( "      \" %11.5f \" = X;" );
    text_file_writer.write_line( "      \" %11.5f \" = Yobs;" );
    text_file_writer.write_line( "      \" %11.5f \" = Ycalc;" );
    text_file_writer.write_line( "      \" %11.5f\\n\" = SigmaYobs;" );
    text_file_writer.write_line( "  }" );
    text_file_writer.write_line( "  phase_out " + FileName( file_name.directory(), file_name.name() + "_tickmarks", "txt" ).full_name() + " load out_record out_fmt out_eqn" );
    text_file_writer.write_line( "  {" );
    text_file_writer.write_line( "      \" %11.5f -200\\n\" = 2.0 * Rad * Th;" );
    text_file_writer.write_line( "  }" );
}

// ********************************************************************************

