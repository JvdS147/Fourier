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

#include "RealisticXRPDSimulatorSettings.h"
#include "BasicMathsFunctions.h"
#include "FileName.h"
#include "StringConversions.h"
#include "StringFunctions.h"
#include "TextFileReader_2.h"
#include "TextFileWriter.h"
#include "Utilities.h"

#include <iostream>
#include <stdexcept>

// ********************************************************************************

RealisticXRPDSimulatorSettings::RealisticXRPDSimulatorSettings():
two_theta_start_(5.0,Angle::DEGREES),
two_theta_end_(35.0,Angle::DEGREES),
two_theta_step_(0.015,Angle::DEGREES),
FWHM_(0.1),
include_zero_point_error_(false),
include_preferred_orientation_(false),
preferred_orientation_direction_( 0, 0, 0 ),
r_(1.0),
include_finger_cox_jephcoat_(false),
A_(0.0001),
B_(0.0001),
include_background_(true),
include_noise_(true),
include_noise_for_zero_background_(true),
noise_for_zero_background_threshold_(20),
Bragg_total_signal_normalisation_(10000.0),
background_total_signal_normalisation_(0.2*10000.0),
highest_peak_(10000.0)
{
}

// ********************************************************************************

RealisticXRPDSimulatorSettings::RealisticXRPDSimulatorSettings( const FileName & file_name ):
two_theta_start_(5.0,Angle::DEGREES),
two_theta_end_(35.0,Angle::DEGREES),
two_theta_step_(0.015,Angle::DEGREES),
FWHM_(0.1),
include_zero_point_error_(false),
include_preferred_orientation_(false),
preferred_orientation_direction_( 0, 0, 0 ),
r_(1.0),
include_finger_cox_jephcoat_(false),
A_(0.0001),
B_(0.0001),
include_background_(true),
include_noise_(true),
include_noise_for_zero_background_(true),
noise_for_zero_background_threshold_(20),
Bragg_total_signal_normalisation_(10000.0),
background_total_signal_normalisation_(0.2*10000.0),
highest_peak_(10000.0)
{
    TextFileReader_2 tfr( file_name );
    tfr.purge_comment_lines( "#" );
    Splitter splitter( ":" );
    splitter.set_merge_delimiters( false );
    size_t iPos;
    std::string input;
    iPos = tfr.find_whole_word( "radiation_source" );
    input = extract_variable_value( tfr.line( iPos ), splitter );
    Wavelength::RadiationSource radiation_source = string_to_radiation_source( input );
    if ( radiation_source == Wavelength::SYNCHROTRON )
    {
        iPos = tfr.find_whole_word( "wavelength" );
        input = extract_variable_value( tfr.line( iPos ), splitter );
        wavelength_ = Wavelength::synchrotron_radiation( string2double( input ) );
    }
    else
    {
        iPos = tfr.find_whole_word( "is_monochromated" );
        bool is_monochromated = string2bool( extract_variable_value( tfr.line( iPos ), splitter ) );
        iPos = tfr.find_whole_word( "use_average" );
        bool use_average = string2bool( extract_variable_value( tfr.line( iPos ), splitter ) );
        wavelength_ = Wavelength( radiation_source, is_monochromated, use_average );
    }
    input = extract_variable_value( tfr.line( tfr.find_whole_word( "two_theta_start" ) ), splitter );
    two_theta_start_ = Angle::from_degrees( string2double( input ) );
    input = extract_variable_value( tfr.line( tfr.find_whole_word( "two_theta_end" ) ), splitter );
    two_theta_end_ = Angle::from_degrees( string2double( input ) );
    input = extract_variable_value( tfr.line( tfr.find_whole_word( "two_theta_step" ) ), splitter );
    two_theta_step_ = Angle::from_degrees( string2double( input ) );
    input = extract_variable_value( tfr.line( tfr.find_whole_word( "FWHM" ) ), splitter );
    FWHM_ = string2double( input );
    input = extract_variable_value( tfr.line( tfr.find_whole_word( "include_zero_point_error" ) ), splitter );
    include_zero_point_error_ = string2bool( input );
    input = extract_variable_value( tfr.line( tfr.find_whole_word( "zero_point_error" ) ), splitter );
    zero_point_error_ = Angle::from_degrees( string2double( input ) );
    input = extract_variable_value( tfr.line( tfr.find_whole_word( "include_preferred_orientation" ) ), splitter );
    include_preferred_orientation_ = string2bool( input );
    input = extract_variable_value( tfr.line( tfr.find_whole_word( "preferred_orientation_direction" ) ), splitter );
    preferred_orientation_direction_ = MillerIndices_from_string( input );
    input = extract_variable_value( tfr.line( tfr.find_whole_word( "r" ) ), splitter );
    r_ = string2double( input );
    input = extract_variable_value( tfr.line( tfr.find_whole_word( "include_finger_cox_jephcoat" ) ), splitter );
    include_finger_cox_jephcoat_ = string2bool( input );
    input = extract_variable_value( tfr.line( tfr.find_whole_word( "A" ) ), splitter );
    A_ = string2double( input );
    input = extract_variable_value( tfr.line( tfr.find_whole_word( "B" ) ), splitter );
    B_ = string2double( input );
    input = extract_variable_value( tfr.line( tfr.find_whole_word( "include_background" ) ), splitter );
    include_background_ = string2bool( input );
    input = extract_variable_value( tfr.line( tfr.find_whole_word( "include_noise" ) ), splitter );
    include_noise_ = string2bool( input );
    input = extract_variable_value( tfr.line( tfr.find_whole_word( "include_noise_for_zero_background" ) ), splitter );
    include_noise_for_zero_background_ = string2bool( input );
    input = extract_variable_value( tfr.line( tfr.find_whole_word( "noise_for_zero_background_threshold" ) ), splitter );
    noise_for_zero_background_threshold_ = string2size_t( input );
    input = extract_variable_value( tfr.line( tfr.find_whole_word( "Bragg_total_signal_normalisation" ) ), splitter );
    Bragg_total_signal_normalisation_ = string2double( input );
    input = extract_variable_value( tfr.line( tfr.find_whole_word( "background_total_signal_normalisation" ) ), splitter );
    background_total_signal_normalisation_ = string2double( input );
    input = extract_variable_value( tfr.line( tfr.find_whole_word( "highest_peak" ) ), splitter );
    highest_peak_ = string2double( input );
}

// ********************************************************************************

void RealisticXRPDSimulatorSettings::set_two_theta_step( const Angle two_theta_step )
{
    two_theta_step_ = two_theta_step;
    if ( two_theta_step_ < Angle::from_degrees( TOLERANCE ) )
         throw std::runtime_error( "RealisticXRPDSimulatorSettings::set_two_theta_step(): Error: value must be positive." );
    // There is an approximation in the calculation of the powder pattern that expects the 2theta step to be small.
    if ( two_theta_step_ > Angle::from_degrees( 0.05 ) )
        std::cout << "RealisticXRPDSimulatorSettings::set_two_theta_step(): Warning: because of an internal approximation, 2theta step is expected to be small." << std::endl;
}

// ********************************************************************************

void RealisticXRPDSimulatorSettings::set_zero_point_error( const Angle zero_point_error )
{
    include_zero_point_error_ = true;
    zero_point_error_ = zero_point_error;
    if ( zero_point_error_.nearly_zero() )
        std::cout << "RealisticXRPDSimulatorSettings::set_zero_point_error(): Warning: zero-point error has been set to 0.0." << std::endl;
}

// ********************************************************************************

void RealisticXRPDSimulatorSettings::unset_zero_point_error()
{
    include_zero_point_error_ = false;
    zero_point_error_ = Angle();
}

// ********************************************************************************

void RealisticXRPDSimulatorSettings::set_preferred_orientation( const MillerIndices & miller_indices, const double r )
{
    include_preferred_orientation_ = true;
    preferred_orientation_direction_ = miller_indices;
    r_ = r;
    if ( preferred_orientation_direction_.is_000() )
        std::cout << "RealisticXRPDSimulatorSettings::set_preferred_orientation(): Warning: Miller indices have been set to (000)." << std::endl;
    if ( nearly_equal( r_, 1.0 ) )
        std::cout << "RealisticXRPDSimulatorSettings::set_preferred_orientation(): Warning: r has been set to 1.0." << std::endl;
//    // Check that the PO direction is commensurate with the space-group symmetry.
//    MillerIndices reflection( 37, -23, 3 );
//    Vector3D PO_vector = reciprocal_lattice_point( preferred_orientation_direction_, crystal_structure_.crystal_lattice() );
//    Vector3D H = reciprocal_lattice_point( reflection, crystal_structure_.crystal_lattice() );
//    double reference_dot_product = absolute( PO_vector * H );
//    for ( size_t i( 0 ); i != Laue_class_.nsymmetry_operators(); ++i )
//    {
//        MillerIndices equivalent_reflection = reflection * Laue_class_.symmetry_operator( i );
//        // Now check that the March-Dollase PO corrections are the same for all of them.
//        Vector3D H = reciprocal_lattice_point( equivalent_reflection, crystal_structure_.crystal_lattice() );
//        double current_dot_product = absolute( PO_vector * H );
//        if ( ! nearly_equal( current_dot_product, reference_dot_product ) )
//        {
//            std::cout << "RealisticXRPDSimulator::set_preferred_orientation(): Warning: PO direction is not commensurate with space-group symmetry." << std::endl;
//            return;
//        }
//    }
}

// ********************************************************************************

void RealisticXRPDSimulatorSettings::set_finger_cox_jephcoat( const double A, const double B )
{
    if ( A < 0.0001 )
        throw std::runtime_error( "RealisticXRPDSimulatorSettings::set_finger_cox_jephcoat(): A < 0.0001." );
    if ( B < 0.0001 )
        throw std::runtime_error( "RealisticXRPDSimulatorSettings::set_finger_cox_jephcoat(): B < 0.0001." );
    include_finger_cox_jephcoat_ = true;
    A_ = A;
    B_ = B;
}

// ********************************************************************************

void RealisticXRPDSimulatorSettings::set_include_noise_for_zero_background( const size_t threshold )
{
    include_noise_for_zero_background_ = true;
    noise_for_zero_background_threshold_ = threshold;
    if ( threshold < 20 )
        std::cout << "RealisticXRPDSimulatorSettings::set_include_noise_for_zero_background(): Warning: threshold < 20." << std::endl;
}

// ********************************************************************************

void RealisticXRPDSimulatorSettings::save( const FileName & file_name ) const
{
    TextFileWriter tfw( file_name );
    tfw.write_line( "#Do not change order." );
    tfw.write_line( "radiation_source : " + radiation_source_to_string( wavelength_.radiation_source() ) );
    if ( wavelength_.radiation_source() == Wavelength::SYNCHROTRON )
        tfw.write_line( "wavelength :" + double2string( wavelength_.wavelength() ) );
    else
        tfw.write_line( "wavelength :" + double2string( wavelength_.wavelength_1() ) );
    tfw.write_line( "is_monochromated : " + bool2string( wavelength_.is_monochromated() ) );
    tfw.write_line( "use_average : " + bool2string( wavelength_.use_average() ) );
    tfw.write_line( "two_theta_start : " + double2string( two_theta_start_.value_in_degrees() ) );
    tfw.write_line( "two_theta_end : " + double2string( two_theta_end_.value_in_degrees() ) );
    tfw.write_line( "two_theta_step : " + double2string( two_theta_step_.value_in_degrees() ) );
    tfw.write_line( "FWHM : " + double2string( FWHM_ ) );
    tfw.write_line( "include_zero_point_error : " + bool2string( include_zero_point_error_ ) );
    tfw.write_line( "zero_point_error : " + double2string( zero_point_error_.value_in_degrees() ) );
    tfw.write_line( "include_preferred_orientation : " + bool2string( include_preferred_orientation_ ) );
    tfw.write_line( "preferred_orientation_direction : " + to_string( preferred_orientation_direction_ ) );
    tfw.write_line( "r : " + double2string( r_ ) );
    tfw.write_line( "include_finger_cox_jephcoat : " + bool2string( include_finger_cox_jephcoat_ ) );
    tfw.write_line( "A : " + double2string( A_ ) );
    tfw.write_line( "B : " + double2string( B_ ) );
    tfw.write_line( "include_background : " + bool2string( include_background_ ) );
    tfw.write_line( "include_noise : " + bool2string( include_noise_ ) );
    tfw.write_line( "include_noise_for_zero_background : " + bool2string( include_noise_for_zero_background_ ) );
    tfw.write_line( "noise_for_zero_background_threshold : " + size_t2string( noise_for_zero_background_threshold_ ) );
    tfw.write_line( "Bragg_total_signal_normalisation : " + double2string( Bragg_total_signal_normalisation_ ) );
    tfw.write_line( "background_total_signal_normalisation : " + double2string( background_total_signal_normalisation_ ) );
    tfw.write_line( "highest_peak : " + double2string( highest_peak_ ) );
}

// ********************************************************************************

