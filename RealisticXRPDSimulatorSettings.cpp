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
#include "TextFileWriter.h"
//#include "FileName.h"

#include <iostream>
#include <stdexcept>

// ********************************************************************************

RealisticXRPDSimulatorSettings::RealisticXRPDSimulatorSettings():
wavelength_(1.54056),
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
wavelength_(1.54056),
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
//    tfw.write_line();
//    Wavelength wavelength_;
//    Angle two_theta_start_;
//    Angle two_theta_end_;
//    Angle two_theta_step_;
//    double FWHM_;
//    bool include_zero_point_error_;
//    Angle zero_point_error_;
//    bool include_preferred_orientation_;
//    MillerIndices preferred_orientation_direction_;
//    double r_;
//    bool include_finger_cox_jephcoat_;
//    double A_; // Finger-Cox-Jephcoat.
//    double B_; // Finger-Cox-Jephcoat.
//    bool include_background_;
//    bool include_noise_;
//    bool include_noise_for_zero_background_;
//    size_t noise_for_zero_background_threshold_;
//    double Bragg_total_signal_normalisation_;
//    double background_total_signal_normalisation_;
//    double highest_peak_;

}

// ********************************************************************************

