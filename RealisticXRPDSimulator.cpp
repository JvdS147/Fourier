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
#include "PowderPatternCalculator.h"

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

