#ifndef POWDERPATTERNCALCULATOR_H
#define POWDERPATTERNCALCULATOR_H

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

#include "Angle.h"
#include "FingerCoxJephcoat.h"
#include "MillerIndices.h"
#include "PointGroup.h"
#include "ReflectionList.h"
#include "Wavelength.h"

class CrystalStructure;
class PowderPattern;

#include <set>

// The mixing parameter for the pseudo-Voigt (eta) cannot be set because originally the peak shape was intended to be flexible.
// But pseudo-Voigt works so well and it is required for Finger-Cox-Jephcoat to work, so we
// might just as well consider it as hard-coded.
class PowderPatternCalculator
{
public:

    explicit PowderPatternCalculator( const CrystalStructure & crystal_structure );

    Wavelength wavelength() const { return wavelength_; }
    void set_wavelength( const Wavelength & wavelength ) { wavelength_ = wavelength; }
    Angle two_theta_start() const { return two_theta_start_; }
    void set_two_theta_start( const Angle two_theta_start ) { two_theta_start_ = two_theta_start; }
    Angle two_theta_end() const { return two_theta_end_; }
    void set_two_theta_end( const Angle two_theta_end ) { two_theta_end_ = two_theta_end; }
    Angle two_theta_step() const { return two_theta_step_; }
    void set_two_theta_step( const Angle two_theta_step );
    double FWHM() const { return FWHM_; }
    void set_FWHM( const double FWHM ) { FWHM_ = FWHM; }

    ReflectionList reflection_list() const { return reflection_list_; }

    // Note that the zero-point error is the error itself, not the correction for it.
    // Sample displacement gives rise to a positive error of the order of, say 0.04,
    // which moves the pattern to the right.
    // I think that the +/- convention is the same as in DASH and TOPAS.
    void set_zero_point_error( const Angle zero_point_error );

    void unset_zero_point_error() { include_zero_point_error_ = false; }

    bool include_zero_point_error() const { return include_zero_point_error_; }
    Angle zero_point_error() const { return zero_point_error_; }

    // A March-Dollase model is used.
    void set_preferred_orientation( const MillerIndices & miller_indices, const double r );

    void unset_preferred_orientation() { include_preferred_orientation_ = false; }

    bool include_preferred_orientation() const { return include_preferred_orientation_; }
    MillerIndices preferred_orientation_direction() const { return preferred_orientation_direction_; }
    double r() const { return r_; }

    void set_finger_cox_jephcoat( const FingerCoxJephcoat & finger_cox_jephcoat );

    void unset_finger_cox_jephcoat() { include_finger_cox_jephcoat_ = false; }

    bool include_finger_cox_jephcoat() const { return include_finger_cox_jephcoat_; }
    FingerCoxJephcoat finger_cox_jephcoat() const { return finger_cox_jephcoat_; }

// Same for eta and/or peak shape

    void calculate( PowderPattern & powder_pattern );

    // Calculates d, multiplicity and h,k,l.
    // Only stores one representative reflection if multiplicity > 1 (which is always the case because of Friedel's law).
    // The structure factors are NOT calculated.
    // Usually some reflections before and after the exact 2theta limits are included to avoid strange cut-off effects.
    // This can be switched off by setting exact to true.
    void calculate_reflection_list( const bool exact = false );

    void calculate_structure_factors();

    // Sets all structure factors to 1, to get an artificial powder pattern to compare lattices.
    // There is no need for this to be in the class, the same could be achieved through a combination
    // of other member functions.
// This is the default after calling calculate_reflection_list(), so there is no need for this member function.
//    void set_structure_factors_to_1();

//    void calculate_powder_pattern( PowderPattern & powder_pattern );
    void calculate( const ReflectionList & reflection_list, PowderPattern & powder_pattern );

private:
    Wavelength wavelength_;
    Angle two_theta_start_;
    Angle two_theta_end_;
    Angle two_theta_step_;
    double FWHM_;
    ReflectionList reflection_list_;
    bool include_zero_point_error_;
    Angle zero_point_error_;
    bool include_preferred_orientation_;
    MillerIndices preferred_orientation_direction_;
    double r_;
    bool include_finger_cox_jephcoat_;
    FingerCoxJephcoat finger_cox_jephcoat_;
    const CrystalStructure & crystal_structure_; // Creating a copy would be too expensive given that we have tens of thousands of atoms.
    // But what if the crystal structure goes out of scope and the destructor is called? We need a smart pointer here.
    PointGroup Laue_class_;
    
    bool is_systematic_absence( const MillerIndices miller_indices ) const;
    std::set< MillerIndices > calculate_equivalent_reflections( const MillerIndices miller_indices ) const;
};

#endif // POWDERPATTERNCALCULATOR_H

