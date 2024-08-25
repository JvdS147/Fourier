#ifndef REALISTICXRPDSIMULATORSETTINGS_H
#define REALISTICXRPDSIMULATORSETTINGS_H

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

#include "Angle.h"
#include "MillerIndices.h"
#include "Wavelength.h"

/*

*/
class RealisticXRPDSimulatorSettings
{
public:

    // Default constructor
    RealisticXRPDSimulatorSettings();

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

    void set_finger_cox_jephcoat( const double A, const double B );
    void unset_finger_cox_jephcoat() { include_finger_cox_jephcoat_ = false; }
    bool include_finger_cox_jephcoat() const { return include_finger_cox_jephcoat_; }
    double A() const { return A_; }
    double B() const { return B_; }

    bool include_background() const { return include_background_; }
    void set_include_background( const bool include_background ) { include_background_ = include_background; }

    bool include_noise() const { return include_noise_; }
    void set_include_noise( const bool include_noise ) { include_noise_ = include_noise; }

    double Bragg_total_signal_normalisation() const { return Bragg_total_signal_normalisation_; }
    void set_Bragg_total_signal_normalisation( const double Bragg_total_signal_normalisation ) { Bragg_total_signal_normalisation_ = Bragg_total_signal_normalisation; }
    double background_total_signal_normalisation() const { return background_total_signal_normalisation_; }
    void set_background_total_signal_normalisation( const double background_total_signal_normalisation ) { background_total_signal_normalisation_ = background_total_signal_normalisation; }

    double highest_peak() const { return highest_peak_; }
    void set_highest_peak( const double highest_peak ) { highest_peak_ = highest_peak; }

private:
    Wavelength wavelength_;
    Angle two_theta_start_;
    Angle two_theta_end_;
    Angle two_theta_step_;
    double FWHM_;
    bool include_zero_point_error_;
    Angle zero_point_error_;
    bool include_preferred_orientation_;
    MillerIndices preferred_orientation_direction_;
    double r_;
    bool include_finger_cox_jephcoat_;
    double A_; // Finger-Cox-Jephcoat.
    double B_; // Finger-Cox-Jephcoat.
    bool include_background_;
    bool include_noise_;
    double Bragg_total_signal_normalisation_;
    double background_total_signal_normalisation_;
    double highest_peak_;
};

#endif // REALISTICXRPDSIMULATORSETTINGS_H

