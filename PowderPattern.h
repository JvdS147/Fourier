#ifndef POWDERPATTERN_H
#define POWDERPATTERN_H

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

class FileName;

#include "Angle.h"
#include "Wavelength.h"

#include <vector>

class PowderPattern
{
public:

    PowderPattern();

    // Initialises the 2theta values. Intensities and ESDs are initialised to 0.0.
    PowderPattern( const Angle two_theta_start, const Angle two_theta_end, const Angle two_theta_step );

    explicit PowderPattern( const FileName & file_name );

    void reserve( const size_t nvalues );

    // ESD is initialised to std::max( sqrt( intensity ), intensity / 100.0 ).
    void push_back( const Angle two_theta, const double intensity );
    
    void push_back( const Angle two_theta, const double intensity, const double estimated_standard_deviation );

    size_t size() const { return two_theta_values_.size(); }

    bool empty() const { return two_theta_values_.empty(); }

    // 2theta and intensity are recalculated as averages, the ESDs are recalculated as the square root of the sum of the squares.
    void rebin( const size_t bin_size );

    // Returns the *nearest* 2theta value.
    size_t find_two_theta( const Angle two_theta_value ) const;

    // Multiplies intensities and ESDs by factor.
    void scale( const double factor );

    Angle two_theta( const size_t i ) const;
    double intensity( const size_t i ) const;
    double estimated_standard_deviation( const size_t i ) const;
    void set_two_theta( const size_t i, const Angle value );
    // ESD is NOT updated.
    void set_intensity( const size_t i, const double value );
    void set_estimated_standard_deviation( const size_t i, const double value );
    Wavelength wavelength() const { return wavelength_; }
    void set_wavelength( const Wavelength & wavelength ) { wavelength_ = wavelength; }

    Angle average_two_theta_step() const;

    Angle two_theta_start() const;
    Angle two_theta_end() const;

// I don't know if the following two are a good idea. Perhaps it is better to have a function that
// temporarily reduces (but never expands) the range, something like "set_effective_two_theta_range".
// I would then be possible to reset it to the old values again.

    // Uses average_two_theta_step() to add new points. Intensities and ESDs are initialised to 0.0.
    void set_two_theta_start( const Angle two_theta_start );

    // Uses average_two_theta_step() to add new points. Intensities and ESDs are initialised to 0.0.
    void set_two_theta_end( const Angle two_theta_end );

    // Range cannot be extended.
    void reduce_range_to( const Angle two_theta_start, const Angle two_theta_end );

    // Area under the pattern.
    double cumulative_intensity() const;

    double cumulative_intensity( const Angle two_theta_start, const Angle two_theta_end ) const;

    void read_xye( const FileName & file_name );
    void read_xrdml( const FileName & file_name );
    void read_raw( const FileName & file_name );
    void read_mdi( const FileName & file_name );
    void read_brml( const FileName & file_name );
    void read_txt( const FileName & file_name );
    void read_cif( const FileName & file_name );
    void save_xye( const FileName & file_name, const bool include_wave_length ) const;

    // Writes to std::cout the code that is necessary to generate the PowderPattern object.
    // Useful for writing test-suite code that does not rely on external files.
    void generate_code( const bool include_estimated_standard_deviation ) const;

    PowderPattern & operator+=( const PowderPattern & rhs );
    PowderPattern & operator-=( const PowderPattern & rhs );

    // Normalises the highest peak.
    // Returns the scale factor.
    double normalise_highest_peak( const double highest_peak = 10000.0 );

    // Normalises the total signal = area under the pattern = cumulative_intensity() .
    // Returns the scale factor.
    double normalise_total_signal( const double total_signal = 10000.0 );

    // Simply subtracts the value from each 2theta value, i.e. shifts to the left.
    // Sample displacement shifts the powder pattern to the right.
    // If the error is e.g. 0.04, which shifts the pattern to the *right*, then this function, which *corrects* the error,
    // shifts the pattern to the *left*.
    // I think that the +/- convention is the same as in DASH and TOPAS.
    void correct_zero_point_error( const Angle two_theta_value );

    // ESD is initialised to std::max( sqrt( intensity ), intensity / 100.0 ), or to 4.4 if intensity < 20.
    void recalculate_estimated_standard_deviations();

    // I(fixed slit) = I(variable slit) / sin(theta) is a very good approximation.
    void convert_to_fixed_slit();
    void convert_to_variable_slit();

    void add_constant_background( const double background );

    // It is recommended to call add_constant_background() because otherwise the background points with an average of 0.0 will remain 0.0.
    // For a maximum of about 10,000 counts, adding a background of at least 20 counts gives realistic Estimated Standard Deviations and
    // makes all points of the pattern behave as Gaussian.
    // Note that Poisson noise is only defined for positive integer values, so after adding noise all intensities are integers.
    void add_Poisson_noise();

    // If the number of counts is less than threshold, adds threshold, then calculates the Poisson noise,
    // then subtracts the threshold, then makes the remainder positive.
    // If the maximum is scaled to be 10,000 counts, a good threshold value is 20.
    void add_Poisson_noise_including_zero( const size_t threshold = 20 );

    void make_counts_integer();

    // Should not be necessary. Introduced to manipulate data from a tool that extracted a powder pattern from a bitmap picture.
    void sort_two_theta();

    // Should not be necessary. Introduced to manipulate data from a tool that extracted a powder pattern from a bitmap picture.
    void average_if_two_theta_equal();

private:
    Wavelength wavelength_;
    std::vector< Angle > two_theta_values_;
    std::vector< double > intensities_;
    std::vector< double > estimated_standard_deviations_;
};

// Assumes uniform 2theta step size.
bool same_range( const PowderPattern & lhs, const PowderPattern & rhs );

// This is NOT Rene de Gelder's normalised weighted cross correlation!
// normalised_weighted_cross_correlation( A, B ) =  weighted_cross_correlation( A, B ) / sqrt( weighted_cross_correlation( A, A ) * weighted_cross_correlation( B, B ) )
// This function is only public to enable speeding up the calculation of multiple weighted cross correlation functions (which would all need the same normalisation constants).
// Assumes uniform 2theta step size.
double weighted_cross_correlation( const PowderPattern & lhs, const PowderPattern & rhs, Angle l = Angle( 3.0, Angle::DEGREES ) );

// Because powder patterns are always positive, returns a value between 0.0 and 1.0.
// Assumes uniform 2theta step size.
double normalised_weighted_cross_correlation( const PowderPattern & lhs, const PowderPattern & rhs, Angle l = Angle( 3.0, Angle::DEGREES ) );

// The first pattern is supposed to be the experimental pattern, and its ESDs are used as "the" weights. The ESDs of the second pattern are ignored.
// Since the background cannot be determined from the input, it cannot be subtracted.
double Rwp( const PowderPattern & lhs, const PowderPattern & rhs );

PowderPattern calculate_Brueckner_background( const PowderPattern & powder_pattern,
                                              const size_t niterations,
                                              const size_t window,
                                              const bool apply_smoothing,
                                              const size_t smoothing_window );

// It is recommended to call add_constant_background() because otherwise the background points with an average of 0.0 will remain 0.0.
// For a maximum of about 10,000 counts, adding a background of at least 20 counts gives realistic Estimated Standard Deviations and
// makes all points of the pattern behave as Gaussian.
// Note that Poisson noise is only defined for positive integer values, so the noise consists of integers.
// The noise can be positive or negative. ESDs are set to 0.0.
// @@ I guess this should just return a std::vector< double > but I probably did not do that because we have an
// @@ PowderPattern & operator+=( const PowderPattern & rhs );
// @@ but we do not have an
// @@ PowderPattern & operator+=( const std::vector< double > & rhs );
PowderPattern calculate_Poisson_noise( const PowderPattern & powder_pattern );

// If the number of counts is less than threshold, adds threshold, then calculates the Poisson noise,
// then subtracts the threshold, then makes the remainder positive.
// If the maximum is scaled to be 10,000 counts, a good threshold value is 20.
PowderPattern calculate_Poisson_noise_including_zero( const PowderPattern & powder_pattern, const size_t threshold = 20 );

// Useful for Variable Count Time schemes.
// The 2theta values must currently be the same.
// @@ Currently does not allow the "monitor" to be passed, so only suitable for laboratory data.
// For each powder pattern the "number of seconds counted per 2theta step" (noscp2ts) must be provided to allow the intensities to be put onto the same scale.
PowderPattern add_powder_patterns( const std::vector< PowderPattern > & powder_patterns, const std::vector< double > & noscp2ts );

// Splits a powder pattern over n powder patterns, taking into account Poisson statistics.
// So splitting an intensity of 100 counts over two patterns may assign, say, 54 counts to one, 46 to the other.
// There are two ways to calculate the ESDs: the old ESD divided by n or recalculate the ESD for each pattern.
// In one case, the cumulative ESD is simply the original ESD, in the other case, the cumulative ESD is
// sqrt(n) times the old ESD.
// The proper algorithm depends on whether negative intensities can be present and whether all intensities are integers and the
// intensities ae on a reasonable scale (e.g. NOT all between 0.0 and 1.0). This algorithm assume that all
// intensities in the original pattern are non-negative integers (i.e. counts). All intensities in the output patterns
// are non-negative integers.
std::vector< PowderPattern > split( const PowderPattern & powder_pattern, const size_t n, const bool recalculate_ESDs = true );

#endif // POWDERPATTERN_H

