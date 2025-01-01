#ifndef WAVELENGTH_H
#define WAVELENGTH_H

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

#include <string>

/*
    This class should be called something like RadiationType instead of Wavelength.
    I guess this could also be "Neutrons" or "Electrons".

    There are four different cases:

    1. Synchrotron radiation. Always one wavelength.
    wavelength() returns the wavelength
    wavelength_1() throws
    wavelength_2() throws
    use_average() returns false
    is_monochromated() returns true
    average_wavelength() throws

    2. Monochromated lab source. Always one wavelength.
    wavelength() returns the wavelength
    wavelength_1() returns the wavelength
    wavelength_2() throws
    use_average() returns false
    is_monochromated() returns true
    average_wavelength() throws

    3. Non-monochromated lab source as used for XRPD, i.e. two separate wavelengths.
    wavelength() throws
    wavelength_1() returns the first wavelength (e.g. CuKalpha_1)
    wavelength_2() returns the second wavelength (e.g. CuKalpha_2)
    use_average() returns false
    is_monochromated() returns false
    average_wavelength() returns ( wavelength_1() + ( 0.5 * wavelength_2() ) ) / 1.5, issues warning.

    4. Non-monochromated lab source as used for SX, i.e. one average wavelength.
    wavelength() returns the average wavelength
    wavelength_1() returns the first wavelength (e.g. CuKalpha_1), issues warning
    wavelength_2() returns the second wavelength (e.g. CuKalpha_2), issues warning
    use_average() returns true
    is_monochromated() returns false
    average_wavelength() returns ( wavelength_1() + ( 0.5 * wavelength_2() ) ) / 1.5.

*/
class Wavelength
{
public:

    enum RadiationSource { SYNCHROTRON, Cu, Cr, Fe, Co, Mo };

    // Throws if radiation_type is SYNCHROTRON, because that would require a wavelength.
    explicit Wavelength( const RadiationSource radiation_source = Cu, const bool is_monochromated = true, const bool use_average = false );

    // Named constructors.
    static Wavelength synchrotron_radiation( const double wavelength );
    static Wavelength determine_from_wavelength( const double wavelength );

    RadiationSource radiation_source() const { return radiation_source_; }

    double wavelength() const;
    double wavelength_1() const;
    double wavelength_2() const;
    double average_wavelength() const;

    bool is_monochromated() const { return is_monochromated_; }
    bool use_average() const { return use_average_; }

    // @@ Functionality like this should really be in a separate file.
    // For "_diffrn_radiation_type", e.g. 'Cu K\a~1~' or 'Synchrotron'.
    std::string cif_style() const;

private:
    RadiationSource radiation_source_;
    double wavelength_;
    bool is_monochromated_;
    bool use_average_;

    // Private constructor for named constructors.
    Wavelength( const RadiationSource radiation_source,
                const double wavelength,
                const bool is_monochromated,
                const bool use_average ): radiation_source_(radiation_source), wavelength_(wavelength), is_monochromated_(is_monochromated), use_average_(use_average) {}
};

// @@ Functionality like this should really be in a separate file.
std::string radiation_source_to_string( const Wavelength::RadiationSource radiation_source );
Wavelength::RadiationSource string_to_radiation_source( std::string input );

// You can only give one wavelength and no ratio. Will recognise an average wavelength (as used in single-crystal diffraction)
// properly, but that information will be lost.
Wavelength::RadiationSource determine_radiation_source_from_wavelength( const double wavelength );

bool nearly_equal( const Wavelength & lhs, const Wavelength & rhs );

#endif // WAVELENGTH_H

