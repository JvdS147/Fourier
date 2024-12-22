#ifndef WAVELENGTH_H
#define WAVELENGTH_H

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

#include <string>

/*
    There is a minor grey zone here: the radiation that is used in single crystal is generally given as the average of the two CuK(alpha)
    wavelengths, so there is only one wavelength, but it is not monochromated. In this class, that counts as "monochromated" in the sense that only
    a single wavelength is considered.
*/
class Wavelength
{
public:

    // Default constructor: CuKa1.
    Wavelength();

    explicit Wavelength( const double wavelength );

    // "CuKa", "CuKa1". "Cu K\a~1~" is also allowed,
//    Wavelength( const std::string & anode_material, const bool is_monochromated );

    // "ratio" is the amount of wavelength_2 and must be given such that the amount of wavelength_1 is 1.0.
    Wavelength( const double wavelength_1, const double wavelength_2, const double ratio );

    double wavelength_1() const { return wavelength_1_; }
    double wavelength_2() const;

    double ratio() const;
    double average_wavelength() const;

    bool is_monochromated() const { return is_monochromated_; }
    void set_wavelength_2( const double wavelength_2, const double ratio );
    void unset_wavelength_2();

    bool is_lab_source() const { return is_lab_source_; }
    void set_is_lab_source( const bool is_lab_source ) { is_lab_source_ = is_lab_source; }

    // For "_diffrn_radiation_type", e.g. 'Cu K\a~1~' or 'Synchrotron'
    std::string cif_style() const;

private:
    double wavelength_1_;
    double wavelength_2_;
    double ratio_;
    bool is_monochromated_;
    bool is_lab_source_;

};

bool nearly_equal( const Wavelength & lhs, const Wavelength & rhs );

#endif // WAVELENGTH_H

