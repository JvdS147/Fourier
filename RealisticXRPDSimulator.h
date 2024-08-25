#ifndef REALISTICXRPDSIMULATOR_H
#define REALISTICXRPDSIMULATOR_H

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
#include "PowderPattern.h"
#include "RealisticXRPDSimulatorSettings.h"
#include "Wavelength.h"

class CrystalStructure;

/*

*/
class RealisticXRPDSimulator
{
public:

    explicit RealisticXRPDSimulator( const CrystalStructure & crystal_structure, const RealisticXRPDSimulatorSettings & settings = RealisticXRPDSimulatorSettings() );

    RealisticXRPDSimulatorSettings settings() const { return settings_; }
    void set_settings( const RealisticXRPDSimulatorSettings & settings ) { settings_ = settings; }

    // In a powder pattern that consists of multiple contributions (Bragg scatter, background, NaCl), the relative
    // amounts of the individual contributions is given by Bragg_total_signal_normalisation and background_total_signal_normalisation.
    // When adjusting the powder pattern to the number of counts in the highest peaks, that introduces another scale factor.
    // You need that scale factor e.g. to add NaCl as an imurity to the pattern. 
    double scale_factor() const { return scale_factor_; }

// Same for eta and/or peak shape

    // Currently not const because the powder patterns (noise, background, Bragg) are saved.
    PowderPattern calculate();

    // The Bragg diffraction pattern is the ideal pattern, without noise or background,
    // but including preferrerd orientation and peak asymmetry.
    PowderPattern Bragg_diffraction() const;
    PowderPattern background() const;
    PowderPattern noise() const;
    PowderPattern powder_pattern() const;

private:
    RealisticXRPDSimulatorSettings settings_;
    double scale_factor_;
    const CrystalStructure & crystal_structure_; // Creating a copy would be too expensive given that we have tens of thousands of atoms.
    // But what if the crystal structure goes out of scope and the destructor is called? We need a smart pointer here.
    PowderPattern Bragg_diffraction_;
    PowderPattern background_;
    PowderPattern noise_;
    PowderPattern powder_pattern_;
    bool pattern_has_been_calculated_;
};

#endif // REALISTICXRPDSIMULATOR_H

