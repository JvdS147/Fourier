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

#include "Wavelength.h"
#include "BasicMathsFunctions.h"
#include "StringFunctions.h"

#include <iostream> // For warnings.
#include <stdexcept>

namespace
{

static const size_t num_anode_materials = 5;

static const std::string anode_material[num_anode_materials] = { "Cu", "Cr", "Fe", "Co", "Mo" };
static const double anode_wavelength_1[num_anode_materials] = { 1.54056, 2.28970, 1.93604, 1.78897, 0.70930 };
static const double anode_wavelength_2[num_anode_materials] = { 1.54439, 2.29361, 1.93998, 1.79285, 0.71359 };
static const double anode_average_wavelength[num_anode_materials] = { 1.54184, 2.29100, 1.93735, 1.79026, 0.71073 };

} // namespace

// ********************************************************************************

// Throws if radiation_type is SYNCHROTRON, because that would require a wavelength.
Wavelength::Wavelength( const RadiationSource radiation_source, const bool is_monochromated, const bool use_average ):
radiation_source_(radiation_source),
is_monochromated_(is_monochromated),
use_average_(use_average)
{
    if ( radiation_source_ == SYNCHROTRON )
        throw std::runtime_error( "Wavelength::Wavelength(): error: radiation_type is SYNCHROTRON." );
    if ( is_monochromated_ && use_average_ )
        throw std::runtime_error( "Wavelength::Wavelength(): error: cannot use average if wavelength is monochromated, ." );
}

// ********************************************************************************

Wavelength Wavelength::synchrotron_radiation( const double wavelength )
{
    return Wavelength( SYNCHROTRON, wavelength, true, false );
}

// ********************************************************************************

Wavelength Wavelength::determine_from_wavelength( const double wavelength )
{
    for ( size_t i( 0 ); i != num_anode_materials; ++i )
    {
        if ( nearly_equal( wavelength, anode_wavelength_1[i], 0.0001 ) )
            return Wavelength( static_cast<RadiationSource>(i), 0.0, true, false );
        if ( nearly_equal( wavelength, anode_average_wavelength[i], 0.0001 ) )
            return Wavelength( static_cast<RadiationSource>(i), 0.0, false, true );
    }
    return Wavelength( SYNCHROTRON, wavelength, true, false );
}

// ********************************************************************************

double Wavelength::wavelength() const
{
    if ( ( ! is_monochromated_ ) && ( ! use_average_ ) )
        throw std::runtime_error( "Wavelength::wavelength(): error: a single wavelength is not defined." );
    if ( radiation_source_ == SYNCHROTRON )
        return wavelength_;
    if ( use_average_ )
        return anode_average_wavelength[radiation_source_];
    return anode_wavelength_1[radiation_source_];
}

// ********************************************************************************

double Wavelength::wavelength_1() const
{
    if ( radiation_source_ == SYNCHROTRON )
        throw std::runtime_error( "Wavelength::wavelength_1(): error: synchrotron radiation only has one wavelength." );
    if ( use_average_ )
        std::cout << "Wavelength::wavelength_1(): warning: use_average() set to true." << std::endl;
    return anode_wavelength_1[radiation_source_];
}

// ********************************************************************************

double Wavelength::wavelength_2() const
{
    if ( radiation_source_ == SYNCHROTRON )
        throw std::runtime_error( "Wavelength::wavelength_2(): error: synchrotron radiation only has one wavelength." );
    if ( is_monochromated_ )
        throw std::runtime_error( "Wavelength::wavelength_2(): error: the wavelength is monochromated." );
    if ( use_average_ )
        std::cout << "Wavelength::wavelength_2(): warning: use_average() set to true." << std::endl;
    return anode_wavelength_2[radiation_source_];
}

// ********************************************************************************

double Wavelength::average_wavelength() const
{
    if ( radiation_source_ == SYNCHROTRON )
        throw std::runtime_error( "Wavelength::average_wavelength(): error: synchrotron radiation only has one wavelength." );
    if ( is_monochromated_ )
        throw std::runtime_error( "Wavelength::average_wavelength(): error: the radiation is monochromated." );
    if ( ! use_average_ )
        std::cout << "Wavelength::average_wavelength(): warning: use_average() set to false." << std::endl;
    return anode_average_wavelength[radiation_source_];
}

// ********************************************************************************

// For "_diffrn_radiation_type", e.g. 'Cu K\a~1~' or 'Synchrotron'.
std::string Wavelength::cif_style() const
{
    if ( radiation_source_ == SYNCHROTRON )
        return "'Synchrotron'";
    std::string result( "'" + anode_material[radiation_source_] + " K\\a");
    if ( is_monochromated_ )
        result += "~1~";
    result += "'";
    return result;
}

// ********************************************************************************

std::string radiation_source_to_string( const Wavelength::RadiationSource radiation_source )
{
    switch ( radiation_source )
    {
        case Wavelength::SYNCHROTRON : return "SYNCHROTRON";
        case Wavelength::Cu : return "Cu";
        case Wavelength::Cr : return "Cr";
        case Wavelength::Fe : return "Fe";
        case Wavelength::Co : return "Co";
        case Wavelength::Mo : return "Mo";
        default : throw std::runtime_error( "radiation_source_to_string(): error: unknown radiation_source." );
    }
}

// ********************************************************************************

Wavelength::RadiationSource string_to_radiation_source( std::string input )
{
    input = to_upper( input );
    if ( input == "SYNCHROTRON" )
        return Wavelength::SYNCHROTRON;
    if ( input == "CU" )
        return Wavelength::Cu;
    if ( input == "CR" )
        return Wavelength::Cr;
    if ( input == "FE" )
        return Wavelength::Fe;
    if ( input == "CO" )
        return Wavelength::Co;
    if ( input == "MO" )
        return Wavelength::Mo;
    throw std::runtime_error( "string_to_radiation_source(): error: unknown radiation_source." );
}

// ********************************************************************************

// You can only give one wavelength and no ratio. Will recognise an average wavelength (as used in single-crystal diffraction)
// properly, but that information will be lost.
Wavelength::RadiationSource determine_radiation_source_from_wavelength( const double wavelength )
{
    for ( size_t i( 0 ); i != num_anode_materials; ++i )
    {
        if ( nearly_equal( wavelength, anode_wavelength_1[i], 0.0001 ) ||
             nearly_equal( wavelength, anode_average_wavelength[i], 0.0001 ) )
            return static_cast<Wavelength::RadiationSource>(i);
    }
    return Wavelength::SYNCHROTRON;
}

// ********************************************************************************

bool nearly_equal( const Wavelength & lhs, const Wavelength & rhs )
{
    if ( lhs.radiation_source() != rhs.radiation_source() )
        return false;
    if ( lhs.is_monochromated() != rhs.is_monochromated() )
        return false;
    if ( lhs.use_average() != rhs.use_average() )
        return false;
    if ( lhs.radiation_source() == Wavelength::SYNCHROTRON )
    {
        if ( ! nearly_equal( lhs.wavelength(), rhs.wavelength() ) )
            return false;
    }
    return true;
}

// ********************************************************************************

