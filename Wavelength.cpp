/* *********************************************
Copyright (c) 2013-2021, Cornelis Jan (Jacco) van de Streek
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
#include "Utilities.h"

namespace
{

static const size_t num_anode_materials = 5;

static const std::string anode_material[num_anode_materials] = { "Cu", "Cr", "Fe", "Co", "Mo" };
static const double anode_wavelength_1[num_anode_materials] = { 1.54056, 2.28970, 1.93604, 1.78897, 0.70930 };
static const double anode_wavelength_2[num_anode_materials] = { 1.54439, 2.29361, 1.93998, 1.79285, 0.71359 };

} // namespace

// ********************************************************************************

// Default constructor: CuKa1
Wavelength::Wavelength() :
wavelength_1_( 1.54056 ),
wavelength_2_( 0.0 ),
monochromated_( true ),
is_lab_source_( true )
{
}

// ********************************************************************************

Wavelength::Wavelength( const double wavelength ) :
wavelength_1_( wavelength ),
wavelength_2_( 0.0 ),
monochromated_( true )
{
    for ( size_t i( 0 ); i != num_anode_materials; ++i )
    {
        if ( nearly_equal( wavelength, (anode_wavelength_1[i]+anode_wavelength_2[i])/2.0, 0.005 ) )
        {
            is_lab_source_ = true;
            return;
        }
    }
    is_lab_source_ = false;
}

// ********************************************************************************

Wavelength::Wavelength( const double wavelength_1, const double wavelength_2 ) :
wavelength_1_( wavelength_1 ),
wavelength_2_( wavelength_2 ),
monochromated_( false ),
is_lab_source_( true )
{
}

// ********************************************************************************

double Wavelength::average_wavelength() const
{
    if ( monochromated_ )
        return wavelength_1_;
    else
        return ( 2.0 * wavelength_1_ + wavelength_2_ ) / 3.0;
}

// ********************************************************************************

// For "_diffrn_radiation_type", e.g. 'Cu K\a~1~' or "Synchrotron"
std::string Wavelength::cif_style() const
{
    if ( ! is_lab_source_ )
        return "'Synchrotron'";
    for ( size_t i( 0 ); i != num_anode_materials; ++i )
    {
        if ( nearly_equal( wavelength_1_, (anode_wavelength_1[i]+anode_wavelength_2[i])/2.0, 0.005 ) )
        {
            std::string result( "'" + anode_material[i] + " K\\a");
            if ( monochromated_ )
                result += "~1~";
            result += "'";
            return result;
        }
    }
    return "Nonsense";
    // @@ When we are here, we are in trouble.
}

// ********************************************************************************

