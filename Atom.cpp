/* *********************************************
Copyright (c) 2013-2020, Cornelis Jan (Jacco) van de Streek
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

#include "Atom.h"

#include <stdexcept>

// ********************************************************************************

Atom::Atom( const Element & element,
            const Vector3D & position,
            const std::string & label ):
element_(element),
position_(position),
label_(label),
charge_(0.0),
ADPs_type_(NONE),
Uiso_(0.0),
occupancy_(1.0)
{
}

// ********************************************************************************

Atom::Atom( const Element & element,
            const Vector3D & position,
            const std::string & label,
            const AnisotropicDisplacementParameters & anisotropic_displacement_parameters ):
element_(element),
position_(position),
label_(label),
charge_(0.0),
ADPs_type_(ANISOTROPIC),
anisotropic_displacement_parameters_(anisotropic_displacement_parameters),
occupancy_(1.0)
{
    Uiso_ = anisotropic_displacement_parameters_.U_iso();
}

// ********************************************************************************

Atom::Atom( const Element & element,
            const Vector3D & position,
            const std::string & label,
            const double charge,
            const AnisotropicDisplacementParameters & anisotropic_displacement_parameters ):
element_(element),
position_(position),
label_(label),
charge_(charge),
ADPs_type_(ANISOTROPIC),
anisotropic_displacement_parameters_(anisotropic_displacement_parameters),
occupancy_(1.0)
{
    Uiso_ = anisotropic_displacement_parameters_.U_iso();
}

// ********************************************************************************

AnisotropicDisplacementParameters Atom::anisotropic_displacement_parameters() const
{
    if ( ADPs_type_ == NONE )
        return AnisotropicDisplacementParameters( 0.0 );
    if ( ADPs_type_ == ISOTROPIC )
        return AnisotropicDisplacementParameters( Uiso_ );
    if ( ADPs_type_ == ANISOTROPIC )
        return anisotropic_displacement_parameters_;
    throw std::runtime_error( "Error." );

}

// ********************************************************************************

void Atom::set_anisotropic_displacement_parameters( const AnisotropicDisplacementParameters & anisotropic_displacement_parameters )
{
    anisotropic_displacement_parameters_ = anisotropic_displacement_parameters;
    Uiso_ = anisotropic_displacement_parameters_.U_iso();
    ADPs_type_ = ANISOTROPIC;
}

// ********************************************************************************

double Atom::Uiso() const
{
    if ( ADPs_type_ == NONE )
        return 0.0;
    if ( ADPs_type_ == ISOTROPIC )
        return Uiso_;
    if ( ADPs_type_ == ANISOTROPIC )
        return anisotropic_displacement_parameters_.U_iso();
    throw std::runtime_error( "Error." );

}

// ********************************************************************************

void Atom::set_Uiso( const double Uiso )
{
    Uiso_ = Uiso;
    anisotropic_displacement_parameters_ = AnisotropicDisplacementParameters( Uiso_ );
    ADPs_type_ = ISOTROPIC;
}

// ********************************************************************************
