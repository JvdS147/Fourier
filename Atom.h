#ifndef ATOM_H
#define ATOM_H

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

#include "AnisotropicDisplacementParameters.h"
#include "Element.h"
#include "Vector3D.h"

class Atom
{
public:

    enum ADPsType { ANISOTROPIC, ISOTROPIC, NONE };

    Atom();

    Atom( const Element & element,
          const Vector3D & position,
          const std::string & label );

    Atom( const Element & element,
          const Vector3D & position,
          const std::string & label,
          const AnisotropicDisplacementParameters & anisotropic_displacement_parameters );

    Atom( const Element & element,
          const Vector3D & position,
          const std::string & label,
          const double charge,
          const AnisotropicDisplacementParameters & anisotropic_displacement_parameters );

    Element element() const { return element_; }

    void set_element( const Element element ) { element_ = element; }

    // Is this position in Cartesian coordinates or in fractional coordinates?
    Vector3D position() const { return position_; }

    // Is this position in Cartesian coordinates or in fractional coordinates?
    void set_position( const Vector3D & position ) { position_ = position; }
    
    std::string label() const { return label_; }

    void set_label( const std::string & label ) { label_ = label; }

    double charge() const { return charge_; }

    void set_charge( const double charge ) { charge_ = charge; }

    double occupancy() const { return occupancy_; }

    void set_occupancy( const double occupancy ) { occupancy_ = occupancy; }

    ADPsType ADPs_type() const { return ADPs_type_; }

    // Set ADP type to NONE.
    void reset_ADPs_type();

    // If ADPs_type_ is NONE or ISOTROPIC, returns diagonal matrix with 0.0 or with Uiso_
    AnisotropicDisplacementParameters anisotropic_displacement_parameters() const;

    void set_anisotropic_displacement_parameters( const AnisotropicDisplacementParameters & anisotropic_displacement_parameters );
    
    double Uiso() const;

    void set_Uiso( const double Uiso );

    // | element | element | number of bonded atoms | number of bonded hydrogen/deuterium atoms | first bonded atom by atomic number, H == D, 0 if none | second bonded atom by atomic number
    // | third bonded atom by atomic number | fourth bonded atom by atomic number | member of three-membered ring | member of four-membered ring
    // | member of five-membered ring | member of six-membered ring | member of seven-membered ring | cyclic
    // All these properties must be independent of the presence of 3D coordinates, they are topological attributes.
    // They cannot be calculated in this class, they can be calculated in e.g. the CrystalStructure class
    std::string topological_attributes() const { return topological_attributes_; }
    
    void set_topological_attributes( const std::string & topological_attributes ) { topological_attributes_ = topological_attributes; }

private:
    Element element_;
    Vector3D position_; // Fractional coordinates
    std::string label_;
    double charge_;
    ADPsType ADPs_type_;
    double Uiso_;
    AnisotropicDisplacementParameters anisotropic_displacement_parameters_;
    double occupancy_;
    std::string disorder_assembly_;
    std::string disorder_group_;
    std::string topological_attributes_;
};

// Throws if element not the same, also averages ADPs, which means that the atoms must have been
// transformed to coincide (e.g. if they are symmetry-related).
// In general, two atoms that need to be averaged were probably related by symmetry
// and the symmetry operator is a fundamental part of the averaging process,
// e.g. because the ADPs must adhere to the site symmetry.
// Currently, this information is spread over several functions,
// but ideally this function takes a SymmetryOperator. But... the translation part would be canonicalised...
// ... and for the ADPs we need the full point group, not just the one symmetry operator.
Atom average( const Atom & lhs, const Atom & rhs );

#endif // ATOM_H

