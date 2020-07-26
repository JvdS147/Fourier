#ifndef MOLECULEINCRYSTAL_H
#define MOLECULEINCRYSTAL_H

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

#include <vector>

/*
    A molecule in a crystal, in other words, the atoms have 3D coordinates and these 3D coordinates are fractional coordinates.
    It has a method "is_on_special_position()", which is also pretty specific to a molecule in a crystal.
    The crystal lattice is not stored with the molecule but must be kept track of by the user.
    
    I guess this means that this molecule can be disordered (whereas a Molecule2D or a Molecule3D cannot be disordered.
    A Molecule2D would be a chemical diagram, just a topology/connectivity, a Molecule3D would be the
    class that you would use in a conformer search, coordinates would be Cartesian and the molecule cannot be disordered.)
*/
class MoleculeInCrystal
{
public:

    // Default constructor
    MoleculeInCrystal();

    size_t natoms() const { return atoms_.size(); }

    void add_atom( const Atom & atom ) { atoms_.push_back( atom ); }
    void add_atoms( const std::vector< Atom > & atoms );

    // This should probably be a member function of CrystalStructure: molecule_is_on_special_position( const size_t i );
  //  bool is_on_special_position( const CrystalStructure & crystal_structure ) const;

    Atom atom( const size_t i ) const { return atoms_[i]; }
    
    void set_atom( const size_t i, const Atom & new_atom ) { atoms_[i] = new_atom; }

private:
    // The atom perception algorithm moves atoms so that they form connected molecules,
    // so we cannot use indices but must really store entire atoms.
    std::vector< Atom > atoms_;
};

#endif // MOLECULEINCRYSTAL_H

