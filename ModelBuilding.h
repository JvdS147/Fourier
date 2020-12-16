#ifndef MODELBUILDING_H
#define MODELBUILDING_H

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

class CrystalStructure;
class Element;
class Vector3D;

#include "Angle.h"

#include <vector>

// The methyl group is attached to atom_1, away from atom_2
// angle allows control over the rotation of the methyl group as it is attached.
// Returns the coordinates of the three hydrogen atoms
// Everything in Cartesian coordinates
std::vector< Vector3D > add_methyl_group( const Vector3D & atom_1, const Vector3D & atom_2, const Angle angle = Angle() );

// The hydrogen atoms are attached to atom_1, away from atom_2
// angle allows control over the rotation of the hydrogen atoms as they are attached.
// Returns the coordinates of the nhydrogens hydrogen atoms
// Everything in Cartesian coordinates
std::vector< Vector3D > add_hydrogen_atoms( const Vector3D & atom_1, const Vector3D & atom_2, const size_t nhydrogens = 6, const Angle angle = Angle() );

// Add two hydrogen atoms to an sp3 nitrogen atom, with atom_H pointing at the lone pair of the nitrogen atom
// (the two hydrogen atoms are added "pointing away" from atom_H)
// It is as if you add three hydrogen atoms to a methyl group, with the direction of one hydrogen atom supplied by the user.
// Returns the coordinates of the two hydrogen atoms
// Everything in Cartesian coordinates
std::vector< Vector3D > add_2_hydrogen_atoms_to_sp3_nitrogen( const Vector3D & atom_N, const Vector3D & atom_2, const Vector3D & atom_H );

// Add two hydrogen atoms to an sp2 nitrogen atom
// Returns the coordinates of the two hydrogen atoms
// atom_1 and atom_2 are required to define the plan in which the two atoms should lie.
// Everything in Cartesian coordinates
std::vector< Vector3D > add_2_hydrogen_atoms_to_sp2_nitrogen( const Vector3D & atom_N, const Vector3D & atom_bound_to_N, const Vector3D & atom_1, const Vector3D & atom_2 );

Vector3D add_hydrogen_atom_to_sp2_atom( const Vector3D & central_atom, const Element element_central_atom, const Vector3D & neighbour_1, const Vector3D & neighbour_2 );

Vector3D add_hydrogen_atom_to_sp3_atom( const Vector3D & central_atom, const Element element_central_atom, const Vector3D & neighbour_1, const Vector3D & neighbour_2, const Vector3D & neighbour_3 );

std::vector< Vector3D > add_2_hydrogen_atoms_to_sp3_atom( const Vector3D & central_atom, const Element element_central_atom, const Vector3D & neighbour_1, const Vector3D & neighbour_2 );

void normalise_X_H_bonds( CrystalStructure & crystal_structure );

void normalise_C_F_bonds( CrystalStructure & crystal_structure );

void normalise_Cu_Cl_bonds( CrystalStructure & crystal_structure );

#endif // MODELBUILDING_H

