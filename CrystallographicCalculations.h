#ifndef CRYSTALLOGRAPHICCALCULATIONS_H
#define CRYSTALLOGRAPHICCALCULATIONS_H

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

// A major reason for the existence of this file is that it collects all functions that combine
// Centring, MillerIndices, SymmetryOperator and SpaceGroup, so that those classes do not need to know about each other.

//#include "BasicMathsFunctions.h"

class AnisotropicDisplacementParameters;
class Centring;
//class CollectionOfPoints;
//class CoordinateFrame;
class CrystalLattice;
class CrystalStructure;
class Matrix3D;
class MillerIndices;
class NormalisedVector3D;
class PointGroup;
class SpaceGroup;
class SymmetricMatrix3D;
class SymmetryOperator;
//class Vector2D;
class Vector3D;

#include <string>
#include <vector>

//struct SpecialPositionsReport
//{
//    size_t number_of_atoms_in_unit_cell_; // Could make it a ChemicalFormula.
//    // We could enumerate the number of atoms found per point group, but two atoms that are on two special positions with the same point group are not necesarily the same atom.
//    std::vector< size_t > point_group_orders_;
//    std::vector< size_t > natoms_per_point_group_order_; // Same size as point_group_orders_
//    size_t natoms_on_special_positions_; // Sum of natoms_per_point_group_order_
//    size_t natoms_with_symmetry_copies; // Benzene in P-1 would return six (three carbon, three hydrogen).
//    size_t nsymmetry_operators_;
//    size_t nmolecules_;
//};

//SpecialPositionsReport special_positions_report()

// In a orthorhombic space groups we can have a couple of permutations for the axes.
// There are three axes, so there should be 3! = 6 permutations, I think.
// They are abc, cab, bca, ba-c, c-ba, -acb.
// This function returns these six transformation matrices.
// The first transformation matrix is the identity matrix.
std::vector< Matrix3D > orthorhombic_unit_cell_axes_permutations();

// Returns abc, cab, bca, ba-c, c-ba, -acb.
std::vector< std::string > orthorhombic_unit_cell_axes_permutation_labels();

Vector3D reciprocal_lattice_point( const MillerIndices miller_indices, const CrystalLattice & crystal_lattice );

NormalisedVector3D reciprocal_lattice_direction( const MillerIndices miller_indices, const CrystalLattice & crystal_lattice );

//MillerIndices operator*( const Matrix3D & matrix, const MillerIndices & miller_indices );
MillerIndices operator*( const MillerIndices & miller_indices, const Matrix3D & matrix );
MillerIndices operator*( const MillerIndices & miller_indices, const SymmetricMatrix3D & matrix );

double operator*( const MillerIndices & miller_indices, const Vector3D & vector_3D );

// The transformation has already been applied to the space group.
// If e.g. the transformation matrix codes P-1 to I-1, it adds [ 1/2, 1/2, 1/2 ] to the symmetry operators.
void add_centring_to_space_group_after_transformation( Matrix3D tranformation_matrix, SpaceGroup & space_group );

std::vector< SymmetryOperator > centring_generators( const Centring & centring );

// @@ Currently we can only cope with whatever was generated by std::vector< SymmetryOperator > centring_generators( const Centring & centring ).
Centring expand_centring_generators( const std::vector< SymmetryOperator > & generators );

bool nearly_equal( const CrystalStructure & lhs, const CrystalStructure & rhs );

// Takes and ADP matrix which may have been affected by e.g. rounding errors and adjusts it to
// comply with the site symmetry.
AnisotropicDisplacementParameters adjust_to_site_symmetry( const AnisotropicDisplacementParameters & adps, const PointGroup & point_group, const CrystalLattice & crystal_lattice );

// Throws if cubic.
MillerIndices select_realistic_preferred_orientation_direction( const CrystalLattice & crystal_lattice );

// Throws if cubic.
MillerIndices select_realistic_preferred_orientation_direction( const CrystalLattice & crystal_lattice, const SpaceGroup & space_group );

#endif // CRYSTALLOGRAPHICCALCULATIONS_H

