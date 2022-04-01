#ifndef DCALCULATIONS_H
#define DCALCULATIONS_H

/* *********************************************
Copyright (c) 2013-2022, Cornelis Jan (Jacco) van de Streek
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
// Matrix3D, SymmetricMatrix3D, NormalisedVector3D, Vector2D and Vector3D, so that those classes do not need to know about each other.

#include "BasicMathsFunctions.h"

class Angle;
class CollectionOfPoints;
class CoordinateFrame;
class CrystalLattice;
class CrystalStructure;
class Matrix3D;
class MillerIndices;
class NormalisedVector3D;
class Plane;
class SpaceGroup;
class SymmetricMatrix3D;
class Vector2D;
class Vector3D;

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

Vector3D reciprocal_lattice_point( const MillerIndices miller_indices, const CrystalLattice & crystal_lattice );

NormalisedVector3D reciprocal_lattice_direction( const MillerIndices miller_indices, const CrystalLattice & crystal_lattice );

// Gram-Schmidt orthogonalisation
NormalisedVector3D orthogonalise( const NormalisedVector3D & n, const Vector3D & r );

Vector3D change_of_basis( const Vector3D & point, const CoordinateFrame & before, const CoordinateFrame & after );

// Ambiguity here: do we use the original points or the points w.r.t. the centre of mass?
Plane plane( const CollectionOfPoints & points );

// The basis in the plane is generated with Plane::coordinate_frame().
Vector2D projection( const Plane & plane, const Vector3D & point );

double root_mean_square_devation_from_mean_plane( const CollectionOfPoints & points, const Plane & plane );

double root_mean_square_devation_from_mean_plane( const std::vector< Vector3D > & points, const Plane & plane );

// Easier: apply rotational part of symmetry operator to Ucif
//AnisotropicDisplacementParameters operator*( const SymmetryOperator & lhs, const AnisotropicDisplacementParameters rhs );

// Try to use the overload with CollectionOfPoints.
SymmetricMatrix3D covariance_matrix( const std::vector< Vector3D > & points );

SymmetricMatrix3D covariance_matrix( const CollectionOfPoints & points );

// This implicitly assumes that the vector is in Cartesian coordinates, not fractional coordinates.
NormalisedVector3D normalised_vector( const Vector3D & vector3D );
Vector3D NormalisedVector3D2Vector3D( const NormalisedVector3D & rhs );

Vector3D operator+( const NormalisedVector3D & lhs, const NormalisedVector3D & rhs );
Vector3D operator-( const NormalisedVector3D & lhs, const NormalisedVector3D & rhs );
Vector3D operator-( const NormalisedVector3D & lhs, const Vector3D & rhs );

Vector3D operator*( const double lhs, const NormalisedVector3D & rhs );
Vector3D operator*( const NormalisedVector3D & lhs, const double rhs );
Vector3D operator/( const NormalisedVector3D & lhs, const double rhs );

// Convert SymmetricMatrix3D to a Matrix3D.
Matrix3D SymmetricMatrix3D2Matrix3D( const SymmetricMatrix3D & matrix );

// Convert Matrix3D to a SymmetricMatrix3D.
// Throws if the matrix was not symmetric within tolerance.
SymmetricMatrix3D Matrix3D2SymmetricMatrix3D( const Matrix3D & matrix, const double tolerance = TOLERANCE );

Vector3D operator*( const Matrix3D & matrix, const Vector3D & vector );
Vector3D operator*( const SymmetricMatrix3D & matrix, const Vector3D & vector );

// The transposition is implied.
Vector3D operator*( const Vector3D & vector, const Matrix3D & matrix );
Vector3D operator*( const Vector3D & vector, const SymmetricMatrix3D & matrix );

Vector3D operator*( const NormalisedVector3D & vector, const Matrix3D & matrix );
Vector3D operator*( const NormalisedVector3D & vector, const SymmetricMatrix3D & matrix );

Matrix3D operator*( const SymmetricMatrix3D & lhs, const SymmetricMatrix3D & rhs );
Matrix3D operator*( const Matrix3D & lhs, const SymmetricMatrix3D & rhs );
Matrix3D operator*( const SymmetricMatrix3D & lhs, const Matrix3D & rhs );
Matrix3D operator+( const Matrix3D & lhs, const SymmetricMatrix3D & rhs );
Matrix3D operator+( const SymmetricMatrix3D & lhs, const Matrix3D & rhs );
Matrix3D operator-( const Matrix3D & lhs, const SymmetricMatrix3D & rhs );
Matrix3D operator-( const SymmetricMatrix3D & lhs, const Matrix3D & rhs );

// The transposition is implied.
double operator*( const NormalisedVector3D & lhs, const Vector3D & rhs );
// The transposition is implied.
double operator*( const Vector3D & lhs, const NormalisedVector3D & rhs );

//MillerIndices operator*( const Matrix3D & matrix, const MillerIndices & miller_indices );
MillerIndices operator*( const MillerIndices & miller_indices, const Matrix3D & matrix );
MillerIndices operator*( const MillerIndices & miller_indices, const SymmetricMatrix3D & matrix );

double operator*( const MillerIndices & miller_indices, const Vector3D & vector_3D );

// Counter-clockwise, rotation axis coming out of the plane of the paper.
Matrix3D rotation_about_x( const Angle angle );
Matrix3D rotation_about_y( const Angle angle );
Matrix3D rotation_about_z( const Angle angle );

// "origin" is a point on the axis, "n" is the direction of the axis.
Vector3D rotate_point_about_axis( Vector3D point, const Vector3D & origin, const NormalisedVector3D & n, const Angle angle );

Vector3D cylindrical2Cartesian( const double r, Angle phi, const double z );

bool are_translationally_equivalent( const double x, const double y );

bool are_translationally_equivalent( const Vector3D & lhs, const Vector3D & rhs );

double adjust_for_translations( const double input );

// This should probably change in place,
// but our Vector3D class does not allow for its elements to be addressed that way.
Vector3D adjust_for_translations( const Vector3D & input );

// Returns the (smaller) angle between two vectors.
// Throws if at least one of the vectors is the zero vector.
Angle angle( const Vector3D & lhs, const Vector3D & rhs );

// Returns the (smaller) angle between two vectors.
Angle angle( const NormalisedVector3D & lhs, const Vector3D & rhs );

// Returns the (smaller) angle between two vectors.
// The length is the length of the vector rhs
Angle angle( const NormalisedVector3D & lhs, const Vector3D & rhs, const double length );

// Returns the (smaller) angle between two vectors.
Angle angle( const NormalisedVector3D & lhs, const NormalisedVector3D & rhs );

// Returns the (smaller) angle between two planes.
Angle angle( const Plane & lhs, const Plane & rhs );

// The order of the points determines the sign of the torsion.
Angle signed_torsion( const Vector3D & r1, const Vector3D & r2, const Vector3D & r3, const Vector3D & r4 );

Matrix3D A_centred_to_primitive();
Matrix3D B_centred_to_primitive();
Matrix3D C_centred_to_primitive();
Matrix3D I_centred_to_primitive();
Matrix3D F_centred_to_primitive();
Matrix3D R_centred_to_primitive();

bool nearly_equal( const CrystalStructure & lhs, const CrystalStructure & rhs );

// The transformation has already been applied to the space group.
// If e.g. the transformation matrix codes P-1 to I-1, it adds [ 1/2, 1/2, 1/2 ] to the symmetry operators.
void add_centring_to_space_group_after_transformation( Matrix3D tranformation_matrix, SpaceGroup & space_group );

#endif // DCALCULATIONS_H

