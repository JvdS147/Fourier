#ifndef CRYSTALLATTICE_H
#define CRYSTALLATTICE_H

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

#include "Angle.h"
#include "Matrix3D.h"
#include "SymmetricMatrix3D.h"
#include "Vector3D.h"

#include <string>

// a along x, b in xy plane, right-handed coordinate frame
// We also abuse this class for any functionality related to parallelepipeds
// The CrystalLattice and the SpaceGroup have a mutual dependence that is currently not captured anywhere.
class CrystalLattice
{
public:

    CrystalLattice();

    CrystalLattice( const double a,
                    const double b,
                    const double c,
                    const Angle alpha,
                    const Angle beta,
                    const Angle gamma );

    // Initialises from a CASTEP LATTICE_CART block in a .cell file.
    void from_CASTEP( const Matrix3D & matrix );

    // We need a() as a vector and a() as a length. We cannot overload by return type,
    // so they must have different names
    double a() const { return a_; }
    double b() const { return b_; }
    double c() const { return c_; }
    Angle alpha() const { return alpha_; }
    Angle beta()  const { return beta_; }
    Angle gamma() const { return gamma_; }

    // Always <= 1.0, the smaller, the more acute the angles.
    double orthogonality_defect() const;

    Vector3D lattice_vector( const size_t index ) const;

    // a along x
    Vector3D a_vector() const { return a_vector_; }
    Vector3D b_vector() const { return b_vector_; }
    Vector3D c_vector() const { return c_vector_; }
    double a_star() const { return a_star_; }
    double b_star() const { return b_star_; }
    double c_star() const { return c_star_; }
    Angle alpha_star() const { return alpha_star_; }
    Angle beta_star()  const { return beta_star_; }
    Angle gamma_star() const { return gamma_star_; }
    Vector3D a_star_vector() const { return a_star_vector_; }
    Vector3D b_star_vector() const { return b_star_vector_; }
    Vector3D c_star_vector() const { return c_star_vector_; }
    double volume() const { return volume_; }
    Matrix3D fractional_to_orthogonal_matrix() const { return fractional_to_orthogonal_matrix_; }
    Matrix3D orthogonal_to_fractional_matrix() const { return orthogonal_to_fractional_matrix_; }

    void enclosing_box( Vector3D & min_min_min, Vector3D & max_max_max ) const;

    Matrix3D metric_matrix() const;

    // Fractional to Cartesian, c along z *and* transposed
    Matrix3D for_CASTEP() const;

    Vector3D orthogonal_to_fractional( const Vector3D & input ) const;
    Vector3D fractional_to_orthogonal( const Vector3D & input ) const;

    // This is the matrix N as used by Grosse-Kunstleve to convert U_cif to U_star
    SymmetricMatrix3D N() const { return N_; }
    SymmetricMatrix3D N_inverse() const { return N_inverse_; }

    // Rescales a, b and c isotropically so that the new unit cell volume becomes
    // equal to the specified target_volume. alpha, beta and gamma are not changed.
    // If Z is specified for the target_volume, tries to guess Z from the current unit-cell volume and
    // adjusts accordingly.
    void rescale_volume( const double target_volume, size_t Z = 0 );

    // Finds shortest distance, in Angstrom, between two positions given in fractional coordinates.
    double shortest_distance( const Vector3D & lhs, const Vector3D & rhs ) const;

    // Finds shortest distance^2, in Angstrom^2, between two positions given in fractional coordinates.
    double shortest_distance2( const Vector3D & lhs, const Vector3D & rhs ) const;

    // Finds shortest distance, in Angstrom, between two positions given in fractional coordinates.
    // Returns the shortest distance (in Angstrom) and the shortest difference vector (defined as rhs - lhs, in fractional coordinates).
    void shortest_distance( const Vector3D & lhs, const Vector3D & rhs, double & distance, Vector3D & difference_vector ) const;

    enum LatticeSystem { TRICLINIC, MONOCLINIC_A, MONOCLINIC_B, MONOCLINIC_C, ORTHORHOMBIC, TRIGONAL, TETRAGONAL, HEXAGONAL, RHOMBOHEDRAL, CUBIC };

    // The lattice system is initialised by deducing it from the unit-cell parameters
    LatticeSystem lattice_system() const { return lattice_system_; }
    void set_lattice_system( const LatticeSystem lattice_system ) { lattice_system_ = lattice_system; }

    // There are many things we can do to make a lattice more "beautiful": order a, b and c by size,
    // make alpha, beta and gamma all greater than or all less than 90.
    // And the transformation matrix that is chosen can also be chosen from
    // many candidates, we can e.g. choose the one that is as close as possible to the identity matrix.
    Matrix3D choose_angles_close_to_90() const;

    void transform( const Matrix3D & transformation_matrix );

    // Prints a, b, c, alpha, beta, gamma
    void print() const;

    // Shows the full vectors and reciprocal vectors
    void show() const; // For debugging

    // Some matrices necessary for the formulae given by R. T. Downs / G. V. Gibbs

    Matrix3D Downs_D() const;
    Matrix3D Downs_D_star() const;
    Matrix3D Downs_G() const;
    Matrix3D Downs_G_star() const;

private:
    double a_;
    double b_;
    double c_;
    Angle alpha_;
    Angle beta_;
    Angle gamma_;
    double a_star_;
    double b_star_;
    double c_star_;
    Angle alpha_star_;
    Angle beta_star_;
    Angle gamma_star_;
    Vector3D a_vector_;
    Vector3D b_vector_;
    Vector3D c_vector_;
    Vector3D a_star_vector_;
    Vector3D b_star_vector_;
    Vector3D c_star_vector_;
    double volume_;
    Matrix3D fractional_to_orthogonal_matrix_;
    Matrix3D orthogonal_to_fractional_matrix_;
    SymmetricMatrix3D N_;
    SymmetricMatrix3D N_inverse_;
    LatticeSystem lattice_system_;
};

// Deduces the lattice system based on the unit-cell parameters.
CrystalLattice::LatticeSystem deduce_lattice_system( const CrystalLattice & crystal_lattice );

std::string LatticeSystem2string( const CrystalLattice::LatticeSystem lattice_system );

CrystalLattice average( const CrystalLattice & lhs, const CrystalLattice & rhs );

CrystalLattice average( const CrystalLattice & cl_1, const CrystalLattice & cl_2, const CrystalLattice & cl_3 );

// This would be an alternative solution. The main question is: does the weight more naturally apply to lhs or to rhs?
// CrystalLattice average( const CrystalLattice & lhs, const CrystalLattice & rhs, const double weight = 1.0 );

// Length tolerance is relative, angle tolerance is absolute
bool nearly_equal( const CrystalLattice & lhs, const CrystalLattice & rhs, double length_tolerance_percentage = 10.0, const Angle angle_tolerance = Angle::from_degrees( 10.0 ) );

#endif // CRYSTALLATTICE_H

