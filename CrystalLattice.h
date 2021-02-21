#ifndef CRYSTALLATTICE_H
#define CRYSTALLATTICE_H

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

#include "Vector3D.h"
#include "Matrix3D.h"
#include "Angle.h"

#include <string>

// a along x, b in xy plane, right-handed coordinate frame
// We also abuse this class for any functionality related to parallelepipeds
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

    // We need a() as a vector and a() as a length. We cannot overload by return type,
    // so they must have different names
    double a() const { return a_; }
    double b() const { return b_; }
    double c() const { return c_; }
    Angle alpha() const { return alpha_; }
    Angle beta()  const { return beta_; }
    Angle gamma() const { return gamma_; }
    
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

    enum LatticeSystem { TRICLINIC, MONOCLINIC, ORTHORHOMBIC, TRIGONAL, TETRAGONAL, HEXAGONAL, RHOMBOHEDRAL, CUBIC };
    
    // The lattice system is initialised by deducing it from the unit-cell parameters
    LatticeSystem lattice_system() const { return lattice_system_; }
    void set_lattice_system( const LatticeSystem lattice_system ) { lattice_system_ = lattice_system; }

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
    LatticeSystem lattice_system_;
};

// Deduces the lattice system based on the unit-cell parameters.
CrystalLattice::LatticeSystem deduce_lattice_system( const CrystalLattice & crystal_lattice );

std::string LatticeSystem2string( const CrystalLattice::LatticeSystem lattice_system );

#endif // CRYSTALLATTICE_H

