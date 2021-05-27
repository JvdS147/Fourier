#ifndef ANISOTROPICDISPLACEMENTPARAMETERS_H
#define ANISOTROPICDISPLACEMENTPARAMETERS_H

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

class CrystalLattice;
class Matrix3D;
class Vector3D;

#include "SymmetricMatrix3D.h"

#include <vector>

/*
  Anisotropic Displacement Parameters
  
  Note that the use of beta and B is deprecated, please use U.
*/
class AnisotropicDisplacementParameters
{
public:

    // Default constructor
    AnisotropicDisplacementParameters();

    // Must be points (atomic coordinates) in Cartesian coordinates, not necessarily centred around the origin.
    AnisotropicDisplacementParameters( const std::vector< Vector3D > & points );

    // Expects U_cart
    AnisotropicDisplacementParameters( const SymmetricMatrix3D & matrix );

    explicit AnisotropicDisplacementParameters( const double u_iso );

    // In keeping with the silly C++ convention: zero-based
    // Returns U_cart, which is independent of the unit cell.
    double value( size_t i, size_t j ) const { return data_.value( i, j ); }

    // Requires U_cart, which is independent of the unit cell.
    void set_value( size_t i, size_t j, const double value ) { data_.set_value( i, j, value ); }

    SymmetricMatrix3D U_star( const CrystalLattice & crystal_lattice ) const;
    
    // beta = 2*pi^2 * U_star
//    SymmetricMatrix3D beta( const CrystalLattice & crystal_lattice ) const;
    SymmetricMatrix3D U_cart() const;
    SymmetricMatrix3D U_cif( const CrystalLattice & crystal_lattice ) const;

    double U_iso() const;
    
    // Returns the square of the average displacement along a vector.
    // The vector is in fractional coordinates.
    // Units: Angstrom^2.
    double average_displacement_squared( const Vector3D & v, const CrystalLattice & crystal_lattice );
    
    void show() const;
    
private:
    // Because this class will be used for analysing MD trajectories, we will have
    // very many instances. It will always take a while to do the analysis,
    // so it does not matter if the class is a bit slower than it could be, but
    // we can easily run out of memory so store as little as possible.
    SymmetricMatrix3D data_; // Stores U_cart, the only one for which no knowledge of the unit cell is needed.

};

// The following is what is needed to transform the ADPs when the unit cell is transformed with a transformation matrix
AnisotropicDisplacementParameters transform_adps( const AnisotropicDisplacementParameters & ADPs, Matrix3D transformation, CrystalLattice crystal_lattice );

// The following is what is needed to rotate the ADPs as part of a symmetry operation
AnisotropicDisplacementParameters rotate_adps( const AnisotropicDisplacementParameters & ADPs, const Matrix3D & rotation, const CrystalLattice & crystal_lattice );

// Of course, this should be done inside the class, but then we must store a CrystalLattice in the class
// and we do not have smart pointers so we must store a copy.
SymmetricMatrix3D U_star_2_U_cart( const SymmetricMatrix3D & U_star, const CrystalLattice & crystal_lattice );
SymmetricMatrix3D U_star_2_U_cif ( const SymmetricMatrix3D & U_star, const CrystalLattice & crystal_lattice );
SymmetricMatrix3D U_cart_2_U_cif ( const SymmetricMatrix3D & U_cart, const CrystalLattice & crystal_lattice );
SymmetricMatrix3D U_cart_2_U_star( const SymmetricMatrix3D & U_cart, const CrystalLattice & crystal_lattice );
SymmetricMatrix3D U_cif_2_U_star ( const SymmetricMatrix3D & U_cif , const CrystalLattice & crystal_lattice );
SymmetricMatrix3D U_cif_2_U_cart ( const SymmetricMatrix3D & U_cif , const CrystalLattice & crystal_lattice );

#endif // ANISOTROPICDISPLACEMENTPARAMETERS_H

