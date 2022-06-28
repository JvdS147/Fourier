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

#include "AnisotropicDisplacementParameters.h"
#include "3DCalculations.h"
#include "CrystalLattice.h"
#include "Matrix3D.h"
#include "Vector3D.h"

#include <iostream>
#include <stdexcept>

// ********************************************************************************

AnisotropicDisplacementParameters::AnisotropicDisplacementParameters()
{
    data_ = SymmetricMatrix3D( 0.0 );
}

// ********************************************************************************

AnisotropicDisplacementParameters::AnisotropicDisplacementParameters( const std::vector< Vector3D > & points )
{
    data_ = covariance_matrix( points );
}

// ********************************************************************************

AnisotropicDisplacementParameters::AnisotropicDisplacementParameters( const SymmetricMatrix3D & matrix )
{
    data_ = matrix;
}

// ********************************************************************************

AnisotropicDisplacementParameters::AnisotropicDisplacementParameters( const double u_iso )
{
    data_ = SymmetricMatrix3D( u_iso );
}

// ********************************************************************************

double AnisotropicDisplacementParameters::value( size_t i, size_t j ) const
{
    if ( ( 2 < i ) || ( 2 < j ) )
        throw std::runtime_error( "AnisotropicDisplacementParameters::value(): index out of bounds." );
    return data_.value( i, j );
}

// ********************************************************************************

void AnisotropicDisplacementParameters::set_value( size_t i, size_t j, const double value )
{
    if ( ( 2 < i ) || ( 2 < j ) )
        throw std::runtime_error( "AnisotropicDisplacementParameters::value(): index out of bounds." );
    data_.set_value( i, j, value );
}

// ********************************************************************************

SymmetricMatrix3D AnisotropicDisplacementParameters::U_star( const CrystalLattice & crystal_lattice ) const
{
    Matrix3D result = crystal_lattice.orthogonal_to_fractional_matrix() * U_cart() * transpose( crystal_lattice.orthogonal_to_fractional_matrix() );
    return Matrix3D2SymmetricMatrix3D( result );
}

// ********************************************************************************

SymmetricMatrix3D AnisotropicDisplacementParameters::U_cart() const
{
    return data_;
}

// ********************************************************************************

SymmetricMatrix3D AnisotropicDisplacementParameters::U_cif( const CrystalLattice & crystal_lattice ) const
{
    return Matrix3D2SymmetricMatrix3D( crystal_lattice.N_inverse() * U_star( crystal_lattice ) * crystal_lattice.N_inverse() );
}

// ********************************************************************************

double AnisotropicDisplacementParameters::U_iso() const
{
    return U_cart().trace() / 3.0;
}

// ********************************************************************************

double AnisotropicDisplacementParameters::average_displacement_squared( const Vector3D & v, const CrystalLattice & crystal_lattice )
{
    return ( v * ( transpose( crystal_lattice.Downs_G() ) * U_star( crystal_lattice ) * crystal_lattice.Downs_G() ) * v ) /
        ( v * crystal_lattice.Downs_G() * v );
}

// ********************************************************************************

void AnisotropicDisplacementParameters::show() const
{
    std::cout << "Uiso = " << U_iso() << std::endl;
    data_.show();
}

// ********************************************************************************

AnisotropicDisplacementParameters transform_adps( const AnisotropicDisplacementParameters & ADPs, Matrix3D transformation, CrystalLattice crystal_lattice )
{
    SymmetricMatrix3D U_star = ADPs.U_star( crystal_lattice );
    crystal_lattice.transform( transformation );
    transformation.invert();
    transformation.transpose();
    // The actual transformation of the ADPs
    U_star = Matrix3D2SymmetricMatrix3D( transformation * U_star * transpose( transformation ) );
    return AnisotropicDisplacementParameters( U_star_2_U_cart( U_star, crystal_lattice ) );
}

// ********************************************************************************

// The following is what is needed to rotate the ADPs as part of a symmetry operation
AnisotropicDisplacementParameters rotate_adps( const AnisotropicDisplacementParameters & ADPs, const Matrix3D & rotation, const CrystalLattice & crystal_lattice )
{
    SymmetricMatrix3D U_star = ADPs.U_star( crystal_lattice );
    U_star = Matrix3D2SymmetricMatrix3D( rotation * U_star * transpose( rotation ) );
    return AnisotropicDisplacementParameters( U_star_2_U_cart( U_star, crystal_lattice ) );
}

// ********************************************************************************

SymmetricMatrix3D U_star_2_U_cart( const SymmetricMatrix3D & U_star, const CrystalLattice & crystal_lattice )
{
    Matrix3D A = crystal_lattice.fractional_to_orthogonal_matrix();
    return Matrix3D2SymmetricMatrix3D( A * U_star * transpose( A ) );
}

// ********************************************************************************

SymmetricMatrix3D U_star_2_U_cif ( const SymmetricMatrix3D & U_star, const CrystalLattice & crystal_lattice )
{
    return Matrix3D2SymmetricMatrix3D( crystal_lattice.N_inverse() * U_star * crystal_lattice.N_inverse() );
}

// ********************************************************************************

SymmetricMatrix3D U_cart_2_U_cif ( const SymmetricMatrix3D & U_cart, const CrystalLattice & crystal_lattice )
{
    return ( U_star_2_U_cif( U_cart_2_U_star( U_cart, crystal_lattice ), crystal_lattice ) );
}

// ********************************************************************************

SymmetricMatrix3D U_cart_2_U_star( const SymmetricMatrix3D & U_cart, const CrystalLattice & crystal_lattice )
{
    Matrix3D A = crystal_lattice.orthogonal_to_fractional_matrix();
    return Matrix3D2SymmetricMatrix3D( A * U_cart * transpose( A ) );
}

// ********************************************************************************

SymmetricMatrix3D U_cif_2_U_star( const SymmetricMatrix3D & U_cif, const CrystalLattice & crystal_lattice )
{
    return Matrix3D2SymmetricMatrix3D( crystal_lattice.N() * U_cif * crystal_lattice.N() );
}

// ********************************************************************************

SymmetricMatrix3D U_cif_2_U_cart( const SymmetricMatrix3D & U_cif, const CrystalLattice & crystal_lattice )
{
    return ( U_star_2_U_cart( U_cif_2_U_star( U_cif, crystal_lattice ), crystal_lattice ) );
}

// ********************************************************************************

