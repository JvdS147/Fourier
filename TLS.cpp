/* *********************************************
Copyright (c) 2013-2023, Cornelis Jan (Jacco) van de Streek
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

#include "TLS.h"
#include "Angle.h"
#include "3DCalculations.h"
#include "AnisotropicDisplacementParameters.h"
#include "BasicMathsFunctions.h"

#include <stdexcept>

// ********************************************************************************

TLS::TLS():
is_on_inversion_(false)
{
}

// ********************************************************************************

TLS::TLS( const SymmetricMatrix3D & T, const SymmetricMatrix3D & L, const Matrix3D & S ):
T_(T),
L_(L),
S_(S),
is_on_inversion_(false)
{
    correction_matrix_ = SymmetricMatrix3D() + 0.5 * ( L_.trace() * SymmetricMatrix3D() - L_ );
}

// ********************************************************************************

void TLS::set_is_on_inversion( const bool value )
{
    is_on_inversion_ = value;
    if ( ! is_on_inversion_ )
        return;
    if ( ! nearly_zero( origin_ ) )
        throw std::runtime_error( "TLS::set_is_on_inversion( true ): origin is not the zero vector." );
    origin_ = Vector3D();
    // S is now the zero matrix
    if ( ! nearly_equal( S_, Matrix3D( 0.0 ) ) )
        throw std::runtime_error( "TLS::set_is_on_inversion( true ): S is not the zero matrix." );
    S_ = Matrix3D( 0.0 );
}

// ********************************************************************************

// We have to be very careful here: both r and the TLS matrices are defined with respect to a Cartesian basis, but
// of course it must be the same Cartesian basis.
AnisotropicDisplacementParameters TLS::U( const Vector3D & r ) const
{
    double x = r.x();
    double y = r.y();
    double z = r.z();
    Matrix3D A( 0.0,   z,  -y,
                 -z, 0.0,   x,
                  y, - x, 0.0 );
    // 1.) and 2.) should give the same result
    // 1.)
    SymmetricMatrix3D U_1 = T_ +  Matrix3D2SymmetricMatrix3D( A*S_ + transpose( A*S_ ) + A*L_*transpose( A ) );
    // 2.)
    double u11 = T_.value( 0, 0 ) + L_.value( 1, 1 )*square(z) + L_.value( 2, 2 )*square(y) - 2.0*y*z*L_.value( 1, 2 ) - 2.0*y*S_.value( 2, 0 ) + 2.0*z*S_.value( 1, 0 ) ;
    double u22 = T_.value( 1, 1 ) + L_.value( 0, 0 )*square(z) + L_.value( 2, 2 )*square(x) - 2.0*x*z*L_.value( 0, 2 ) - 2.0*z*S_.value( 0, 1 ) + 2.0*x*S_.value( 2, 1 ) ;
    double u33 = T_.value( 2, 2 ) + L_.value( 0, 0 )*square(y) + L_.value( 1, 1 )*square(x) - 2.0*x*y*L_.value( 0, 1 ) - 2.0*x*S_.value( 1, 2 ) + 2.0*y*S_.value( 0, 2 ) ;
    double u23 = T_.value( 1, 2 ) + -y*z*L_.value( 0, 0 ) - square(x)*L_.value( 1, 2 ) + x*y*L_.value( 0, 2 ) + x*z*L_.value( 0, 1 ) - x*S_.value( 1, 1 ) + x*S_.value( 2, 2 ) + y*S_.value( 0, 1 ) - z*S_.value( 0, 2 ) ;
    double u13 = T_.value( 0, 2 ) + -x*z*L_.value( 1, 1 ) + x*y*L_.value( 1, 2 ) - square(y)*L_.value( 0, 2 ) + y*z*L_.value( 0, 1 ) + y*S_.value( 0, 0 ) - y*S_.value( 2, 2 ) + z*S_.value( 1, 2 ) - x*S_.value( 1, 0 ) ;
    double u12 = T_.value( 0, 1 ) + -x*y*L_.value( 2, 2 ) + x*z*L_.value( 1, 2 ) + y*z*L_.value( 0, 2 ) - square(z)*L_.value( 0, 1 ) - z*S_.value( 0, 0 ) + z*S_.value( 1, 1 ) + x*S_.value( 2, 0 ) - y*S_.value( 2, 1 ) ;
    SymmetricMatrix3D U_2( u11, u22, u33, u12, u13, u23 );
    if ( ! nearly_equal( U_1, U_2 ) )
        throw std::runtime_error( "TLS::U( Vector3D ): Error U_1 != U_2." );
    return AnisotropicDisplacementParameters( U_1 );
}

// ********************************************************************************

Vector3D TLS::centre_of_reaction() const
{
    Matrix3D L( L_.value( 1, 1 ) + L_.value( 2, 2 ), -L_.value( 0, 1 ),                    -L_.value( 0, 2 ),
               -L_.value( 0, 1 ),                     L_.value( 0, 0 ) + L_.value( 2, 2 ), -L_.value( 1, 2 ),
               -L_.value( 0, 2 ),                    -L_.value( 1, 2 ),                     L_.value( 0, 0 ) + L_.value( 1, 1 ) );
    L.invert();
    Vector3D result = L * Vector3D( S_.value( 1, 2 ) - S_.value( 2, 1 ), S_.value( 2, 0 ) - S_.value( 0, 2 ), S_.value( 0, 1 ) - S_.value( 1, 0 ) );
    return result;
}

// ********************************************************************************

Angle TLS::libration() const
{
    return Angle::from_radians( sqrt( L_.trace() ) );
}

// ********************************************************************************

// @@ Could move this to the TOPAS class?
std::vector< std::string > TLS::TOPAS_lines( const Vector3D & r ) const
{
    std::vector< std::string > result;
    if ( is_on_inversion_ )
    {
    }
    else
    {
    }
    return result;
//    u11" + label + " = L22*rz" + label + "^2 + L33*ry" + label + "^2 - 2*ry" + label + "*rz" + label + "*L23 - 2*ry" + label + "*S31 + 2*rz" + label + "*S21 + T11; : 0.0" );
//    u22" + label + " = L11*rz" + label + "^2 + L33*rx" + label + "^2 - 2*rx" + label + "*rz" + label + "*L13 - 2*rz" + label + "*S12 + 2*rx" + label + "*S32 + T22; : 0.0" );
//    u33" + label + " = L11*ry" + label + "^2 + L22*rx" + label + "^2 - 2*rx" + label + "*ry" + label + "*L12 - 2*rx" + label + "*S23 + 2*ry" + label + "*S13 + T33; : 0.0" );
//    u12" + label + " = -rx" + label + "*ry" + label + "*L33 + rx" + label + "*rz" + label + "*L23 + ry" + label + "*rz" + label + "*L13 - rz" + label + "^2*L12 - rz" + label +
//                                    "*S11 + rz" + label + "*S22 + rx" + label + "*S31 - ry" + label + "*S32 + T12; : 0.0" );
//    u13" + label + " = -rx" + label + "*rz" + label + "*L22 + rx" + label + "*ry" + label + "*L23 - ry" + label + "^2*L13 + ry" + label + "*rz" + label + "*L12 + ry" + label +
//                                     "*S11 - ry" + label + "*S33 + rz" + label + "*S23 - rx" + label + "*S21 + T13; : 0.0" );
//    u23" + label + " = -ry" + label + "*rz" + label + "*L11 - rx" + label + "^2*L23 + rx" + label + "*ry" + label + "*L13 + rx" + label + "*rz" + label + "*L12 - rx" + label +
//                                     "*S22 + rx" + label + "*S33 + ry" + label + "*S12 - rz" + label + "*S13 + T23; : 0.0" );
}

// ********************************************************************************

Vector3D TLS::corrected_relative_Cartesian_coordinate( const Vector3D & r ) const
{
    return correction_matrix_ * r;
}

// ********************************************************************************

//double TLS::corrected_distance( const Vector3D & lhs, const Vector3D & rhs ) const
//{
//
//    Vector3D dr( 0.5*( ( L_.value( 1, 1 ) + L_.value( 2, 2 ) )*rx                         -L_.value( 0, 1 )*ry                         -L_.value( 0, 2 )*rz ),
//                 0.5*(                       -L_.value( 0, 1 )*rx + ( L_.value( 0, 0 ) + L_.value( 2, 2 ) )*ry                         -L_.value( 1, 2 )*rz ),
//                 0.5*(                       -L_.value( 0, 2 )*rx                         -L_.value( 1, 2 )*ry + ( L_.value( 0, 0 ) + L_.value( 1, 1 ) )*rz ) );
//    return r + dr;
//}

// ********************************************************************************
