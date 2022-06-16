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

#include "3DCalculations.h"
#include "Angle.h"
#include "BasicMathsFunctions.h"
#include "Centring.h"
#include "CollectionOfPoints.h"
#include "CoordinateFrame.h"
#include "CrystalLattice.h"
#include "Matrix3D.h"
#include "MillerIndices.h"
//#include "NormalisedVector3D.h"
#include "Plane.h"
#include "SpaceGroup.h"
#include "SymmetricMatrix3D.h"
#include "Vector2D.h"
#include "Vector3DCalculations.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

// ********************************************************************************

Vector3D reciprocal_lattice_point( const MillerIndices miller_indices, const CrystalLattice & crystal_lattice )
{
    return ( miller_indices.h() * crystal_lattice.a_star_vector() +
             miller_indices.k() * crystal_lattice.b_star_vector() +
             miller_indices.l() * crystal_lattice.c_star_vector() );
}

// ********************************************************************************

NormalisedVector3D reciprocal_lattice_direction( const MillerIndices miller_indices, const CrystalLattice & crystal_lattice )
{
    return normalised_vector( miller_indices.h() * crystal_lattice.a_star_vector() +
                              miller_indices.k() * crystal_lattice.b_star_vector() +
                              miller_indices.l() * crystal_lattice.c_star_vector() );
}

// ********************************************************************************

// Gram-Schmidt orthogonalisation
NormalisedVector3D orthogonalise( const NormalisedVector3D & n, const Vector3D & r )
{
    return normalised_vector( r - (n*r) * n );
}

// ********************************************************************************

Vector3D change_of_basis( const Vector3D & point, const CoordinateFrame & before, const CoordinateFrame & after )
{
    return Vector3D( point.x() * before.x_axis() * after.x_axis() + point.y() * before.y_axis() * after.x_axis() + point.z() * before.z_axis() * after.x_axis(),
                     point.x() * before.x_axis() * after.y_axis() + point.y() * before.y_axis() * after.y_axis() + point.z() * before.z_axis() * after.y_axis(),
                     point.x() * before.x_axis() * after.z_axis() + point.y() * before.y_axis() * after.z_axis() + point.z() * before.z_axis() * after.z_axis() );
}

// ********************************************************************************

Plane plane( const CollectionOfPoints & points )
{
    // @@ This suggests that we must make CollectionOfPoints, renamed as SetOfPoints, a low-level class.
    Plane result( points.points_wrt_com() );
    return result;
}

// ********************************************************************************

Vector2D projection( const Plane & plane, const Vector3D & point )
{
    Vector3D projected_point = plane.projection3D( point );
    projected_point = change_of_basis( projected_point, CoordinateFrame(), plane.coordinate_frame() );
    return Vector2D( projected_point.x(), projected_point.y() );
}

// ********************************************************************************

double root_mean_square_devation_from_mean_plane( const CollectionOfPoints & points, const Plane & plane )
{
    double result( 0.0 );
    for ( size_t i(0); i != points.size(); ++i )
        result += square( plane.normal() * points.point_wrt_com( i ) );
    result /= static_cast<double>(points.size());
    return sqrt( result );
}

// ********************************************************************************

double root_mean_square_devation_from_mean_plane( const std::vector< Vector3D > & points, const Plane & plane )
{
    double result( 0.0 );
    for ( size_t i(0); i != points.size(); ++i )
        result += square( plane.constant() - plane.normal() * points[i] );
    result /= static_cast<double>(points.size());
    return sqrt( result );
}

// ********************************************************************************

SymmetricMatrix3D covariance_matrix( const std::vector< Vector3D > & points )
{
    SymmetricMatrix3D result( 0.0 );
    Vector3D centre_of_mass = average( points );
    // Calculate covariance matrix
    for ( size_t iPoint( 0 ); iPoint != points.size(); ++iPoint )
    {
        Vector3D difference = points[iPoint] - centre_of_mass;
        for ( size_t i( 0 ); i != 3; ++i )
        {
            for ( size_t j( i ); j != 3; ++j )
            {
                result.set_value( i, j, result.value( i, j ) + difference.value( i ) * difference.value( j ) );
            }
        }
    }
    return ( 1.0 / static_cast<double>( points.size() ) ) * result;
}

// ********************************************************************************

SymmetricMatrix3D covariance_matrix( const CollectionOfPoints & points )
{
    SymmetricMatrix3D result( 0.0 );
    // Calculate covariance matrix
    for ( size_t iPoint( 0 ); iPoint != points.size(); ++iPoint )
    {
        for ( size_t i( 0 ); i != 3; ++i )
        {
            for ( size_t j( i ); j != 3; ++j )
            {
                result.set_value( i, j, result.value( i, j ) + points.point_wrt_com( iPoint ).value( i ) * points.point_wrt_com( iPoint ).value( j ) );
            }
        }
    }
    return ( 1.0 / static_cast<double>( points.size() ) ) * result;
}

// ********************************************************************************

NormalisedVector3D normalised_vector( const Vector3D & rhs )
{
    return NormalisedVector3D( rhs.x(), rhs.y(), rhs.z() );
}

// ********************************************************************************

Vector3D NormalisedVector3D2Vector3D( const NormalisedVector3D & rhs )
{
    return Vector3D( rhs.x(), rhs.y(), rhs.z() );
}

// ********************************************************************************

Vector3D operator+( const NormalisedVector3D & lhs, const NormalisedVector3D & rhs )
{
    return Vector3D( lhs.x() + rhs.x(), lhs.y() + rhs.y(), lhs.z() + rhs.z() );
}

// ********************************************************************************

Vector3D operator-( const NormalisedVector3D & lhs, const NormalisedVector3D & rhs )
{
    return Vector3D( lhs.x() - rhs.x(), lhs.y() - rhs.y(), lhs.z() - rhs.z() );
}

// ********************************************************************************

Vector3D operator-( const NormalisedVector3D & lhs, const Vector3D & rhs )
{
    return Vector3D( lhs.x() - rhs.x(), lhs.y() - rhs.y(), lhs.z() - rhs.z() );
}

// ********************************************************************************

Vector3D operator*( const double lhs, const NormalisedVector3D & rhs )
{
    return Vector3D( rhs.x()*lhs, rhs.y()*lhs, rhs.z()*lhs );
}

// ********************************************************************************

Vector3D operator*( const NormalisedVector3D & lhs, const double rhs )
{
    return Vector3D( lhs.x()*rhs, lhs.y()*rhs, lhs.z()*rhs );
}

// ********************************************************************************

Matrix3D SymmetricMatrix3D2Matrix3D( const SymmetricMatrix3D & matrix )
{
    return Matrix3D( matrix.value( 0, 0 ), matrix.value( 0, 1 ), matrix.value( 0, 2 ),
                     matrix.value( 1, 0 ), matrix.value( 1, 1 ), matrix.value( 1, 2 ),
                     matrix.value( 2, 0 ), matrix.value( 2, 1 ), matrix.value( 2, 2 )
                   );
}

// ********************************************************************************

SymmetricMatrix3D Matrix3D2SymmetricMatrix3D( const Matrix3D & matrix, const double tolerance )
{
    if ( nearly_equal( matrix.value( 0, 1 ), matrix.value( 1, 0 ), tolerance ) &&
         nearly_equal( matrix.value( 0, 2 ), matrix.value( 2, 0 ), tolerance ) &&
         nearly_equal( matrix.value( 1, 2 ), matrix.value( 2, 1 ), tolerance ) )
        return SymmetricMatrix3D( matrix.value( 0, 0 ), matrix.value( 1, 1 ), matrix.value( 2, 2 ),
                                  ( matrix.value( 0, 1 ) + matrix.value( 1, 0 ) ) / 2.0,
                                  ( matrix.value( 0, 2 ) + matrix.value( 2, 0 ) ) / 2.0,
                                  ( matrix.value( 1, 2 ) + matrix.value( 2, 1 ) ) / 2.0 );
    throw std::runtime_error( "Matrix3D2SymmetricMatrix3D() : input matrix is not symmetric." );
}

// ********************************************************************************

Vector3D operator*( const Matrix3D & matrix, const Vector3D & vector )
{
    return Vector3D( matrix.value( 0, 0 ) * vector.x() + matrix.value( 0, 1 ) * vector.y() + matrix.value( 0, 2 ) * vector.z(),
                     matrix.value( 1, 0 ) * vector.x() + matrix.value( 1, 1 ) * vector.y() + matrix.value( 1, 2 ) * vector.z(),
                     matrix.value( 2, 0 ) * vector.x() + matrix.value( 2, 1 ) * vector.y() + matrix.value( 2, 2 ) * vector.z()
                   );
}

// ********************************************************************************

Vector3D operator*( const SymmetricMatrix3D & matrix, const Vector3D & vector )
{
    return Vector3D( matrix.value( 0, 0 ) * vector.x() + matrix.value( 0, 1 ) * vector.y() + matrix.value( 0, 2 ) * vector.z(),
                     matrix.value( 1, 0 ) * vector.x() + matrix.value( 1, 1 ) * vector.y() + matrix.value( 1, 2 ) * vector.z(),
                     matrix.value( 2, 0 ) * vector.x() + matrix.value( 2, 1 ) * vector.y() + matrix.value( 2, 2 ) * vector.z()
                   );
}

// ********************************************************************************

Vector3D operator*( const Vector3D & vector, const Matrix3D & matrix )
{
    return Vector3D( matrix.value( 0, 0 ) * vector.x() + matrix.value( 1, 0 ) * vector.y() + matrix.value( 2, 0 ) * vector.z(),
                     matrix.value( 0, 1 ) * vector.x() + matrix.value( 1, 1 ) * vector.y() + matrix.value( 2, 1 ) * vector.z(),
                     matrix.value( 0, 2 ) * vector.x() + matrix.value( 1, 2 ) * vector.y() + matrix.value( 2, 2 ) * vector.z()
                   );
}

// ********************************************************************************

Vector3D operator*( const Vector3D & vector, const SymmetricMatrix3D & matrix )
{
    return vector * SymmetricMatrix3D2Matrix3D( matrix );
//    return Vector3D( matrix.value( 0, 0 ) * vector.x() + matrix.value( 1, 0 ) * vector.y() + matrix.value( 2, 0 ) * vector.z(),
//                     matrix.value( 0, 1 ) * vector.x() + matrix.value( 1, 1 ) * vector.y() + matrix.value( 2, 1 ) * vector.z(),
//                     matrix.value( 0, 2 ) * vector.x() + matrix.value( 1, 2 ) * vector.y() + matrix.value( 2, 2 ) * vector.z()
//                   );
}

// ********************************************************************************

Vector3D operator*( const NormalisedVector3D & vector, const Matrix3D & matrix )
{
    return NormalisedVector3D2Vector3D( vector ) * matrix;
//    return Vector3D( matrix.value( 0, 0 ) * vector.x() + matrix.value( 1, 0 ) * vector.y() + matrix.value( 2, 0 ) * vector.z(),
//                     matrix.value( 0, 1 ) * vector.x() + matrix.value( 1, 1 ) * vector.y() + matrix.value( 2, 1 ) * vector.z(),
//                     matrix.value( 0, 2 ) * vector.x() + matrix.value( 1, 2 ) * vector.y() + matrix.value( 2, 2 ) * vector.z()
//                   );
}

// ********************************************************************************

Vector3D operator*( const NormalisedVector3D & vector, const SymmetricMatrix3D & matrix )
{
    return NormalisedVector3D2Vector3D( vector ) * SymmetricMatrix3D2Matrix3D( matrix );
//    return Vector3D( matrix.value( 0, 0 ) * vector.x() + matrix.value( 1, 0 ) * vector.y() + matrix.value( 2, 0 ) * vector.z(),
//                     matrix.value( 0, 1 ) * vector.x() + matrix.value( 1, 1 ) * vector.y() + matrix.value( 2, 1 ) * vector.z(),
//                     matrix.value( 0, 2 ) * vector.x() + matrix.value( 1, 2 ) * vector.y() + matrix.value( 2, 2 ) * vector.z()
//                   );
}

// ********************************************************************************

Matrix3D operator*( const SymmetricMatrix3D & lhs, const SymmetricMatrix3D & rhs )
{
    Matrix3D result;
    result.set_value( 0, 0, lhs.value(0,0) * rhs.value(0,0) + lhs.value(0,1) * rhs.value(1,0) + lhs.value(0,2) * rhs.value(2,0) );
    result.set_value( 0, 1, lhs.value(0,0) * rhs.value(0,1) + lhs.value(0,1) * rhs.value(1,1) + lhs.value(0,2) * rhs.value(2,1) );
    result.set_value( 0, 2, lhs.value(0,0) * rhs.value(0,2) + lhs.value(0,1) * rhs.value(1,2) + lhs.value(0,2) * rhs.value(2,2) );
    result.set_value( 1, 0, lhs.value(1,0) * rhs.value(0,0) + lhs.value(1,1) * rhs.value(1,0) + lhs.value(1,2) * rhs.value(2,0) );
    result.set_value( 1, 1, lhs.value(1,0) * rhs.value(0,1) + lhs.value(1,1) * rhs.value(1,1) + lhs.value(1,2) * rhs.value(2,1) );
    result.set_value( 1, 2, lhs.value(1,0) * rhs.value(0,2) + lhs.value(1,1) * rhs.value(1,2) + lhs.value(1,2) * rhs.value(2,2) );
    result.set_value( 2, 0, lhs.value(2,0) * rhs.value(0,0) + lhs.value(2,1) * rhs.value(1,0) + lhs.value(2,2) * rhs.value(2,0) );
    result.set_value( 2, 1, lhs.value(2,0) * rhs.value(0,1) + lhs.value(2,1) * rhs.value(1,1) + lhs.value(2,2) * rhs.value(2,1) );
    result.set_value( 2, 2, lhs.value(2,0) * rhs.value(0,2) + lhs.value(2,1) * rhs.value(1,2) + lhs.value(2,2) * rhs.value(2,2) );
    return result;
}

// ********************************************************************************

Matrix3D operator*( const Matrix3D & lhs, const SymmetricMatrix3D & rhs )
{
    Matrix3D result;
    result.set_value( 0, 0, lhs.value(0,0) * rhs.value(0,0) + lhs.value(0,1) * rhs.value(1,0) + lhs.value(0,2) * rhs.value(2,0) );
    result.set_value( 0, 1, lhs.value(0,0) * rhs.value(0,1) + lhs.value(0,1) * rhs.value(1,1) + lhs.value(0,2) * rhs.value(2,1) );
    result.set_value( 0, 2, lhs.value(0,0) * rhs.value(0,2) + lhs.value(0,1) * rhs.value(1,2) + lhs.value(0,2) * rhs.value(2,2) );
    result.set_value( 1, 0, lhs.value(1,0) * rhs.value(0,0) + lhs.value(1,1) * rhs.value(1,0) + lhs.value(1,2) * rhs.value(2,0) );
    result.set_value( 1, 1, lhs.value(1,0) * rhs.value(0,1) + lhs.value(1,1) * rhs.value(1,1) + lhs.value(1,2) * rhs.value(2,1) );
    result.set_value( 1, 2, lhs.value(1,0) * rhs.value(0,2) + lhs.value(1,1) * rhs.value(1,2) + lhs.value(1,2) * rhs.value(2,2) );
    result.set_value( 2, 0, lhs.value(2,0) * rhs.value(0,0) + lhs.value(2,1) * rhs.value(1,0) + lhs.value(2,2) * rhs.value(2,0) );
    result.set_value( 2, 1, lhs.value(2,0) * rhs.value(0,1) + lhs.value(2,1) * rhs.value(1,1) + lhs.value(2,2) * rhs.value(2,1) );
    result.set_value( 2, 2, lhs.value(2,0) * rhs.value(0,2) + lhs.value(2,1) * rhs.value(1,2) + lhs.value(2,2) * rhs.value(2,2) );
    return result;
}

// ********************************************************************************

Matrix3D operator*( const SymmetricMatrix3D & lhs, const Matrix3D & rhs )
{
    Matrix3D result;
    result.set_value( 0, 0, lhs.value(0,0) * rhs.value(0,0) + lhs.value(0,1) * rhs.value(1,0) + lhs.value(0,2) * rhs.value(2,0) );
    result.set_value( 0, 1, lhs.value(0,0) * rhs.value(0,1) + lhs.value(0,1) * rhs.value(1,1) + lhs.value(0,2) * rhs.value(2,1) );
    result.set_value( 0, 2, lhs.value(0,0) * rhs.value(0,2) + lhs.value(0,1) * rhs.value(1,2) + lhs.value(0,2) * rhs.value(2,2) );
    result.set_value( 1, 0, lhs.value(1,0) * rhs.value(0,0) + lhs.value(1,1) * rhs.value(1,0) + lhs.value(1,2) * rhs.value(2,0) );
    result.set_value( 1, 1, lhs.value(1,0) * rhs.value(0,1) + lhs.value(1,1) * rhs.value(1,1) + lhs.value(1,2) * rhs.value(2,1) );
    result.set_value( 1, 2, lhs.value(1,0) * rhs.value(0,2) + lhs.value(1,1) * rhs.value(1,2) + lhs.value(1,2) * rhs.value(2,2) );
    result.set_value( 2, 0, lhs.value(2,0) * rhs.value(0,0) + lhs.value(2,1) * rhs.value(1,0) + lhs.value(2,2) * rhs.value(2,0) );
    result.set_value( 2, 1, lhs.value(2,0) * rhs.value(0,1) + lhs.value(2,1) * rhs.value(1,1) + lhs.value(2,2) * rhs.value(2,1) );
    result.set_value( 2, 2, lhs.value(2,0) * rhs.value(0,2) + lhs.value(2,1) * rhs.value(1,2) + lhs.value(2,2) * rhs.value(2,2) );
    return result;
}

// ********************************************************************************

Matrix3D operator+( const Matrix3D & lhs, const SymmetricMatrix3D & rhs )
{
    Matrix3D result;
    result.set_value( 0, 0, lhs.value(0,0) + rhs.value(0,0) );
    result.set_value( 0, 1, lhs.value(0,1) + rhs.value(0,1) );
    result.set_value( 0, 2, lhs.value(0,2) + rhs.value(0,2) );
    result.set_value( 1, 0, lhs.value(1,0) + rhs.value(1,0) );
    result.set_value( 1, 1, lhs.value(1,1) + rhs.value(1,1) );
    result.set_value( 1, 2, lhs.value(1,2) + rhs.value(1,2) );
    result.set_value( 2, 0, lhs.value(2,0) + rhs.value(2,0) );
    result.set_value( 2, 1, lhs.value(2,1) + rhs.value(2,1) );
    result.set_value( 2, 2, lhs.value(2,2) + rhs.value(2,2) );
    return result;
}

// ********************************************************************************

Matrix3D operator+( const SymmetricMatrix3D & lhs, const Matrix3D & rhs )
{
    return ( rhs + lhs );
}

// ********************************************************************************

Matrix3D operator-( const Matrix3D & lhs, const SymmetricMatrix3D & rhs )
{
    return ( lhs + ( -1.0 * rhs ) );
}

// ********************************************************************************

Matrix3D operator-( const SymmetricMatrix3D & lhs, const Matrix3D & rhs )
{
    return ( lhs + ( -1.0 * rhs ) );
}

// ********************************************************************************

// The transposition is implied
double operator*( const NormalisedVector3D & lhs, const Vector3D & rhs )
{
    return ( lhs.x() * rhs.x() + lhs.y() * rhs.y() + lhs.z() * rhs.z() );
}

// ********************************************************************************

double operator*( const Vector3D & lhs, const NormalisedVector3D & rhs )
{
    return ( lhs.x() * rhs.x() + lhs.y() * rhs.y() + lhs.z() * rhs.z() );
}

// ********************************************************************************

//MillerIndices operator*( const Matrix3D & matrix, const MillerIndices & Miller_indices )
//{
//    return MillerIndices( round_to_int( matrix.value( 0, 0 ) * Miller_indices.h() + matrix.value( 0, 1 ) * Miller_indices.k() + matrix.value( 0, 2 ) * Miller_indices.l() ),
//                          round_to_int( matrix.value( 1, 0 ) * Miller_indices.h() + matrix.value( 1, 1 ) * Miller_indices.k() + matrix.value( 1, 2 ) * Miller_indices.l() ),
//                          round_to_int( matrix.value( 2, 0 ) * Miller_indices.h() + matrix.value( 2, 1 ) * Miller_indices.k() + matrix.value( 2, 2 ) * Miller_indices.l() )
//                        );
//}

// ********************************************************************************

MillerIndices operator*( const MillerIndices & Miller_indices, const Matrix3D & matrix )
{
    return MillerIndices( round_to_int( Miller_indices.h() * matrix.value( 0, 0 ) + Miller_indices.k() * matrix.value( 1, 0 ) + Miller_indices.l() * matrix.value( 2, 0 ) ),
                          round_to_int( Miller_indices.h() * matrix.value( 0, 1 ) + Miller_indices.k() * matrix.value( 1, 1 ) + Miller_indices.l() * matrix.value( 2, 1 ) ),
                          round_to_int( Miller_indices.h() * matrix.value( 0, 2 ) + Miller_indices.k() * matrix.value( 1, 2 ) + Miller_indices.l() * matrix.value( 2, 2 ) )
                        );
}

// ********************************************************************************

MillerIndices operator*( const MillerIndices & Miller_indices, const SymmetricMatrix3D & matrix )
{
    return MillerIndices( round_to_int( Miller_indices.h() * matrix.value( 0, 0 ) + Miller_indices.k() * matrix.value( 1, 0 ) + Miller_indices.l() * matrix.value( 2, 0 ) ),
                          round_to_int( Miller_indices.h() * matrix.value( 0, 1 ) + Miller_indices.k() * matrix.value( 1, 1 ) + Miller_indices.l() * matrix.value( 2, 1 ) ),
                          round_to_int( Miller_indices.h() * matrix.value( 0, 2 ) + Miller_indices.k() * matrix.value( 1, 2 ) + Miller_indices.l() * matrix.value( 2, 2 ) )
                        );
}

// ********************************************************************************

double operator*( const MillerIndices & miller_indices, const Vector3D & vector_3D )
{
    return miller_indices.h() * vector_3D.x() + miller_indices.k() * vector_3D.y() + miller_indices.l() * vector_3D.z();
}

// ********************************************************************************

Matrix3D rotation_about_x( const Angle angle )
{
    return Matrix3D( 1.0, 0.0           , 0.0           ,
                     0.0, angle.cosine(), -angle.sine() ,
                     0.0, angle.sine()  , angle.cosine() );
}

// ********************************************************************************

Matrix3D rotation_about_y( const Angle angle )
{
    return Matrix3D( angle.cosine(), 0.0, angle.sine()  ,
                     0.0           , 1.0, 0.0           ,
                     -angle.sine() , 0.0, angle.cosine() );
}

// ********************************************************************************

Matrix3D rotation_about_z( const Angle angle )
{
    return Matrix3D( angle.cosine(), -angle.sine() ,  0.0,
                     angle.sine()  , angle.cosine(),  0.0,
                     0.0           , 0.0           ,  1.0 );
}

// ********************************************************************************

Vector3D rotate_point_about_axis( Vector3D point, const Vector3D & origin, const NormalisedVector3D & n, const Angle angle )
{
    CoordinateFrame coordination_frame( n );
    NormalisedVector3D basis_vector_2 = coordination_frame.y_axis();
    NormalisedVector3D basis_vector_3 = coordination_frame.z_axis();
    // Calculate the projection of "point" onto the axis (i.e. the point on the axis closest to "point".
    double omega = (n*point)-(n*origin);
    Vector3D new_origin = origin + (omega * n);
    point -= new_origin;
    // Measure current angle with respect to basis_vector_2.
    double x = point * basis_vector_2;
    double y = point * basis_vector_3;
    Angle current_angle = ATAN2( y, x );
    point = ( point.length() * (current_angle+angle).cosine() * basis_vector_2 ) + ( point.length() * (current_angle+angle).sine() * basis_vector_3 );
    point += new_origin;
    return point;
}

// ********************************************************************************

// Converts the orientation of a vector to Eulerian angles alpha and beta.
// Eulerian angles have many problems with definitions of origin and domain and ambiguous values.
void Eulerian_angles( const Vector3D & rhs, Angle & alpha, Angle & beta )
{
    // beta is the angle between the vector and the z-axis, domain 0.0 <= beta <= 180.0.
    if ( nearly_zero( rhs.x() ) && nearly_zero( rhs.y() ) )
    {
        if ( nearly_zero( rhs.z() ) )
            throw std::runtime_error( "Eulerian_angles( Vector3D, Angle, Angle ): error: the zero vector does not have an orientation." );
        alpha = Angle();
        if ( rhs.z() > 0.0 )
            beta = Angle();
        else
            beta = Angle::angle_180_degrees();
        return;
    }
    beta = arccosine( rhs.z() / rhs.length() );
    // alpha is the orientation in the xy-plane, the phase is measured w.r.t. the x-axis, i.e. the x-axis is the origin.
    alpha = ATAN2( rhs.y(), rhs.x() );
}

// ********************************************************************************

// Converts the orientation of a vector to Eulerian angles alpha and beta.
// Eulerian angles have many problems with definitions of origin and domain and ambiguous values.
void Eulerian_angles( const NormalisedVector3D & rhs, Angle & alpha, Angle & beta )
{
    // beta is the angle between the vector and the z-axis, domain 0.0 <= beta <= 180.0.
    if ( nearly_zero( rhs.x() ) && nearly_zero( rhs.y() ) )
    {
        alpha = Angle();
        if ( rhs.z() > 0.0 )
            beta = Angle();
        else
            beta = Angle::angle_180_degrees();
        return;
    }
    beta = arccosine( rhs.z() );
    // alpha is the orientation in the xy-plane, the phase is measured w.r.t. the x-axis, i.e. the x-axis is the origin.
    alpha = ATAN2( rhs.y(), rhs.x() );
}

// ********************************************************************************

Vector3D cylindrical2Cartesian( const double r, Angle phi, const double z )
{
    return Vector3D( r * phi.cosine(), r * phi.sine(), z );
}

// ********************************************************************************

// Takes two fractional coordinates and determines if they are the same, taking into acccount translations
bool are_translationally_equivalent( const double x, const double y  )
{
    double tolerance = 0.0001;
    double distance = std::abs( x - y );
    double dummy;
    distance = modf( distance, &dummy ); // modf() returns the fractional part
    // One special case left: if distance is now 0.9999, then it should be 0.0001
    return ( ( ( distance - tolerance ) < 0.0 ) ||
             ( ( distance + tolerance ) > 1.0 ) );
}

// ********************************************************************************

bool are_translationally_equivalent( const Vector3D & lhs, const Vector3D & rhs  )
{
    return ( are_translationally_equivalent( lhs.x(), rhs.x() ) &&
             are_translationally_equivalent( lhs.y(), rhs.y() ) &&
             are_translationally_equivalent( lhs.z(), rhs.z() ) );
}

// ********************************************************************************

// Atoms with coordinates like 0.999999: keep at 0.999999 or move to 0.0 or move to -0.000001?
double adjust_for_translations( const double input )
{
    double integer_part;
    double result = modf( input, &integer_part );
    if ( result < 0.0 )
        result += 1.0;
    if ( nearly_zero( result ) || nearly_equal( result, 1.0 ) )
        result = 0.0;
    return result;
}

// ********************************************************************************

// Atoms with coordinates like 0.999999: keep at 0.999999 or move to 0.0 or move to -0.000001?
// This should probably change in place,
// but our Vector3D class does not allow for its elements to be addressed that way
Vector3D adjust_for_translations( const Vector3D & input )
{
    Vector3D result;
    for ( size_t i( 0 ); i != 3; ++i )
        result.set_value( i, adjust_for_translations( input.value( i ) ) );
    return result;
}

// ********************************************************************************

// Returns the (smaller) angle between two vectors.
// Throws if at least one of the vectors is the zero vector.
Angle angle( const Vector3D & lhs, const Vector3D & rhs )
{
    if ( nearly_equal( ( lhs * rhs ) / ( lhs.length() * rhs.length() ), 1.0 ) )
        return Angle::from_degrees( 0.0 );
    if ( nearly_equal( ( lhs * rhs ) / ( lhs.length() * rhs.length() ), -1.0 ) )
        return Angle::from_degrees( 180.0 );
    return arccosine( ( lhs * rhs ) / ( lhs.length() * rhs.length() ) );
}

// ********************************************************************************

// Returns the (smaller) angle between two vectors.
Angle angle( const NormalisedVector3D & lhs, const Vector3D & rhs )
{
    if ( nearly_equal( ( lhs * rhs ) / rhs.length(), 1.0 ) )
        return Angle::from_degrees( 0.0 );
    if ( nearly_equal( ( lhs * rhs ) / rhs.length(), -1.0 ) )
        return Angle::from_degrees( 180.0 );
    return arccosine( ( lhs * rhs ) / rhs.length() );
}

// ********************************************************************************

Angle angle( const NormalisedVector3D & lhs, const Vector3D & rhs, const double length )
{
    if ( nearly_equal( ( lhs * rhs ) / length, 1.0 ) )
        return Angle::from_degrees( 0.0 );
    if ( nearly_equal( ( lhs * rhs ) / length, -1.0 ) )
        return Angle::from_degrees( 180.0 );
    return arccosine( ( lhs * rhs ) / length );
}

// ********************************************************************************

// Returns the (smaller) angle between two vectors.
Angle angle( const NormalisedVector3D & lhs, const NormalisedVector3D & rhs )
{
    if ( nearly_equal( lhs * rhs, 1.0 ) )
        return Angle::from_degrees( 0.0 );
    if ( nearly_equal( lhs * rhs, -1.0 ) )
        return Angle::from_degrees( 180.0 );
    return arccosine( lhs * rhs );
}

// ********************************************************************************

Angle angle( const Plane & lhs, const Plane & rhs )
{
    return arccosine( std::abs( lhs.normal() * rhs.normal() ) );
}

// ********************************************************************************

// The order of the points determines the sign of the torsion.
Angle signed_torsion( const Vector3D & r1, const Vector3D & r2, const Vector3D & r3, const Vector3D & r4 )
{
    Plane plane1( r1, r2, r3 );
    Plane plane2( r2, r3, r4 );
    Angle torsion = angle( plane1.normal(), plane2.normal() );
    double d = plane1.signed_distance( r4 );
    double s;
    if ( nearly_zero( d ) ) // Our sign() function can return 0...
        s = 1.0;
    else
        s = -sign( d ); // The minus is to ensure we get the same results as Mercury
    return s * torsion;
}

// ********************************************************************************

Matrix3D centred_to_primitive( const Centring & centring )
{
    if ( centring.is_primitive() )
    {
        std::cout << "centred_to_primitive( Centring ): warning: centring is primitive." << std::endl;
        return Matrix3D();
    }
    if ( centring.centring_type() == Centring::A )
        return Matrix3D(  1.0,  0.0,  0.0,
                          0.0,  0.5,  0.5,
                          0.0, -0.5,  0.5 );
    if ( centring.centring_type() == Centring::B )
        return Matrix3D(  0.5,  0.0,  0.5,
                          0.0,  1.0,  0.0,
                         -0.5,  0.0,  0.5 );
    if ( centring.centring_type() == Centring::C )
        return Matrix3D(  0.5,  0.5,  0.0,
                         -0.5,  0.5,  0.0,
                          0.0,  0.0,  1.0 );
    if ( centring.centring_type() == Centring::D )
        throw std::runtime_error( "centred_to_primitive( Centring ): centring D not yet implemented." );
    if ( centring.centring_type() == Centring::F )
        return Matrix3D(  0.0,  0.5,  0.5,
                          0.5,  0.0,  0.5,
                          0.5,  0.5,  0.0 );
    if ( centring.centring_type() == Centring::I )
        return Matrix3D( -0.5,  0.5,  0.5,
                          0.5, -0.5,  0.5,
                          0.5,  0.5, -0.5 );
    if ( centring.centring_type() == Centring::R_OBVERSE )
        return Matrix3D( 2.0/3.0, 1.0/3.0, 1.0/3.0,
                         1.0/3.0, 2.0/3.0, 2.0/3.0,
                           0.0,     0.0,     1.0 );
    if ( centring.centring_type() == Centring::R_REVERSE )
        return Matrix3D( 1.0/3.0, 2.0/3.0, 2.0/3.0,
                         2.0/3.0, 1.0/3.0, 1.0/3.0,
                           0.0,     0.0,     1.0 );
    if ( centring.centring_type() == Centring::U )
        throw std::runtime_error( "centred_to_primitive( Centring ): error: no transformation matrix for centring U." );
    if ( centring.centring_type() == Centring::J )
        return Matrix3D(  0.1,  0.0,  0.0,
                          0.0,  0.1,  0.0,
                          0.0,  0.0,  0.1 );
    throw std::runtime_error( "centred_to_primitive( Centring ): centring not yet implemented." );
}

// ********************************************************************************

bool nearly_equal( const CrystalStructure & lhs, const CrystalStructure & rhs )
{
    return true;
}

// ********************************************************************************

void add_centring_to_space_group_after_transformation( Matrix3D tranformation_matrix, SpaceGroup & space_group )
{
    double d = tranformation_matrix.determinant();
    if ( ! nearly_integer( d ) )
        throw std::runtime_error( "add_centring_to_space_group_after_transformation() : determinant is not an integer." );
    if ( nearly_zero( d ) )
        throw std::runtime_error( "add_centring_to_space_group_after_transformation() : determinant is zero." );
    if ( d < 0.0 )
        throw std::runtime_error( "add_centring_to_space_group_after_transformation() : determinant is negative." );
    size_t D = round_to_int( tranformation_matrix.determinant() );
    if ( D == 1 )
        return;
    // I have not been able to find a smart way to extract the possible additional lattice points
    // from the transformation matrix, so we simply try them all.
    // We want to find the points [ f/D, g/D, h/D ], with D the determinant, that lie within the unit cell
    // but that is not one of the current lattice points. So f/D, g/D and h/D are not allowed all to be integers at once.
    // f/D must be in the range [ 0, 1 >. 
    tranformation_matrix.transpose();
    // tranformation_matrix /= d;
    std::vector< Vector3D > centring_vectors;
    for ( size_t f( 0 ); f != D; ++f )
    {
        for ( size_t g( 0 ); g != D; ++g )
        {
            for ( size_t h( 0 ); h != D; ++h )
            {
                Vector3D trial_vector( f/d, g/d, h/d );
                // I guess it would be more efficient to divide the transformation matrix by d and
                // to construct the trial vector as [ f, g, h ].
                Vector3D lp = tranformation_matrix * trial_vector; // lp = (trial) lattice point in the old coordinate frame.
                if ( ! nearly_integer( lp.x() ) )
                    continue;
                if ( ! nearly_integer( lp.y() ) )
                    continue;
                if ( ! nearly_integer( lp.z() ) )
                    continue;
                centring_vectors.push_back( trial_vector );
            }
        }
    }
    if ( centring_vectors.size() != D )
        throw std::runtime_error( "add_centring_to_space_group_after_transformation() : centring_vectors.size() != D." );
    space_group.add_centring( Centring( centring_vectors ) );
}

// ********************************************************************************

