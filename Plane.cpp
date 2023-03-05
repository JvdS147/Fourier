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

#include "Plane.h"
#include "3DCalculations.h"
#include "Eigenvalue.h"
#include "SymmetricMatrix3D.h"
#include "Utilities.h"
#include "Vector3D.h"
#include "Vector3DCalculations.h"

#include <cmath>
#include <iostream> // For debugging

// ********************************************************************************

Plane::Plane() :
c_(0.0)
{
    n_ = NormalisedVector3D();
    calculate_coordinate_frame();
}

// ********************************************************************************

Plane::Plane( const NormalisedVector3D & n, const double c ) :
n_(n),
c_(c)
{
    calculate_coordinate_frame();
}

// ********************************************************************************

Plane::Plane( const Vector3D & r1, const Vector3D & r2, const Vector3D & r3 )
{
    Vector3D r2r1 = r2-r1;
    r2r1.throw_if_zero_vector();
    Vector3D r3r1 = r3-r1;
    r3r1.throw_if_zero_vector();
    Vector3D normal = cross_product( r2r1, r3r1 );
    n_ = normalised_vector( normal );
    c_ = n_ * ( (r1 + r2 + r3 ) / 3.0 );
    calculate_coordinate_frame();
}

// ********************************************************************************

Plane::Plane( const std::vector< Vector3D > & points )
{
    SymmetricMatrix3D covariance_matrix2 = covariance_matrix( points );
    std::vector< double > eigenvalues;
    std::vector< NormalisedVector3D > eigenvectors;
    calculate_eigenvalues( covariance_matrix2, eigenvalues, eigenvectors );
    // The eigenvectors are sorted by eigenvalue, the smallest is the first
    n_ = eigenvectors[0];
    // The centre of mass of the points lies in the plane
    c_ = n_ * average( points );
    calculate_coordinate_frame();
}

// ********************************************************************************

double Plane::signed_distance( const Vector3D & r ) const
{
    return c_ - n_*r;
}

// ********************************************************************************

double Plane::distance( const Vector3D & r ) const
{
    return std::abs( signed_distance( r ) );
}

// ********************************************************************************

// Projection along the normal onto the plane, turning a 3D point into a 2D point.
Vector3D Plane::projection3D( const Vector3D & point ) const
{
    return Vector3D( point + signed_distance( point ) * n_ );
}

// ********************************************************************************

// Projection along the normal onto the plane, turning a 3D point into a 2D point.
// Maybe we need a Vector2D class.
//Vector2D Plane::projection2D( const Vector3D & point ) const
//{
//    
//}

// ********************************************************************************

std::string Plane::to_string() const
{
   return "n = " + n_.to_string() + " c = " + double2string( c_ );
}

// ********************************************************************************

void Plane::show() const
{
    std::cout << this->to_string() << std::endl;
}

// ********************************************************************************

void Plane::calculate_coordinate_frame()
{
    coordinate_frame_ = CoordinateFrame( n_ );
    coordinate_frame_ = CoordinateFrame( coordinate_frame_.y_axis(), coordinate_frame_.z_axis(), coordinate_frame_.x_axis() );
}

// ********************************************************************************

