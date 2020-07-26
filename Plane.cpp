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

#include "Plane.h"
#include "Eigenvalue.h"
#include "SymmetricMatrix3D.h"
#include "Utilities.h"
#include "Vector3D.h"
#include "Vector3DCalculations.h"
#include "3DCalculations.h"

#include <cmath>
#include <iostream> // For debugging

// ********************************************************************************

Plane::Plane() :
c_(0.0)
{
    n_ = NormalisedVector3D();
}

// ********************************************************************************

Plane::Plane( const NormalisedVector3D & n, const double c ) :
n_(n),
c_(c)
{
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

