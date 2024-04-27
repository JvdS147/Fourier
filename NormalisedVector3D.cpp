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

#include "NormalisedVector3D.h"
#include "BasicMathsFunctions.h"
#include "Utilities.h"

#include <cmath>
#include <stdexcept>
#include <iostream> // For debugging

// ********************************************************************************
NormalisedVector3D::NormalisedVector3D()
{
    data_[0] = 1.0;
    data_[1] = 0.0;
    data_[2] = 0.0;
}
// ********************************************************************************

NormalisedVector3D::NormalisedVector3D( const double x, const double y, const double z )
{
    data_[0] = x;
    data_[1] = y;
    data_[2] = z;
    normalise_2();
}

// ********************************************************************************

NormalisedVector3D NormalisedVector3D::operator-() const
{
    return NormalisedVector3D( -data_[0], -data_[1], -data_[2] );
}

// ********************************************************************************

NormalisedVector3D NormalisedVector3D::operator+() const
{
    return NormalisedVector3D( data_[0], data_[1], data_[2] );
}

// ********************************************************************************

std::string NormalisedVector3D::to_string() const
{
   return double2string( x() ) + " " + double2string( y() ) + " " + double2string( z() );
}

// ********************************************************************************

void NormalisedVector3D::show() const
{
    std::cout << "x = " << x() << ", y = " << y() << ", z = " << z() << std::endl;
}

// ********************************************************************************

void NormalisedVector3D::normalise_2()
{
    double l = data_[0]*data_[0] + data_[1]*data_[1] + data_[2]*data_[2];
    if ( nearly_zero( l ) )
        throw std::runtime_error( "NormalisedVector3D::normalise_2(): zero vector" );
    l = sqrt( l );
    data_[0] /= l;
    data_[1] /= l;
    data_[2] /= l;
}

// ********************************************************************************

double operator*( const NormalisedVector3D& lhs, const NormalisedVector3D& rhs )
{
    return ( lhs.x() * rhs.x() + lhs.y() * rhs.y() + lhs.z() * rhs.z() );
}

// ********************************************************************************

NormalisedVector3D orthogonalise( const NormalisedVector3D & n, const NormalisedVector3D & r )
{
    // The expression is ( r - (n*r) * n ), but the intermediate value is a Vector3D, not a NormalisedVector3D,
    // which leads to dependencies between files.
    double t = n*r;
    return NormalisedVector3D( r.x() - t*n.x(), r.y() - t*n.y(), r.z() - t*n.z() );
}
// ********************************************************************************

