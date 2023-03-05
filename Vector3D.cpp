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

#include "Vector3D.h"
#include "BasicMathsFunctions.h"
#include "Utilities.h"

#include <cmath>
#include <stdexcept>
#include <iostream>

// ********************************************************************************

Vector3D::Vector3D()
{
    data_[0] = 0.0;
    data_[1] = 0.0;
    data_[2] = 0.0;
}

// ********************************************************************************

Vector3D::Vector3D( const double x, const double y, const double z )
{
    data_[0] = x;
    data_[1] = y;
    data_[2] = z;
}

// ********************************************************************************

double Vector3D::value( const size_t i ) const
{
    if ( i < 3 )
        return data_[i];
    throw std::runtime_error( "Vector3D::value(): index out of bounds." );
}

// ********************************************************************************

void Vector3D::set_value( const size_t i, const double value )
{
    if ( i < 3 )
        data_[i] = value;
    else
        throw std::runtime_error( "Vector3D::set_value(): index out of bounds." );
}

// ********************************************************************************

void Vector3D::set_length( const double value )
{
    throw_if_zero_vector();
    double l = value / length();
    data_[0] *= l;
    data_[1] *= l;
    data_[2] *= l;
}

// ********************************************************************************

bool Vector3D::nearly_zero( const double tolerance ) const
{
    return ( ::nearly_zero( x(), tolerance ) && ::nearly_zero( y(), tolerance ) && ::nearly_zero( z(), tolerance ) );
}

// ********************************************************************************

void Vector3D::throw_if_zero_vector( const double tolerance ) const
{
    if ( nearly_zero( tolerance ) )
        throw std::runtime_error( "Vector3D::throw_if_zero_vector()." );
}

// ********************************************************************************

Vector3D & Vector3D::operator+=( const Vector3D & rhs )
{
    *this = Vector3D( x() + rhs.x(), y() + rhs.y(), z() + rhs.z() );
    return *this;
}

// ********************************************************************************

Vector3D & Vector3D::operator-=( const Vector3D & rhs )
{
    *this = Vector3D( x() - rhs.x(), y() - rhs.y(), z() - rhs.z() );
    return *this;
}

// ********************************************************************************

Vector3D & Vector3D::operator/=( const double rhs )
{
    *this = Vector3D( x() / rhs, y() / rhs, z() / rhs );
    return *this;
}

// ********************************************************************************

Vector3D & Vector3D::operator*=( const double rhs )
{
    *this = Vector3D( x() * rhs, y() * rhs, z() * rhs );
    return *this;
}

// ********************************************************************************

Vector3D Vector3D::operator-() const
{
    return Vector3D( -data_[0], -data_[1], -data_[2] );
}

// ********************************************************************************

Vector3D Vector3D::operator+() const
{
    return Vector3D( data_[0], data_[1], data_[2] );
}

// ********************************************************************************

double Vector3D::norm2() const
{
    return data_[0] * data_[0] + data_[1] * data_[1] + data_[2] * data_[2];
}

// ********************************************************************************

double Vector3D::length() const
{
    return sqrt( norm2() );
}

// ********************************************************************************

void Vector3D::orthogonalise( const Vector3D & other )
{
    double c = -(*this*other)/other.norm2();
    *this += c * other;
}

// ********************************************************************************

std::string Vector3D::to_string() const
{
   return double2string( x() ) + " " + double2string( y() ) + " " + double2string( z() );
}

// ********************************************************************************

void Vector3D::show() const
{
    std::cout << "x = " << x() << ", y = " << y() << ", z = " << z() << std::endl;
}

// ********************************************************************************

std::string Vector3D::index2string( const size_t index )
{
    switch ( index )
    {
        case 0: return "x";
        case 1: return "y";
        case 2: return "z";
        default: throw std::runtime_error( "Vector3D::index2string( const size_t ): out of range." );
    }
}

// ********************************************************************************

std::ostream & operator<<( std::ostream & os, const Vector3D & vector3d )
{
    os << "[ " << vector3d.x() << ", " << vector3d.y() << ", " << vector3d.z() << " ]" << std::endl;
    return os;
}

// ********************************************************************************

bool nearly_equal( const Vector3D & lhs, const Vector3D & rhs, const double tolerance )
{
    return ( nearly_equal( lhs.x(), rhs.x(), tolerance ) && nearly_equal( lhs.y(), rhs.y(), tolerance ) && nearly_equal( lhs.z(), rhs.z(), tolerance ) );
}

// ********************************************************************************

bool nearly_zero( const Vector3D & lhs, const double tolerance )
{
    return lhs.nearly_zero( tolerance );
}

// ********************************************************************************

Vector3D operator+( const Vector3D & lhs, const Vector3D & rhs )
{
    return Vector3D( lhs.x() + rhs.x(), lhs.y() + rhs.y(), lhs.z() + rhs.z() );
}

// ********************************************************************************

Vector3D operator-( const Vector3D & lhs, const Vector3D & rhs )
{
    return Vector3D( lhs.x() - rhs.x(), lhs.y() - rhs.y(), lhs.z() - rhs.z() );
}

// ********************************************************************************

double operator*( const Vector3D & lhs, const Vector3D & rhs )
{
    return ( lhs.x() * rhs.x() + lhs.y() * rhs.y() + lhs.z() * rhs.z() );
}

// ********************************************************************************

Vector3D operator*( const Vector3D & lhs, const double rhs )
{
    return Vector3D( lhs.x() * rhs, lhs.y() * rhs, lhs.z() * rhs );
}

// ********************************************************************************

Vector3D operator/( const Vector3D & lhs, const double rhs )
{
    return Vector3D( lhs.x() / rhs, lhs.y() / rhs, lhs.z() / rhs );
}

// ********************************************************************************

Vector3D operator*( const double lhs, const Vector3D & rhs )
{
    return Vector3D( rhs.x() * lhs, rhs.y() * lhs, rhs.z() * lhs );
}

// ********************************************************************************

Vector3D cross_product( const Vector3D & lhs, const Vector3D & rhs )
{
    return Vector3D( lhs.y() * rhs.z() - lhs.z() * rhs.y(),
                     lhs.z() * rhs.x() - lhs.x() * rhs.z(),
                     lhs.x() * rhs.y() - lhs.y() * rhs.x() );
}

// ********************************************************************************

Vector3D square( const Vector3D & vector3d )
{
    return Vector3D( vector3d.x() * vector3d.x(), vector3d.y() * vector3d.y(), vector3d.z() * vector3d.z() );
}

// ********************************************************************************

Vector3D sqrt( const Vector3D & vector3d )
{
    return Vector3D( std::sqrt( vector3d.x() ), std::sqrt( vector3d.y() ), std::sqrt( vector3d.z() ) );
}

// ********************************************************************************

