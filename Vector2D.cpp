/* *********************************************
Copyright (c) 2013-2025, Cornelis Jan (Jacco) van de Streek
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

#include "Vector2D.h"
#include "BasicMathsFunctions.h"
#include "Utilities.h"

#include <cmath>
#include <stdexcept>
#include <iostream>

// ********************************************************************************

Vector2D::Vector2D()
{
    data_[0] = 0.0;
    data_[1] = 0.0;
}

// ********************************************************************************

Vector2D::Vector2D( const double x, const double y )
{
    data_[0] = x;
    data_[1] = y;
}

// ********************************************************************************

double Vector2D::value( const size_t i ) const
{
    if ( i < 2 )
        return data_[i];
    throw std::runtime_error( "Vector2D::value(): index out of bounds." );
}

// ********************************************************************************

void Vector2D::set_value( const size_t i, const double value )
{
    if ( i < 2 )
        data_[i] = value;
    else
        throw std::runtime_error( "Vector2D::set_value(): index out of bounds." );
}

// ********************************************************************************

void Vector2D::set_length( const double value )
{
    throw_if_zero_vector();
    double l = value / length();
    data_[0] *= l;
    data_[1] *= l;
}

// ********************************************************************************

bool Vector2D::nearly_zero( const double tolerance ) const
{
    return ( ::nearly_zero( x(), tolerance ) && ::nearly_zero( y(), tolerance ) );
}

// ********************************************************************************

void Vector2D::throw_if_zero_vector( const double tolerance ) const
{
    if ( nearly_zero( tolerance ) )
        throw std::runtime_error( "Vector2D::throw_if_zero_vector()." );
}

// ********************************************************************************

Vector2D & Vector2D::operator+=( const Vector2D & rhs )
{
    *this = Vector2D( x() + rhs.x(), y() + rhs.y() );
    return *this;
}

// ********************************************************************************

Vector2D & Vector2D::operator-=( const Vector2D & rhs )
{
    *this = Vector2D( x() - rhs.x(), y() - rhs.y() );
    return *this;
}

// ********************************************************************************

Vector2D & Vector2D::operator/=( const double rhs )
{
    *this = Vector2D( x() / rhs, y() / rhs );
    return *this;
}

// ********************************************************************************

Vector2D & Vector2D::operator*=( const double rhs )
{
    *this = Vector2D( x() * rhs, y() * rhs );
    return *this;
}

// ********************************************************************************

Vector2D Vector2D::operator-() const
{
    return Vector2D( -data_[0], -data_[1] );
}

// ********************************************************************************

Vector2D Vector2D::operator+() const
{
    return Vector2D( data_[0], data_[1] );
}

// ********************************************************************************

double Vector2D::norm2() const
{
    return data_[0] * data_[0] + data_[1] * data_[1];
}

// ********************************************************************************

double Vector2D::length() const
{
    return sqrt( norm2() );
}

// ********************************************************************************

void Vector2D::orthogonalise( const Vector2D & other )
{
    double c = -(*this*other)/other.norm2();
    *this += c * other;
}

// ********************************************************************************

std::string Vector2D::to_string() const
{
   return double2string( x() ) + " " + double2string( y() );
}

// ********************************************************************************

void Vector2D::show() const
{
    std::cout << "x = " << x() << ", y = " << y() << std::endl;
}

// ********************************************************************************

std::string Vector2D::index2string( const size_t index )
{
    switch ( index )
    {
        case 0: return "x";
        case 1: return "y";
        default: throw std::runtime_error( "Vector2D::index2string( const size_t ): out of range." );
    }
}

// ********************************************************************************

std::ostream & operator<<( std::ostream & os, const Vector2D & vector2d )
{
    os << "[ " << vector2d.x() << ", " << vector2d.y() << " ]" << std::endl;
    return os;
}

// ********************************************************************************

bool nearly_equal( const Vector2D & lhs, const Vector2D & rhs, const double tolerance )
{
    return ( nearly_equal( lhs.x(), rhs.x(), tolerance ) && nearly_equal( lhs.y(), rhs.y(), tolerance ) );
}

// ********************************************************************************

bool nearly_zero( const Vector2D & lhs, const double tolerance )
{
    return lhs.nearly_zero( tolerance );
}

// ********************************************************************************

Vector2D operator+( const Vector2D & lhs, const Vector2D & rhs )
{
    return Vector2D( lhs.x() + rhs.x(), lhs.y() + rhs.y() );
}

// ********************************************************************************

Vector2D operator-( const Vector2D & lhs, const Vector2D & rhs )
{
    return Vector2D( lhs.x() - rhs.x(), lhs.y() - rhs.y() );
}

// ********************************************************************************

double operator*( const Vector2D & lhs, const Vector2D & rhs )
{
    return ( lhs.x() * rhs.x() + lhs.y() * rhs.y() );
}

// ********************************************************************************

Vector2D operator*( const Vector2D & lhs, const double rhs )
{
    return Vector2D( lhs.x() * rhs, lhs.y() * rhs );
}

// ********************************************************************************

Vector2D operator/( const Vector2D & lhs, const double rhs )
{
    return Vector2D( lhs.x() / rhs, lhs.y() / rhs );
}

// ********************************************************************************

Vector2D operator*( const double lhs, const Vector2D & rhs )
{
    return Vector2D( rhs.x() * lhs, rhs.y() * lhs );
}

// ********************************************************************************

//Vector2D cross_product( const Vector2D & lhs, const Vector2D & rhs )
//{
//    return Vector2D( lhs.y() * rhs.z() - lhs.z() * rhs.y(),
//                     lhs.z() * rhs.x() - lhs.x() * rhs.z() );
//}

// ********************************************************************************

Vector2D square( const Vector2D & vector2d )
{
    return Vector2D( vector2d.x() * vector2d.x(), vector2d.y() * vector2d.y() );
}

// ********************************************************************************

Vector2D sqrt( const Vector2D & vector2d )
{
    return Vector2D( std::sqrt( vector2d.x() ), std::sqrt( vector2d.y() ) );
}

// ********************************************************************************

