#ifndef VECTOR2D_H
#define VECTOR2D_H

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

#include "BasicMathsFunctions.h"

#include <cstddef> // For definition of size_t
#include <iosfwd>

/*
  A vector of dimension 3.
  
  There is no member function normalise(), instead use NormalisedVector2D n = normalised_vector( r );
*/
class Vector2D
{
public:

    // Default constructor: zero vector
    Vector2D();

    Vector2D( const double x, const double y );

    double x() const { return data_[0]; }
    double y() const { return data_[1]; }
    void set_x( const double value ) { data_[0] = value; }
    void set_y( const double value ) { data_[1] = value; }
    double value( const size_t i ) const;
    void set_value( const size_t i, const double value );

    void set_length( const double value );

    // This could either test if each of the three components is close to 0.0
    // or it could test if length() or norm2() is close to 0.0. I have decided to do the
    // first.
    bool nearly_zero( const double tolerance = TOLERANCE ) const;

    // This could either test if each of the three components is close to 0.0
    // or it could test if length() or norm2() is close to 0.0. I have decided to do the
    // first.
    void throw_if_zero_vector( const double tolerance = TOLERANCE ) const;

    Vector2D & operator+=( const Vector2D & rhs );
    Vector2D & operator-=( const Vector2D & rhs );
    Vector2D & operator/=( const double rhs );
    Vector2D & operator*=( const double rhs );

    Vector2D operator-() const;
    Vector2D operator+() const;

    double norm2() const;

    double length() const;

    void orthogonalise( const Vector2D & other );
    
    std::string to_string() const;

    void show() const; // For debugging
    
    static std::string index2string( const size_t index );

private:
    double data_[2];
};

std::ostream & operator<<( std::ostream & os, const Vector2D & vector3d );

bool nearly_equal( const Vector2D & lhs, const Vector2D & rhs, const double tolerance = TOLERANCE );
bool nearly_zero( const Vector2D & lhs, const double tolerance = TOLERANCE );

Vector2D operator+( const Vector2D & lhs, const Vector2D & rhs );
Vector2D operator-( const Vector2D & lhs, const Vector2D & rhs );
// The transposition is implied
double operator*( const Vector2D & lhs, const Vector2D & rhs );
Vector2D operator*( const Vector2D & lhs, const double rhs );
Vector2D operator/( const Vector2D & lhs, const double rhs );
Vector2D operator*( const double lhs, const Vector2D & rhs );

//Vector2D cross_product( const Vector2D & lhs, const Vector2D & rhs );

Vector2D square( const Vector2D & vector2d );
Vector2D sqrt( const Vector2D & vector2d );

#endif // VECTOR2D_H

