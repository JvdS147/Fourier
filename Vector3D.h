#ifndef VECTOR3D_H
#define VECTOR3D_H

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

#include <cstddef> // For definition of size_t
#include <iosfwd>

/*
  A vector of dimension 3.
  
  There is no member function normalise(), instead use NormalisedVector3D n = normalised_vector( r );
*/
class Vector3D
{
public:

    // Default constructor: zero vector
    Vector3D();

    Vector3D( const double x, const double y, const double z );

    double x() const { return data_[0]; }
    double y() const { return data_[1]; }
    double z() const { return data_[2]; }
    void set_x( const double value ) { data_[0] = value; }
    void set_y( const double value ) { data_[1] = value; }
    void set_z( const double value ) { data_[2] = value; }
    double value( const size_t i ) const { return data_[i]; }
    void set_value( const size_t i, const double value ) { data_[i] = value; }
    
    void set_length( const double value );

    // Checks if norm2() is smaller than tolerance
    // (To avoid the sqrt of length() .)
    bool is_zero_vector( const double tolerance = 0.000001 ) const;
    
    // Checks if norm2() is smaller than tolerance
    // (To avoid the sqrt of length() .)
    void throw_if_zero_vector( const double tolerance = 0.000001 ) const;

    Vector3D & operator+=( const Vector3D & rhs );
    Vector3D & operator-=( const Vector3D & rhs );
    Vector3D & operator/=( const double rhs );
    Vector3D & operator*=( const double rhs );

    Vector3D operator-() const;
    Vector3D operator+() const;

    double norm2() const;

    double length() const;

    void orthogonalise( const Vector3D & other );
    
    std::string to_string() const;

    void show() const; // For debugging
    
    static std::string index2string( const size_t index );

private:
    double data_[3];
};

std::ostream & operator<<( std::ostream & os, const Vector3D & vector3d );

bool nearly_equal( const Vector3D & lhs, const Vector3D & rhs, const double tolerance = 0.0000001 );

Vector3D operator+( const Vector3D & lhs, const Vector3D & rhs );
Vector3D operator-( const Vector3D & lhs, const Vector3D & rhs );
// The transposition is implied
double operator*( const Vector3D & lhs, const Vector3D & rhs );
Vector3D operator*( const Vector3D & lhs, const double rhs );
Vector3D operator/( const Vector3D & lhs, const double rhs );
Vector3D operator*( const double lhs, const Vector3D & rhs );

Vector3D cross_product( const Vector3D & lhs, const Vector3D & rhs );

Vector3D square( const Vector3D & vector3d );
Vector3D sqrt( const Vector3D & vector3d );

#endif // VECTOR3D_H

