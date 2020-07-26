#ifndef MATRIX3D_H
#define MATRIX3D_H

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

#include <cstddef> // For definition of size_t
#include <iosfwd>

/*
  A 3x3 matrix.
*/
class Matrix3D
{
public:

    // Default constructor: identity matrix
    Matrix3D();

    // Constructs a scalar matrix
    explicit Matrix3D( const double scalar );

    Matrix3D( const double a00, const double a01, const double a02,
              const double a10, const double a11, const double a12,
              const double a20, const double a21, const double a22
            );

    // In keeping with the silly C++ convention: zero-based
    double value( const size_t i, const size_t j ) const { return data_[i][j]; }

    void set_value( const size_t i, const size_t j, const double value ) { data_[i][j] = value; }

    double sum_of_elements() const;
    double sum_of_absolute_elements() const;

    void invert();

    void transpose();

    double determinant() const;

    double trace() const;

    // Returns the determinant of the minor matrix defined by element i,j.
    // In keeping with the silly C++ convention: zero-based
    double minor_matrix_determinant( const size_t i, const size_t j ) const;

    Matrix3D & operator+=( const Matrix3D & rhs );
    Matrix3D & operator-=( const Matrix3D & rhs );
    Matrix3D & operator/=( const double value );

    bool operator==( const Matrix3D & rhs ) const;

    void show() const; // For debugging, should be something like #include <iosfwd> print( std::os & ) const;

private:
    double data_[3][3];
};

std::ostream & operator<<( std::ostream & os, const Matrix3D & matrix3d );

bool nearly_equal( const Matrix3D & lhs, const Matrix3D & rhs, const double tolerance = 0.0000001 );

Matrix3D operator+( const Matrix3D & lhs, const Matrix3D & rhs );
Matrix3D operator-( const Matrix3D & lhs, const Matrix3D & rhs );

Matrix3D operator*( const Matrix3D & lhs, const Matrix3D & rhs );

Matrix3D operator*( const double lhs, const Matrix3D & rhs );

Matrix3D inverse( const Matrix3D & matrix3d );
Matrix3D transpose( const Matrix3D & matrix3d );

#endif // MATRIX3D_H
