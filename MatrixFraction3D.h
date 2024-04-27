#ifndef MATRIXFRACTION3D_H
#define MATRIXFRACTION3D_H

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

#include "Fraction.h"

#include <cstddef> // For definition of size_t
#include <iosfwd>

/*
  A 3x3 matrix.
*/

class MatrixFraction3D
{
public:

    // Constructs a scalar matrix
    explicit MatrixFraction3D( const Fraction & scalar = Fraction( 1 ) );

    MatrixFraction3D( const Fraction & a00, const Fraction & a01, const Fraction & a02,
                      const Fraction & a10, const Fraction & a11, const Fraction & a12,
                      const Fraction & a20, const Fraction & a21, const Fraction & a22
                    );

    // In keeping with the silly C++ convention: zero-based.
    Fraction value( const size_t i, const size_t j ) const;

    void set_value( const size_t i, const size_t j, const Fraction & value );

    Fraction sum_of_elements() const;

    Fraction sum_of_absolute_elements() const;

    void invert();

    void transpose();

    Fraction determinant() const;
    
    Fraction trace() const;
    
    void swap_rows( const size_t i, const size_t j );

    void swap_columns( const size_t i, const size_t j );

    Fraction maximum_absolute_value_in_row( const size_t i ) const;

    Fraction maximum_absolute_value_in_column( const size_t i ) const;

    bool is_diagonal() const;

    // Returns the determinant of the minor matrix defined by element i,j.
    // In keeping with the silly C++ convention: zero-based.
    Fraction minor_matrix_determinant( const size_t i, const size_t j ) const;

    void convert_to_row_echelon_form( MatrixFraction3D & T );

    size_t number_of_zero_rows() const;
    
    bool is_the_identity() const;

    MatrixFraction3D & operator+=( const MatrixFraction3D & rhs );
    MatrixFraction3D & operator-=( const MatrixFraction3D & rhs );
    MatrixFraction3D & operator/=( const Fraction value );

    bool operator==( const MatrixFraction3D & rhs ) const;

    void show() const; // For debugging, should be something like #include <iosfwd> print( std::os & ) const;

private:
    Fraction data_[3][3];
};

std::ostream & operator<<( std::ostream & os, const MatrixFraction3D & matrix3d );

MatrixFraction3D operator+( const MatrixFraction3D & lhs, const MatrixFraction3D & rhs );
MatrixFraction3D operator-( const MatrixFraction3D & lhs, const MatrixFraction3D & rhs );

MatrixFraction3D operator*( const MatrixFraction3D & lhs, const MatrixFraction3D & rhs );

MatrixFraction3D operator*( const Fraction & lhs, const MatrixFraction3D & rhs );
MatrixFraction3D operator*( const MatrixFraction3D & lhs, const Fraction & rhs );

MatrixFraction3D operator/( const MatrixFraction3D & lhs, Fraction rhs );

MatrixFraction3D inverse( const MatrixFraction3D & matrix3d );
MatrixFraction3D transpose( const MatrixFraction3D & matrix3d );


#endif // MATRIXFRACTION3D_H

