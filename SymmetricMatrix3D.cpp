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

#include "SymmetricMatrix3D.h"
#include "BasicMathsFunctions.h"

#include <cmath>
#include <stdexcept>
#include <iostream>

/*  Layout of data_[6]:
       0,0               [0]
       1,0  1,1          [1]  [2]
       2,0  2,1  2,2     [3]  [4]  [5]
*/

// ********************************************************************************

SymmetricMatrix3D::SymmetricMatrix3D()
{
    *this = SymmetricMatrix3D( 1.0 );
}

// ********************************************************************************

// Constructs a scalar matrix
SymmetricMatrix3D::SymmetricMatrix3D( const double scalar )
{
    data_[0] = scalar;
    data_[1] = 0.0; data_[2] = scalar;
    data_[3] = 0.0; data_[4] = 0.0; data_[5] = scalar;
}

// ********************************************************************************

SymmetricMatrix3D::SymmetricMatrix3D( const double a00, const double a11, const double a22,
                                      const double a01, const double a02, const double a12 )
{
    data_[0] = a00;
    data_[1] = a01; data_[2] = a11;
    data_[3] = a02; data_[4] = a12; data_[5] = a22;
}

// ********************************************************************************

double SymmetricMatrix3D::value( size_t i, size_t j ) const
{
    if ( i < j )
        std::swap( i, j );
    if ( i > 2 )
        throw std::runtime_error( "SymmetricMatrix3D::value(): out of bounds." );
    return data_[ ((i*(i+1))/2) + j ];
}

// ********************************************************************************

void SymmetricMatrix3D::set_value( size_t i, size_t j, const double value )
{
    if ( i < j )
        std::swap( i, j );
    if ( i > 2 )
        throw std::runtime_error( "SymmetricMatrix3D::set_value(): out of bounds." );
    data_[ ((i*(i+1))/2) + j ] = value;
}

// ********************************************************************************

void SymmetricMatrix3D::invert()
{
    double D = determinant();
    if ( fabs( D ) < 0.000001 )
        throw std::runtime_error( "SymmetricMatrix3D::invert(): determinant = 0." );
    transpose();
    SymmetricMatrix3D adjoint;
    adjoint.data_[0] = +minor_matrix_determinant(0,0);
    adjoint.data_[1] = -minor_matrix_determinant(1,0); adjoint.data_[2] = +minor_matrix_determinant(1,1);
    adjoint.data_[3] = +minor_matrix_determinant(2,0); adjoint.data_[4] = -minor_matrix_determinant(2,1); adjoint.data_[5] = +minor_matrix_determinant(2,2);
    adjoint /= D;
    *this = adjoint;
}

// ********************************************************************************

double SymmetricMatrix3D::determinant() const
{
    return data_[0] * minor_matrix_determinant(0,0) -
           data_[1] * minor_matrix_determinant(0,1) +
           data_[3] * minor_matrix_determinant(0,2);
}

// ********************************************************************************

double SymmetricMatrix3D::trace() const
{
    return data_[0] + data_[2] + data_[5];
}

// ********************************************************************************

bool SymmetricMatrix3D::is_diagonal() const
{
    return ( nearly_zero( data_[1] ) && nearly_zero( data_[3] ) && nearly_zero( data_[4] ) );
}

// ********************************************************************************

double SymmetricMatrix3D::minor_matrix_determinant( const size_t row, const size_t col ) const
{
    double minor_matrix[2][2];
    size_t minor_row( 0 );
    size_t minor_col( 0 );
    for ( size_t i( 0 ); i != 3; ++i )
    {
        if ( i == row )
            continue;
        minor_col = 0;
        for ( size_t j( 0 ); j != 3; ++j )
        {
            if ( j == col )
                continue;
            minor_matrix[minor_row][minor_col] = value( i, j );
            ++minor_col;
        }
        ++minor_row;
    }
    return minor_matrix[0][0]*minor_matrix[1][1] - minor_matrix[0][1]*minor_matrix[1][0];
}

// ********************************************************************************

SymmetricMatrix3D & SymmetricMatrix3D::operator+=( const SymmetricMatrix3D & rhs )
{
    *this = SymmetricMatrix3D( value( 0, 0 ) + rhs.value( 0, 0 ), value( 1, 1 ) + rhs.value( 1, 1 ), value( 2, 2 ) + rhs.value( 2, 2 ),
                               value( 0, 1 ) + rhs.value( 0, 1 ), value( 0, 2 ) + rhs.value( 0, 2 ), value( 1, 2 ) + rhs.value( 1, 2 )
                             );
    return *this;
}

// ********************************************************************************

SymmetricMatrix3D & SymmetricMatrix3D::operator-=( const SymmetricMatrix3D & rhs )
{
    *this = SymmetricMatrix3D( value( 0, 0 ) - rhs.value( 0, 0 ), value( 1, 1 ) - rhs.value( 1, 1 ), value( 2, 2 ) - rhs.value( 2, 2 ),
                               value( 0, 1 ) - rhs.value( 0, 1 ), value( 0, 2 ) - rhs.value( 0, 2 ), value( 1, 2 ) - rhs.value( 1, 2 )
                             );
    return *this;
}

// ********************************************************************************

SymmetricMatrix3D & SymmetricMatrix3D::operator/=( const double value )
{
    data_[0] /= value;
    data_[1] /= value; data_[2] /= value;
    data_[3] /= value; data_[4] /= value; data_[5] /= value;
    return *this;
}

// ********************************************************************************

bool SymmetricMatrix3D::operator==( const SymmetricMatrix3D & rhs ) const
{
    return nearly_equal( *this, rhs );
}

// ********************************************************************************

void SymmetricMatrix3D::show() const
{
    std::cout << "| " << data_[0] << " " << data_[1] << " " << data_[3] << " |" << std::endl;
    std::cout << "| " << data_[1] << " " << data_[2] << " " << data_[4] << " |" << std::endl;
    std::cout << "| " << data_[3] << " " << data_[4] << " " << data_[5] << " |" << std::endl;
}

// ********************************************************************************

std::ostream & operator<<( std::ostream & os, const SymmetricMatrix3D & matrix3d )
{
    os << "| " << matrix3d.value( 0, 0 ) << " " << matrix3d.value( 0, 1 ) << " " << matrix3d.value( 0, 2 ) << " |" << std::endl;
    os << "| " << matrix3d.value( 1, 0 ) << " " << matrix3d.value( 1, 1 ) << " " << matrix3d.value( 1, 2 ) << " |" << std::endl;
    os << "| " << matrix3d.value( 2, 0 ) << " " << matrix3d.value( 2, 1 ) << " " << matrix3d.value( 2, 2 ) << " |" << std::endl;
    return os;
}

// ********************************************************************************

bool nearly_equal( const SymmetricMatrix3D & lhs, const SymmetricMatrix3D & rhs, const double tolerance )
{
    return ( nearly_equal( lhs.value( 0, 0 ), rhs.value( 0, 0 ), tolerance ) &&
             nearly_equal( lhs.value( 1, 0 ), rhs.value( 1, 0 ), tolerance ) && nearly_equal( lhs.value( 1, 1 ), rhs.value( 1, 1 ), tolerance ) &&
             nearly_equal( lhs.value( 2, 0 ), rhs.value( 2, 0 ), tolerance ) && nearly_equal( lhs.value( 2, 1 ), rhs.value( 2, 1 ), tolerance ) && nearly_equal( lhs.value( 2, 2 ), rhs.value( 2, 2 ), tolerance ) );
}

// ********************************************************************************

SymmetricMatrix3D operator+( const SymmetricMatrix3D & lhs, const SymmetricMatrix3D & rhs )
{
    return SymmetricMatrix3D( lhs.value( 0, 0 ) + rhs.value( 0, 0 ), lhs.value( 1, 1 ) + rhs.value( 1, 1 ), lhs.value( 2, 2 ) + rhs.value( 2, 2 ),
                              lhs.value( 0, 1 ) + rhs.value( 0, 1 ), lhs.value( 0, 2 ) + rhs.value( 0, 2 ), lhs.value( 1, 2 ) + rhs.value( 1, 2 )
                            );
}

// ********************************************************************************

SymmetricMatrix3D operator-( const SymmetricMatrix3D & lhs, const SymmetricMatrix3D & rhs )
{
    return SymmetricMatrix3D( lhs.value( 0, 0 ) - rhs.value( 0, 0 ), lhs.value( 1, 1 ) - rhs.value( 1, 1 ), lhs.value( 2, 2 ) - rhs.value( 2, 2 ),
                              lhs.value( 0, 1 ) - rhs.value( 0, 1 ), lhs.value( 0, 2 ) - rhs.value( 0, 2 ), lhs.value( 1, 2 ) - rhs.value( 1, 2 )
                            );
}

// ********************************************************************************

SymmetricMatrix3D operator*( const double lhs, const SymmetricMatrix3D & rhs )
{
    return SymmetricMatrix3D( lhs * rhs.value( 0, 0 ), lhs * rhs.value( 1, 1 ), lhs * rhs.value( 2, 2 ),
                              lhs * rhs.value( 0, 1 ), lhs * rhs.value( 0, 2 ), lhs * rhs.value( 1, 2 ) );
}

// ********************************************************************************

SymmetricMatrix3D operator/( const SymmetricMatrix3D & lhs, const double rhs )
{
    return SymmetricMatrix3D( lhs.value( 0, 0 ) / rhs, lhs.value( 1, 1 ) / rhs, lhs.value( 2, 2 ) / rhs,
                              lhs.value( 0, 1 ) / rhs, lhs.value( 0, 2 ) / rhs, lhs.value( 1, 2 ) / rhs );
}

// ********************************************************************************

SymmetricMatrix3D inverse( const SymmetricMatrix3D & matrix3d )
{
    SymmetricMatrix3D result( matrix3d );
    result.invert();
    return result;
}

// ********************************************************************************

SymmetricMatrix3D transpose( const SymmetricMatrix3D & matrix3d )
{
    return matrix3d;
}
// ********************************************************************************

