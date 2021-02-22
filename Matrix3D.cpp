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

#include "Matrix3D.h"
#include "Utilities.h"

#include <cmath>
#include <stdexcept>
#include <iostream>

// ********************************************************************************

Matrix3D::Matrix3D()
{
    *this = Matrix3D( 1.0 );
}

// ********************************************************************************

// Constructs a scalar matrix
Matrix3D::Matrix3D( const double scalar )
{
    data_[0][0] = scalar; data_[0][1] = 0.0   ; data_[0][2] = 0.0   ;
    data_[1][0] = 0.0   ; data_[1][1] = scalar; data_[1][2] = 0.0   ;
    data_[2][0] = 0.0   ; data_[2][1] = 0.0   ; data_[2][2] = scalar;
}

// ********************************************************************************

Matrix3D::Matrix3D( const double a00, const double a01, const double a02,
                    const double a10, const double a11, const double a12,
                    const double a20, const double a21, const double a22
                  )
{
    data_[0][0] = a00; data_[0][1] = a01; data_[0][2] = a02;
    data_[1][0] = a10; data_[1][1] = a11; data_[1][2] = a12;
    data_[2][0] = a20; data_[2][1] = a21; data_[2][2] = a22;
}

// ********************************************************************************

double Matrix3D::sum_of_elements() const
{
    return data_[0][0] + data_[0][1] + data_[0][2] + data_[1][0] + data_[1][1] + data_[1][2] + data_[2][0] + data_[2][1] + data_[2][2];
}

// ********************************************************************************

double Matrix3D::sum_of_absolute_elements() const
{
    return std::abs( data_[0][0] ) + std::abs( data_[0][1] ) + std::abs( data_[0][2] ) + std::abs( data_[1][0] ) + std::abs( data_[1][1] ) + std::abs( data_[1][2] ) + std::abs( data_[2][0] ) + std::abs( data_[2][1] ) + std::abs( data_[2][2] );
}

// ********************************************************************************

void Matrix3D::invert()
{
    double D = determinant();
    if ( nearly_zero( D ) )
        std::runtime_error( "Matrix3D::invert(): determinant = 0" );
    transpose();
    Matrix3D adjoint;
    adjoint.data_[0][0] = +minor_matrix_determinant(0,0); adjoint.data_[0][1] = -minor_matrix_determinant(0,1); adjoint.data_[0][2] = +minor_matrix_determinant(0,2);
    adjoint.data_[1][0] = -minor_matrix_determinant(1,0); adjoint.data_[1][1] = +minor_matrix_determinant(1,1); adjoint.data_[1][2] = -minor_matrix_determinant(1,2);
    adjoint.data_[2][0] = +minor_matrix_determinant(2,0); adjoint.data_[2][1] = -minor_matrix_determinant(2,1); adjoint.data_[2][2] = +minor_matrix_determinant(2,2);
    adjoint /= D;
    *this = adjoint;
}

// ********************************************************************************

void Matrix3D::transpose()
{
    std::swap( data_[1][0], data_[0][1] );
    std::swap( data_[2][0], data_[0][2] );
    std::swap( data_[2][1], data_[1][2] );
}

// ********************************************************************************

double Matrix3D::determinant() const
{
    return data_[0][0] * minor_matrix_determinant(0,0) -
           data_[0][1] * minor_matrix_determinant(0,1) +
           data_[0][2] * minor_matrix_determinant(0,2);
}

// ********************************************************************************

double Matrix3D::trace() const
{
    return data_[0][0] + data_[1][1] + data_[2][2];
}

// ********************************************************************************

void Matrix3D::swap_rows( const size_t i, const size_t j )
{
    std::swap( data_[i][0], data_[j][0] );
    std::swap( data_[i][1], data_[j][1] );
    std::swap( data_[i][2], data_[j][2] );
}

// ********************************************************************************

void Matrix3D::swap_columns( const size_t i, const size_t j )
{
    std::swap( data_[0][i], data_[0][j] );
    std::swap( data_[1][i], data_[1][j] );
    std::swap( data_[2][i], data_[2][j] );    
}

// ********************************************************************************

double Matrix3D::maximum_absolute_value_in_row( const size_t i ) const
{
    return std::max( std::abs( data_[i][0] ), std::max( std::abs( data_[i][1] ), std::abs( data_[i][2] ) ) );
}

// ********************************************************************************

double Matrix3D::maximum_absolute_value_in_column( const size_t i ) const
{
    return std::max( std::abs( data_[0][i] ), std::max( std::abs( data_[1][i] ), std::abs( data_[2][i] ) ) );    
}

// ********************************************************************************

bool Matrix3D::is_diagonal() const
{
    return ( nearly_zero( data_[0][1] ) && nearly_zero( data_[0][2] ) && 
             nearly_zero( data_[1][0] ) && nearly_zero( data_[1][2] ) && 
             nearly_zero( data_[2][0] ) && nearly_zero( data_[2][1] ) );
}

// ********************************************************************************

double Matrix3D::minor_matrix_determinant( const size_t row, const size_t col ) const
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
            minor_matrix[minor_row][minor_col] = data_[i][j];
            ++minor_col;
        }
        ++minor_row;
    }
    return minor_matrix[0][0]*minor_matrix[1][1] - minor_matrix[0][1]*minor_matrix[1][0];
}

// ********************************************************************************

void Matrix3D::convert_to_row_echelon_form( Matrix3D & T )
{
    T = Matrix3D();
    size_t iRow( 0 );
    size_t iCol( 0 );
    while ( ( iRow < 2 ) && ( iCol < 3 ) )
    {
        // Find the pivot in column iCol.
        size_t iMax( iRow );
        for ( size_t i( iRow + 1 ); i != 3; ++i )
        {
            if ( std::abs( data_[i][iCol] ) > std::abs( data_[iMax][iCol] ) )
                iMax = i;
        }
        if ( ! nearly_zero( data_[iMax][iCol] ) )
        {
            swap_rows( iRow, iMax );
            T.swap_rows( iRow, iMax );
            // Do for all rows below pivot.
            for ( size_t i( iRow + 1 ); i != 3; ++i )
            {
                double factor = -1.0 * ( data_[i][iCol] / data_[iRow][iCol] );
                data_[i][iCol] = 0.0;
                for ( size_t j( iCol + 1 ); j != 3; ++j )
                {
                    data_[i][j] += factor * data_[iRow][j];
                    T.set_value( i, j, T.value( i, j ) + factor * value( iRow, j ) );
                }
            }
            ++iRow;
        }
        ++iCol;
    }
}

// ********************************************************************************

size_t Matrix3D::number_of_zero_rows( const double tolerance ) const
{
    size_t result( 0 );
    for ( size_t i( 0 ); i != 3; ++i )
    {
        if ( nearly_zero( data_[i][0] ) && nearly_zero( data_[i][1] ) && nearly_zero( data_[i][2] ) )
            ++result;
    }
    return result;
}

// ********************************************************************************

Matrix3D & Matrix3D::operator+=( const Matrix3D & rhs )
{
    *this = Matrix3D( value( 0, 0 ) + rhs.value( 0, 0 ), value( 0, 1 ) + rhs.value( 0, 1 ), value( 0, 2 ) + rhs.value( 0, 2 ),
                      value( 1, 0 ) + rhs.value( 1, 0 ), value( 1, 1 ) + rhs.value( 1, 1 ), value( 1, 2 ) + rhs.value( 1, 2 ),
                      value( 2, 0 ) + rhs.value( 2, 0 ), value( 2, 1 ) + rhs.value( 2, 1 ), value( 2, 2 ) + rhs.value( 2, 2 )
                    );
    return *this;
}

// ********************************************************************************

Matrix3D & Matrix3D::operator-=( const Matrix3D & rhs )
{
    *this = Matrix3D( value( 0, 0 ) - rhs.value( 0, 0 ), value( 0, 1 ) - rhs.value( 0, 1 ), value( 0, 2 ) - rhs.value( 0, 2 ),
                      value( 1, 0 ) - rhs.value( 1, 0 ), value( 1, 1 ) - rhs.value( 1, 1 ), value( 1, 2 ) - rhs.value( 1, 2 ),
                      value( 2, 0 ) - rhs.value( 2, 0 ), value( 2, 1 ) - rhs.value( 2, 1 ), value( 2, 2 ) - rhs.value( 2, 2 )
                    );
    return *this;
}

// ********************************************************************************

Matrix3D & Matrix3D::operator/=( const double value )
{
    data_[0][0] /= value; data_[0][1] /= value; data_[0][2] /= value;
    data_[1][0] /= value; data_[1][1] /= value; data_[1][2] /= value;
    data_[2][0] /= value; data_[2][1] /= value; data_[2][2] /= value;
    return *this;
}

// ********************************************************************************

bool Matrix3D::operator==( const Matrix3D & rhs ) const
{
    return nearly_equal( *this, rhs );
}

// ********************************************************************************

void Matrix3D::show() const
{
    std::cout << "| " << data_[0][0] << " " << data_[0][1] << " " << data_[0][2] << " |" << std::endl;
    std::cout << "| " << data_[1][0] << " " << data_[1][1] << " " << data_[1][2] << " |" << std::endl;
    std::cout << "| " << data_[2][0] << " " << data_[2][1] << " " << data_[2][2] << " |" << std::endl;
}

// ********************************************************************************

std::ostream & operator<<( std::ostream & os, const Matrix3D & matrix3d )
{
    os << "| " << matrix3d.value( 0, 0 ) << " " << matrix3d.value( 0, 1 ) << " " << matrix3d.value( 0, 2 ) << " |" << std::endl;
    os << "| " << matrix3d.value( 1, 0 ) << " " << matrix3d.value( 1, 1 ) << " " << matrix3d.value( 1, 2 ) << " |" << std::endl;
    os << "| " << matrix3d.value( 2, 0 ) << " " << matrix3d.value( 2, 1 ) << " " << matrix3d.value( 2, 2 ) << " |" << std::endl;
    return os;
}

// ********************************************************************************

bool nearly_equal( const Matrix3D & lhs, const Matrix3D & rhs, const double tolerance )
{
    return ( nearly_equal( lhs.value( 0, 0 ), rhs.value( 0, 0 ), tolerance ) && nearly_equal( lhs.value( 0, 1 ), rhs.value( 0, 1 ), tolerance ) && nearly_equal( lhs.value( 0, 2 ), rhs.value( 0, 2 ), tolerance ) &&
             nearly_equal( lhs.value( 1, 0 ), rhs.value( 1, 0 ), tolerance ) && nearly_equal( lhs.value( 1, 1 ), rhs.value( 1, 1 ), tolerance ) && nearly_equal( lhs.value( 1, 2 ), rhs.value( 1, 2 ), tolerance ) &&
             nearly_equal( lhs.value( 2, 0 ), rhs.value( 2, 0 ), tolerance ) && nearly_equal( lhs.value( 2, 1 ), rhs.value( 2, 1 ), tolerance ) && nearly_equal( lhs.value( 2, 2 ), rhs.value( 2, 2 ), tolerance ) );
}

// ********************************************************************************

Matrix3D operator+( const Matrix3D & lhs, const Matrix3D & rhs )
{
    return Matrix3D( lhs.value( 0, 0 ) + rhs.value( 0, 0 ), lhs.value( 0, 1 ) + rhs.value( 0, 1 ), lhs.value( 0, 2 ) + rhs.value( 0, 2 ),
                     lhs.value( 1, 0 ) + rhs.value( 1, 0 ), lhs.value( 1, 1 ) + rhs.value( 1, 1 ), lhs.value( 1, 2 ) + rhs.value( 1, 2 ),
                     lhs.value( 2, 0 ) + rhs.value( 2, 0 ), lhs.value( 2, 1 ) + rhs.value( 2, 1 ), lhs.value( 2, 2 ) + rhs.value( 2, 2 )
                   );
}

// ********************************************************************************

Matrix3D operator-( const Matrix3D & lhs, const Matrix3D & rhs )
{
    return Matrix3D( lhs.value( 0, 0 ) - rhs.value( 0, 0 ), lhs.value( 0, 1 ) - rhs.value( 0, 1 ), lhs.value( 0, 2 ) - rhs.value( 0, 2 ),
                     lhs.value( 1, 0 ) - rhs.value( 1, 0 ), lhs.value( 1, 1 ) - rhs.value( 1, 1 ), lhs.value( 1, 2 ) - rhs.value( 1, 2 ),
                     lhs.value( 2, 0 ) - rhs.value( 2, 0 ), lhs.value( 2, 1 ) - rhs.value( 2, 1 ), lhs.value( 2, 2 ) - rhs.value( 2, 2 )
                   );
}

// ********************************************************************************

Matrix3D operator*( const Matrix3D & lhs, const Matrix3D & rhs )
{
    Matrix3D result;
    result.set_value( 0, 0, lhs.value(0,0) * rhs.value(0,0) + lhs.value(0,1) * rhs.value(1,0) + lhs.value(0,2) * rhs.value(2,0) );
    result.set_value( 0, 1, lhs.value(0,0) * rhs.value(0,1) + lhs.value(0,1) * rhs.value(1,1) + lhs.value(0,2) * rhs.value(2,1) );
    result.set_value( 0, 2, lhs.value(0,0) * rhs.value(0,2) + lhs.value(0,1) * rhs.value(1,2) + lhs.value(0,2) * rhs.value(2,2) );
    result.set_value( 1, 0, lhs.value(1,0) * rhs.value(0,0) + lhs.value(1,1) * rhs.value(1,0) + lhs.value(1,2) * rhs.value(2,0) );
    result.set_value( 1, 1, lhs.value(1,0) * rhs.value(0,1) + lhs.value(1,1) * rhs.value(1,1) + lhs.value(1,2) * rhs.value(2,1) );
    result.set_value( 1, 2, lhs.value(1,0) * rhs.value(0,2) + lhs.value(1,1) * rhs.value(1,2) + lhs.value(1,2) * rhs.value(2,2) );
    result.set_value( 2, 0, lhs.value(2,0) * rhs.value(0,0) + lhs.value(2,1) * rhs.value(1,0) + lhs.value(2,2) * rhs.value(2,0) );
    result.set_value( 2, 1, lhs.value(2,0) * rhs.value(0,1) + lhs.value(2,1) * rhs.value(1,1) + lhs.value(2,2) * rhs.value(2,1) );
    result.set_value( 2, 2, lhs.value(2,0) * rhs.value(0,2) + lhs.value(2,1) * rhs.value(1,2) + lhs.value(2,2) * rhs.value(2,2) );
    return result;
}

// ********************************************************************************

Matrix3D operator*( const double lhs, const Matrix3D & rhs )
{
    return Matrix3D( lhs * rhs.value( 0, 0 ), lhs * rhs.value( 0, 1 ), lhs * rhs.value( 0, 2 ),
                     lhs * rhs.value( 1, 0 ), lhs * rhs.value( 1, 1 ), lhs * rhs.value( 1, 2 ),
                     lhs * rhs.value( 2, 0 ), lhs * rhs.value( 2, 1 ), lhs * rhs.value( 2, 2 ) );
}

// ********************************************************************************

Matrix3D operator/( const Matrix3D & lhs, const double rhs )
{
    return Matrix3D( lhs.value( 0, 0 ) / rhs, lhs.value( 0, 1 ) / rhs, lhs.value( 0, 2 ) / rhs,
                     lhs.value( 1, 0 ) / rhs, lhs.value( 1, 1 ) / rhs, lhs.value( 1, 2 ) / rhs,
                     lhs.value( 2, 0 ) / rhs, lhs.value( 2, 1 ) / rhs, lhs.value( 2, 2 ) / rhs );
}

// ********************************************************************************

Matrix3D inverse( const Matrix3D & matrix3d )
{
    Matrix3D result( matrix3d );
    result.invert();
    return result;
}

// ********************************************************************************

Matrix3D transpose( const Matrix3D & matrix3d )
{
    Matrix3D result( matrix3d );
    result.transpose();
    return result;
}
// ********************************************************************************

