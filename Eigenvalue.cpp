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

#include "Eigenvalue.h"

#include "3DCalculations.h"
#include "MathFunctions.h"
#include "Matrix3D.h"
#include "NormalisedVector3D.h"
#include "SymmetricMatrix3D.h"
#include "Vector3D.h"

#include <stdexcept>
#include <cmath>

// ********************************************************************************

/** Eigenvalues and eigenvectors of a real matrix.
    If A is symmetric, then A = V*D*V' where the eigenvalue matrix D is
    diagonal and the eigenvector matrix V is orthogonal.
    I.e. A = V.times(D.times(V.transpose())) and
    V.times(V.transpose()) equals the identity matrix.
**/
void calculate_eigenvalues( const SymmetricMatrix3D & input, std::vector< double > & eigenvalues, std::vector< NormalisedVector3D > & eigenvectors )
{
    Matrix3D V = SymmetricMatrix3D2Matrix3D( input ); // Array for internal storage of eigenvectors.
    Vector3D d; // Arrays for internal storage of eigenvalues.
    Vector3D e; // Arrays for internal storage of eigenvalues.
    // Tridiagonalise.
    tred2( V, d, e );
    // Diagonalise.
    tql2( V, d, e );
    eigenvalues.clear();
    eigenvectors.clear();
    for ( size_t i( 0 ); i != 3; ++i )
        eigenvalues.push_back( d.value( i ) );
    eigenvectors.push_back( NormalisedVector3D( V.value(0,0), V.value(1,0), V.value(2,0) ) );
    eigenvectors.push_back( NormalisedVector3D( V.value(0,1), V.value(1,1), V.value(2,1) ) );
    eigenvectors.push_back( NormalisedVector3D( V.value(0,2), V.value(1,2), V.value(2,2) ) );
}

// ********************************************************************************

// Symmetric Householder reduction to tridiagonal form.
void tred2( Matrix3D & V, Vector3D & d, Vector3D & e )
{
    for ( size_t j = 0; j < 3; ++j )
        d.set_value( j, V.value( 2, j ) );
    // Householder reduction to tridiagonal form.
    for ( int i = 2; i > 0; --i )
    {
        // Scale to avoid under/overflow.
        double scale = 0.0;
        double h = 0.0;
        for ( int k = 0; k < i; ++k )
            scale += std::fabs( d.value( k ) );
        if ( scale == 0.0 )
        {
            e.set_value( i, d.value( i-1 ) );
            for ( int j = 0; j < i; ++j )
            {
                d.set_value( j, V.value( i-1, j ) );
                V.set_value( i, j, 0.0 );
                V.set_value( j, i, 0.0 );
            }
        }
        else
        {
            // Generate Householder vector.
            for ( int k = 0; k < i; ++k )
            {
                d.set_value( k, d.value( k ) / scale );
                h += square( d.value( k ) );
            }
            double f = d.value( i-1 );
            double g = sqrt( h );
            if ( f > 0.0 )
                g = -g;
            e.set_value( i, scale * g );
            h = h - f * g;
            d.set_value( i-1, f - g );
            for ( int j = 0; j < i; ++j )
                e.set_value( j, 0.0 );
            // Apply similarity transformation to remaining columns.
            for ( int j = 0; j < i; ++j )
            {
                f = d.value( j );
                V.set_value( j, i, f );
                g = e.value( j ) + V.value( j, j ) * f;
                for ( int k = j+1; k <= i-1; ++k )
                {
                    g += V.value( k, j ) * d.value( k );
                    e.set_value( k, e.value( k ) + V.value( k, j ) * f );
                }
                e.set_value( j, g );
            }
            f = 0.0;
            for ( int j = 0; j < i; ++j )
            {
                e.set_value( j, e.value( j ) / h );
                f += e.value( j ) * d.value( j );
            }
            double hh = f / (h + h);
            for ( int j = 0; j < i; ++j )
                e.set_value( j, e.value( j ) - hh * d.value( j ) );
            for ( int j = 0; j < i; ++j )
            {
                f = d.value( j );
                g = e.value( j );
                for ( int k = j; k <= i-1; ++k )
                    V.set_value( k, j, V.value(k, j ) - ( f * e.value( k ) + g * d.value( k ) ) );
                d.set_value( j, V.value( i-1, j ) );
                V.set_value( i, j, 0.0 );
            }
        }
        d.set_value( i, h );
    }
    // Accumulate transformations.
    for ( int i = 0; i < 2; ++i )
    {
         V.set_value( 2, i, V.value( i, i ) );
         V.set_value( i, i, 1.0 );
         double h = d.value( i+1 );
         if ( h != 0.0 )
         {
            for ( int k = 0; k <= i; ++k )
               d.set_value( k, V.value( k, i+1 ) / h );
            for ( int j = 0; j <= i; ++j )
            {
               double g = 0.0;
               for ( int k = 0; k <= i; ++k )
                  g += V.value( k, i+1 ) * V.value( k, j );
               for ( int k = 0; k <= i; ++k )
                  V.set_value( k, j, V.value( k, j ) - g * d.value( k ) );
            }
         }
         for ( int k = 0; k <= i; ++k )
            V.set_value( k, i+1, 0.0 );
      }
      for ( int j = 0; j < 3; ++j )
      {
         d.set_value( j, V.value( 2, j ) );
         V.set_value( 2, j, 0.0 );
      }
      V.set_value( 2, 2, 1.0 );
      e.set_value( 0, 0.0 );
}

// ********************************************************************************

// Symmetric tridiagonal QL algorithm.
void tql2( Matrix3D & V, Vector3D & d, Vector3D & e )
{

   //  This is derived from the Algol procedures tql2, by
   //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
   //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
   //  Fortran subroutine in EISPACK.

    for ( int i = 1; i < 3; ++i )
        e.set_value( i-1, e.value( i ) );
    e.set_value( 2, 0.0 );
    double f = 0.0;
    double tst1 = 0.0;
    double eps = std::pow( 2.0, -52.0 );
    for ( int l = 0; l < 3; ++l )
    {
        // Find small subdiagonal element
        tst1 = std::max( tst1, fabs( d.value( l ) ) + fabs( e.value( l ) ) );
        int m = l;
        while ( m < 3 )
        {
            if ( fabs( e.value( m ) ) <= eps * tst1 )
                break;
            ++m;
        }
        // If m == l, d[l] is an eigenvalue,
        // otherwise, iterate.
        if ( m > l )
        {
            int iter = 0;
            do
            {
                ++iter;
                if ( iter > 100 )
                    throw std::runtime_error( "tql2() : too many iterations." );
                // Compute implicit shift
                double g = d.value( l );
                double p = ( d.value( l+1 ) - g ) / ( 2.0 * e.value( l ) );
                double r = hypothenuse( p, 1.0 );
                if ( p < 0 )
                    r = -r;
                d.set_value( l, e.value( l ) / (p + r) );
                d.set_value( l+1, e.value( l ) * (p + r) );
                double dl1 = d.value( l+1 );
                double h = g - d.value( l );
                for ( int i = l+2; i < 3; ++i )
                    d.set_value( i, d.value( i ) - h );
                f += h;
                // Implicit QL transformation.
                p = d.value( m );
                double c = 1.0;
                double c2 = c;
                double c3 = c;
                double el1 = e.value( l+1 );
                double s = 0.0;
                double s2 = 0.0;
                for ( int i = m-1; i >= l; --i )
                {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e.value( i );
                    h = c * p;
                    r = hypothenuse( p, e.value( i ) );
                    e.set_value( i+1, s * r );
                    s = e.value( i ) / r;
                    c = p / r;
                    p = c * d.value( i ) - s * g;
                    d.set_value( i+1, h + s * ( c * g + s * d.value( i ) ) );
                    // Accumulate transformation.
                    for ( int k = 0; k < 3; ++ k )
                    {
                        h = V.value( k, i+1 );
                        V.set_value( k, i+1, s * V.value( k, i) + c * h );
                        V.set_value( k, i, c * V.value( k, i ) - s * h );
                    }
                }
                p = -s * s2 * c3 * el1 * e.value( l ) / dl1;
                e.set_value( l, s * p );
                d.set_value( l, c * p );
                // Check for convergence.
            }
            while ( fabs( e.value( l ) ) > eps * tst1 );
        }
        d.set_value( l, d.value( l ) + f );
        e.set_value( l, 0.0 );
    }
    // Sort eigenvalues and corresponding vectors.
    for ( int i = 0; i < 2; ++i )
    {
        int k = i;
        double p = d.value( i );
        for ( int j = i+1; j < 3; ++j )
        {
            if ( d.value( j ) < p )
            {
                k = j;
                p = d.value( j );
            }
        }
        if ( k != i )
        {
            d.set_value( k, d.value( i ) );
            d.set_value( i, p );
            for ( int j = 0; j < 3; ++j )
            {
                p = V.value( j, i );
                V.set_value( j, i, V.value( j, k ) );
                V.set_value( j, k, p );
            }
        }
    }
}

// ********************************************************************************

