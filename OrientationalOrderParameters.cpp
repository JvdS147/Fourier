/* *********************************************
Copyright (c) 2013-2023, Cornelis Jan (Jacco) van de Streek
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

#include "OrientationalOrderParameters.h"
#include "3DCalculations.h"
#include "CollectionOfPoints.h"
#include "Complex.h"
#include "MathsFunctions.h"
#include "NormalisedVector3D.h"
#include "Plane.h"
#include "SphericalHarmonics.h"
#include "Vector2D.h"

#include <iostream>
#include <stdexcept>

// ********************************************************************************

double orientational_order_parameter( const size_t l, const std::vector< NormalisedVector3D > & orientation_vectors )
{
    if ( orientation_vectors.empty() )
        throw std::runtime_error( "orientational_order_parameter(): must give at least one orientation." );
    std::vector< Angle > alphas;
    std::vector< Angle > betas;
    for ( size_t i( 0 ); i != orientation_vectors.size(); ++i )
    {
        Angle alpha;
        Angle beta;
        Eulerian_angles( orientation_vectors[i], alpha, beta );
        alphas.push_back( alpha );
        betas.push_back( beta );                
    }
    // m = 0
    double a_l_0 = 0.0;
    for ( size_t i( 0 ); i != orientation_vectors.size(); ++i )
        a_l_0 += Legendre_polynomial( l, betas[i].cosine() );
    double S = square( a_l_0 );
    // + 2 * |m|
    for ( int m( 1 ); m != (l+1); ++m )
    {
        Complex a_l_m;
        for ( size_t i( 0 ); i != orientation_vectors.size(); ++i )
        {
            Complex Y_l_m = Racah_spherical_harmonics( l, m, alphas[i], betas[i] );
            a_l_m += Y_l_m;
        }
        S += 2.0 * a_l_m.norm2();
    }
    S /= square( orientation_vectors.size() );
    if ( nearly_zero( S ) )
        S = 0.0;
    return S;
}

// ********************************************************************************

double orientational_order_parameter_S6( const std::vector< Vector3D > & points )
{
    if ( points.size() < 3 )
        throw std::runtime_error( "orientational_order_parameter_S6(): must give at least three points." );
    if ( points.size() != 6 )
        std::cout << "orientational_order_parameter_S6(): warning: number of points is not six." << std::endl;
    CollectionOfPoints points_2( points );
    points_2.move_to_centre_of_mass();
    Plane plane = ::plane( points_2 );
    Complex result;
    for ( size_t i( 0 ); i != points.size(); ++i )
    {
        Vector2D point2 = projection( plane, points_2.point( i ) );
        result += exponential( Complex( 0.0, -6.0 * ATAN2( point2.y(), point2.x() ).value_in_radians() ) );
    }
    return result.norm2() / square( points.size() );
}

// ********************************************************************************

