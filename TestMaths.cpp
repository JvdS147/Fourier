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

#include "CrystalLattice.h"
#include "MathsFunctions.h"

#include "TestSuite.h"

#include <iostream>
#include <string>

double function_01( const double x )
{
    return exp( x );
}

void test_maths( TestSuite & test_suite )
{
    std::cout << "Now running tests for Maths." << std::endl;
    {
    double a = 6.884;
    double b = 9.569;
    double c = 7.093;
    Angle beta = Angle::from_degrees( 126.61 );
    CrystalLattice crystal_lattice_non_reduced( sqrt( 9.0*square( a ) + 4.0*square( c ) + 12.0*a*c*beta.cosine() ),
                                                b,
                                                c,
                                                Angle::angle_90_degrees(),
                                                arccotangent( beta.cotangent() + 2.0*c/( 3.0*a*beta.sine() ) ),
                                                Angle::angle_90_degrees() );
    test_suite.test_equality_double( crystal_lattice_non_reduced.a(), 16.6828, "test_maths() a", 0.0001 );
    test_suite.test_equality_double( crystal_lattice_non_reduced.b(), 9.569, "test_maths() b", 0.0001 );
    test_suite.test_equality_double( crystal_lattice_non_reduced.c(), 7.093, "test_maths() c", 0.0001 );
    test_suite.test_equality_double( crystal_lattice_non_reduced.alpha().value_in_degrees(), 90.0, "test_maths() alpha" );
    test_suite.test_equality_double( crystal_lattice_non_reduced.beta().value_in_degrees(),  83.5645, "test_maths() beta", 0.0001 );
    test_suite.test_equality_double( crystal_lattice_non_reduced.gamma().value_in_degrees(), 90.0, "test_maths() gamma" );
    }
    {
    test_suite.test_equality_double( integral( function_01, -2.0, 1.5, 30 ), 4.35163, "integral 01", 0.00001 );
 //   test_suite.integral_Monte_Carlo( integral( function_01, -2.0, 1.5, 30 ), 6.0, "integral 02" );
    test_suite.test_equality_double( integral_Gauss_Legendre_quadrature( function_01, -2.0, 1.5, 30 ), 4.34635, "integral 03", 0.00001 );
    }
    {
    std::vector< double > x;
    std::vector< double > w;
    Gauss_Legendre_quadrature( -1.0,  1.0, 3, x, w );
    test_suite.test_equality_double( x[0], -sqrt(3.0/5.0), "Gauss_Legendre_quadrature 01" );
    test_suite.test_equality_double( x[1],            0.0, "Gauss_Legendre_quadrature 02" );
    test_suite.test_equality_double( x[2],  sqrt(3.0/5.0), "Gauss_Legendre_quadrature 03" );
    test_suite.test_equality_double( w[0], 5.0/9.0, "Gauss_Legendre_quadrature 04" );
    test_suite.test_equality_double( w[1], 8.0/9.0, "Gauss_Legendre_quadrature 05" );
    test_suite.test_equality_double( w[2], 5.0/9.0, "Gauss_Legendre_quadrature 06" );
    }
    {
        double x = -3.0/17.0;
        double y1 = ( 231.0 * std::pow( x, 6 ) - 315.0 * std::pow( x, 4 ) + 105.0 * square( x ) - 5.0 ) / 16.0;
        double y2;
        double dydx;
        Legendre_polynomial_and_derivative( 6, x, y2, dydx );
        test_suite.test_equality_double( y1, y2, "Legendre_polynomial_and_derivative 01" );
    }
}

