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
    return 1.0;
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
    test_suite.test_equality_double( integral( &function_01, -1.0, 5.0, 0.1 ), 6.0, "integral 01" );
    }
}

