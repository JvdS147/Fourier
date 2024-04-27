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

#include "LinearRegression.h"

#include "TestSuite.h"

#include <iostream>
//#include <string>

void test_linear_regression( TestSuite & test_suite )
{
    std::cout << "Now running tests for LinearRegression." << std::endl;
    {
    std::vector< double > x;
    std::vector< double > y;
    std::vector< double > s;
    x.push_back( 1.0 );
    y.push_back( 1.0 );
    s.push_back( 1.0 );
    x.push_back( 2.0 );
    y.push_back( 2.0 );
    s.push_back( 1.0 );
    x.push_back( 3.0 );
    y.push_back( 3.0 );
    s.push_back( 1.0 );
    x.push_back( 4.0 );
    y.push_back( 4.0 );
    s.push_back( 1.0 );
    x.push_back( 5.0 );
    y.push_back( 5.0 );
    s.push_back( 1.0 );
    double a;
    double b;
    linear_regression( x, y, s, a, b );
    test_suite.test_equality_double( a, 0.0, "test_linear_regression() a" );
    test_suite.test_equality_double( b, 1.0, "test_linear_regression() b" );
    }
    {
    std::vector< double > x;
    std::vector< double > y;
    x.push_back( 0.00278741 );
    x.push_back( 0.00836223 );
    x.push_back( 0.0139371 );
    x.push_back( 0.0195119 );
    x.push_back( 0.0250867 );
    x.push_back( 0.0306615 );
    x.push_back( 0.0362363 );
    x.push_back( 0.0418112 );
    x.push_back( 0.047386 );
    x.push_back( 0.0529608 );
    y.push_back( 6 );
    y.push_back( 15 );
    y.push_back( 38 );
    y.push_back( 104 );
    y.push_back( 227 );
    y.push_back( 397 );
    y.push_back( 754 );
    y.push_back( 1351 );
    y.push_back( 2671 );
    y.push_back( 4437 );
    double a;
    double b;
    fit_exponential( x, y, a, b );
    test_suite.test_equality_double( a, 6.0032, "fit_exponential() a", 0.0001 );
    test_suite.test_equality_double( b, 130.55, "fit_exponential() b", 0.01 );
    }
}

