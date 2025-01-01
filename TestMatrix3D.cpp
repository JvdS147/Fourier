/* *********************************************
Copyright (c) 2013-2025, Cornelis Jan (Jacco) van de Streek
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

#include "TestSuite.h"

#include <iostream>
#include <string>

void test_one_matrix3D( TestSuite & test_suite,
                        const Matrix3D & matrix3D,
                        const double value00,
                        const double value01,
                        const double value02,
                        const double value10,
                        const double value11,
                        const double value12,
                        const double value20,
                        const double value21,
                        const double value22,
                        const std::string& error_message )
{
    test_suite.test_equality_double( matrix3D.value( 0, 0 ), value00, error_message + " element 0 0" );
    test_suite.test_equality_double( matrix3D.value( 0, 1 ), value01, error_message + " element 0 1" );
    test_suite.test_equality_double( matrix3D.value( 0, 2 ), value02, error_message + " element 0 2" );
    test_suite.test_equality_double( matrix3D.value( 1, 0 ), value10, error_message + " element 1 0" );
    test_suite.test_equality_double( matrix3D.value( 1, 1 ), value11, error_message + " element 1 1" );
    test_suite.test_equality_double( matrix3D.value( 1, 2 ), value12, error_message + " element 1 2" );
    test_suite.test_equality_double( matrix3D.value( 2, 0 ), value20, error_message + " element 2 0" );
    test_suite.test_equality_double( matrix3D.value( 2, 1 ), value21, error_message + " element 2 1" );
    test_suite.test_equality_double( matrix3D.value( 2, 2 ), value22, error_message + " element 2 2" );
}

void test_matrix3D( TestSuite & test_suite )
{
    std::cout << "Now running tests for Matrix3D." << std::endl;
    {
    Matrix3D matrix3D;
    test_one_matrix3D( test_suite, matrix3D, 1.0, 0.0, 0.0,
                                             0.0, 1.0, 0.0,
                                             0.0, 0.0, 1.0, "Matrix3D::Matrix3D() 1" );
    }
    {
    Matrix3D matrix3D( 2.0, 0.0, 0.0,
                       0.0, 3.0, 0.0,
                       0.0, 0.0, 4.0 );
    test_one_matrix3D( test_suite, matrix3D, 2.0, 0.0, 0.0,
                                             0.0, 3.0, 0.0,
                                             0.0, 0.0, 4.0, "Matrix3D::Matrix3D() 2" );
    }
    {
    Matrix3D matrix3D( 2.0, 0.0, 0.0,
                       0.0, 3.0, 0.0,
                       0.0, 0.0, 4.0 );
    Matrix3D matrix3D_2( matrix3D );
    matrix3D.invert();
    test_one_matrix3D( test_suite, matrix3D, 1.0/2.0, 0.0,     0.0,
                                             0.0,     1.0/3.0, 0.0,
                                             0.0,     0.0,     1.0/4.0, "Matrix3D::invert() 1" );
    test_one_matrix3D( test_suite, matrix3D * matrix3D_2, 1.0, 0.0, 0.0,
                                                          0.0, 1.0, 0.0,
                                                          0.0, 0.0, 1.0, "Matrix3D::invert() 2" );

    }
    {
    Matrix3D matrix3D( 2.0, 0.0, 0.0,
                       0.0, 3.0, 0.0,
                       0.0, 0.0, 4.0 );
    Matrix3D T;
    matrix3D.convert_to_row_echelon_form( T );
    test_one_matrix3D( test_suite, matrix3D, 2.0, 0.0, 0.0,
                                             0.0, 3.0, 0.0,
                                             0.0, 0.0, 4.0, "Matrix3D::convert_to_row_echelon_form() 1" );
    }
    {
    Matrix3D matrix3D( 1.0, 0.0, 0.0,
                       0.0, 0.0, 0.0,
                       0.0, 0.0, -1.0 );
    Matrix3D T;
    matrix3D.convert_to_row_echelon_form( T );
    test_one_matrix3D( test_suite, matrix3D, 1.0, 0.0, 0.0,
                                             0.0, 0.0, -1.0,
                                             0.0, 0.0, 0.0, "Matrix3D::convert_to_row_echelon_form() 2" );
    }
    {
    Matrix3D matrix3D( 2.0, 0.0, 0.0,
                       1.0, 0.0, 0.0,
                       0.0, -1.0, 0.0 );
    Matrix3D T;
    matrix3D.convert_to_row_echelon_form( T );
    test_one_matrix3D( test_suite, matrix3D, 2.0, 0.0, 0.0,
                                             0.0, -1.0, 0.0,
                                             0.0, 0.0, 0.0, "Matrix3D::convert_to_row_echelon_form() 3" );
    }
    
}

