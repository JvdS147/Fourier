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

#include "Quaternion.h"

#include "TestSuite.h"

#include <string>
#include <iostream>

void test_one_quaternion( TestSuite & test_suite,
                          const Quaternion & actual_value,
                          const double a,
                          const double b,
                          const double c,
                          const double d,
                          const std::string & error_message )
{
    test_suite.test_equality_double( actual_value.a(), a, error_message + " (a)" );
    test_suite.test_equality_double( actual_value.b(), b, error_message + " (b)" );
    test_suite.test_equality_double( actual_value.c(), c, error_message + " (c)" );
    test_suite.test_equality_double( actual_value.d(), d, error_message + " (d)" );
}

void test_quaternion( TestSuite & test_suite )
{
    std::cout << "Now running tests for Quaternion." << std::endl;
    {
    Quaternion q;
    test_one_quaternion( test_suite, q, 1.0, 0.0, 0.0, 0.0, "Quaternion 1" );
    }
    {
    Quaternion q1( 0.6, 0.7, 0.3, -0.1 );
    Quaternion q2 = q1;
    q1.reciprocal();
    test_one_quaternion( test_suite, q1*q2, 1.0, 0.0, 0.0, 0.0, "Quaternion reciprocal()" );
    }
    {
    Quaternion q1( 0.6, 0.7, 0.3, -0.1 );
    Quaternion q2 = q1*q1;
    q1.square();
    if ( ! nearly_equal( q1, q2 ) )
        test_suite.log_error( "Quaternion square()" );
    }
    {
    Quaternion q1( 0.6, 0.7, 0.3, -0.1 );
    Quaternion q2 = q1*q1*q1;
    q1.power(3);
    if ( ! nearly_equal( q1, q2 ) )
        test_suite.log_error( "Quaternion power()" );
    }

    
}

