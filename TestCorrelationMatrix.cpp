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

#include "CorrelationMatrix.h"
#include "TestSuite.h"

#include <iostream>

void test_correlation_matrix( TestSuite & test_suite )
{
    std::cout << "Now running tests for CorrelationMatrix." << std::endl;
    {
    CorrelationMatrix correlation_matrix( 5 );
    correlation_matrix.set_value( 0, 0,  0.0 );
    correlation_matrix.set_value( 0, 1,  1.0 );
    correlation_matrix.set_value( 0, 2,  2.0 );
    correlation_matrix.set_value( 0, 3,  3.0 );
    correlation_matrix.set_value( 0, 4,  4.0 );

    correlation_matrix.set_value( 1, 0, 10.0 );
    correlation_matrix.set_value( 1, 1, 11.0 );
    correlation_matrix.set_value( 1, 2, 12.0 );
    correlation_matrix.set_value( 1, 3, 13.0 );
    correlation_matrix.set_value( 1, 4, 14.0 );
    
    correlation_matrix.set_value( 2, 0, 20.0 );
    correlation_matrix.set_value( 2, 1, 21.0 );
    correlation_matrix.set_value( 2, 2, 22.0 );
    correlation_matrix.set_value( 2, 3, 23.0 );
    correlation_matrix.set_value( 2, 4, 24.0 );

    correlation_matrix.set_value( 3, 0, 30.0 );
    correlation_matrix.set_value( 3, 1, 31.0 );
    correlation_matrix.set_value( 3, 2, 32.0 );
    correlation_matrix.set_value( 3, 3, 33.0 );
    correlation_matrix.set_value( 3, 4, 34.0 );

    correlation_matrix.set_value( 4, 0, 40.0 );
    correlation_matrix.set_value( 4, 1, 41.0 );
    correlation_matrix.set_value( 4, 2, 42.0 );
    correlation_matrix.set_value( 4, 3, 43.0 );
    correlation_matrix.set_value( 4, 4, 44.0 );

    test_suite.test_equality( correlation_matrix.value( 0, 0 ),  1.0, "CorrelationMatrix::(set_)value() 01" );
    test_suite.test_equality( correlation_matrix.value( 0, 1 ), 10.0, "CorrelationMatrix::(set_)value() 02" );
    test_suite.test_equality( correlation_matrix.value( 0, 2 ), 20.0, "CorrelationMatrix::(set_)value() 03" );
    test_suite.test_equality( correlation_matrix.value( 0, 3 ), 30.0, "CorrelationMatrix::(set_)value() 04" );
    test_suite.test_equality( correlation_matrix.value( 0, 4 ), 40.0, "CorrelationMatrix::(set_)value() 05" );
    test_suite.test_equality( correlation_matrix.value( 1, 0 ), 10.0, "CorrelationMatrix::(set_)value() 06" );
    test_suite.test_equality( correlation_matrix.value( 1, 1 ),  1.0, "CorrelationMatrix::(set_)value() 07" );
    test_suite.test_equality( correlation_matrix.value( 1, 2 ), 21.0, "CorrelationMatrix::(set_)value() 08" );
    test_suite.test_equality( correlation_matrix.value( 1, 3 ), 31.0, "CorrelationMatrix::(set_)value() 09" );
    test_suite.test_equality( correlation_matrix.value( 1, 4 ), 41.0, "CorrelationMatrix::(set_)value() 10" );
    test_suite.test_equality( correlation_matrix.value( 2, 0 ), 20.0, "CorrelationMatrix::(set_)value() 11" );
    test_suite.test_equality( correlation_matrix.value( 2, 1 ), 21.0, "CorrelationMatrix::(set_)value() 12" );
    test_suite.test_equality( correlation_matrix.value( 2, 2 ),  1.0, "CorrelationMatrix::(set_)value() 13" );
    test_suite.test_equality( correlation_matrix.value( 2, 3 ), 32.0, "CorrelationMatrix::(set_)value() 14" );
    test_suite.test_equality( correlation_matrix.value( 2, 4 ), 42.0, "CorrelationMatrix::(set_)value() 15" );
    test_suite.test_equality( correlation_matrix.value( 3, 0 ), 30.0, "CorrelationMatrix::(set_)value() 16" );
    test_suite.test_equality( correlation_matrix.value( 3, 1 ), 31.0, "CorrelationMatrix::(set_)value() 17" );
    test_suite.test_equality( correlation_matrix.value( 3, 2 ), 32.0, "CorrelationMatrix::(set_)value() 18" );
    test_suite.test_equality( correlation_matrix.value( 3, 3 ),  1.0, "CorrelationMatrix::(set_)value() 19" );
    test_suite.test_equality( correlation_matrix.value( 3, 4 ), 43.0, "CorrelationMatrix::(set_)value() 20" );
    test_suite.test_equality( correlation_matrix.value( 4, 0 ), 40.0, "CorrelationMatrix::(set_)value() 21" );
    test_suite.test_equality( correlation_matrix.value( 4, 1 ), 41.0, "CorrelationMatrix::(set_)value() 22" );
    test_suite.test_equality( correlation_matrix.value( 4, 2 ), 42.0, "CorrelationMatrix::(set_)value() 23" );
    test_suite.test_equality( correlation_matrix.value( 4, 3 ), 43.0, "CorrelationMatrix::(set_)value() 24" );
    test_suite.test_equality( correlation_matrix.value( 4, 4 ),  1.0, "CorrelationMatrix::(set_)value() 25" );
    }

{
    CorrelationMatrix correlation_matrix( 5 );
    try
    {
        correlation_matrix.set_value( 5, 1, 40.0 );
        test_suite.log_error( "CorrelationMatrix::set_value() should have thrown 01." );
    }
    catch ( std::exception & e ) {}
}
{
    CorrelationMatrix correlation_matrix( 5 );
    try
    {
        correlation_matrix.set_value( 1, 5, 40.0 );
        test_suite.log_error( "CorrelationMatrix::set_value() should have thrown 02." );
    }
    catch ( std::exception & e ) {}
}
{
    CorrelationMatrix correlation_matrix( 5 );
    try
    {
        correlation_matrix.set_value( 67, 1, 40.0 );
        test_suite.log_error( "CorrelationMatrix::set_value() should have thrown 03." );
    }
    catch ( std::exception & e ) {}
}
{
    CorrelationMatrix correlation_matrix( 5 );
    try
    {
        correlation_matrix.set_value( 1, 67, 40.0 );
        test_suite.log_error( "CorrelationMatrix::set_value() should have thrown 04." );
    }
    catch ( std::exception & e ) {}
}

}

