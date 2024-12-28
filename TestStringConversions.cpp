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

#include "StringConversions.h"
#include "TestSuite.h"
#include "Angle.h"
#include "Matrix3D.h"
#include "MillerIndices.h"
#include "Vector3D.h"

#include <iostream>
#include <string>

void test_StringConversions( TestSuite & test_suite )
{
    std::cout << "Now running tests for StringConversions." << std::endl;

    {
        Angle input( 181.4, Angle::DEGREES );
        std::string result = to_string( input );
        Angle output = Angle_from_string( result );
        if ( ! nearly_equal( input, output ) )
            test_suite.log_error( "Angle_from_string( to_string( Angle ) )" );
    }
    {
        Matrix3D input( 1.2556, -3.4778 , 5.6992,
                        0.0   ,  9.47325, 2.453 ,
                        3.562 ,  4.279  , 3.4805 );
        std::string result = to_string( input );
        Matrix3D output = Matrix3D_from_string( result );
        if ( ! nearly_equal( input, output ) )
            test_suite.log_error( "Matrix3D_from_string( to_string( Matrix3d ) )" );
    }
    {
        MillerIndices input( -1, 2, 0 );
        std::string result = to_string( input );
        MillerIndices output = MillerIndices_from_string( result );
        if ( input != output )
            test_suite.log_error( "MillerIndices_from_string( to_string( MillerIndices ) )" );
    }
    {
        Vector3D input( 1.2556, -3.4778, 5.6992 );
        std::string result = to_string( input );
        Vector3D output = Vector3D_from_string( result );
        if ( ! nearly_equal( input, output ) )
            test_suite.log_error( "Vector3D_from_string( to_string( Vector3D ) )" );
    }

}

