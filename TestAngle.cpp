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

#include "Angle.h"

#include "TestSuite.h"

#include <iostream>

void test_angle( TestSuite & test_suite )
{
    std::cout << "Now running tests for Angle." << std::endl;
    {
    Angle angle;
    test_suite.test_equality_double( angle.value_in_radians(), 0.0, "Angle 1" );
    }
    {
        double x;
        double y;
        Angle angle;
        x = 1.0;
        y = 1.0;
        angle = ATAN2( y, x );
        test_suite.test_equality_double( angle.value_in_degrees(), 45.0, "Angle ATAN2() 1" );
        x = -1.0;
        y = 1.0;
        angle = ATAN2( y, x );
        test_suite.test_equality_double( angle.value_in_degrees(), 135.0, "Angle ATAN2() 2" );
        x = -1.0;
        y = -1.0;
        angle = ATAN2( y, x );
        test_suite.test_equality_double( angle.value_in_degrees(), 225.0, "Angle ATAN2() 3" );
        x = 1.0;
        y = -1.0;
        angle = ATAN2( y, x );
        test_suite.test_equality_double( angle.value_in_degrees(), 315.0, "Angle ATAN2() 4" );
    }
    
}

