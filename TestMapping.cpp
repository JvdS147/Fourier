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

#include "Mapping.h"

#include "TestSuite.h"

#include <iostream>

void test_mapping( TestSuite & test_suite )
{
    std::cout << "Now running tests for Mapping." << std::endl;

    {
    Mapping dummy;
    dummy.push_back();
    test_suite.test_equality( dummy.size(), 1, "Mapping() 01" );
    test_suite.test_equality( dummy[0], 0, "Mapping() 02" );
    }
    {
    Mapping dummy( 5 );
    test_suite.test_equality( dummy[1], 1, "Mapping() 03" );
    test_suite.test_equality( dummy[4], 4, "Mapping() 04" );
    dummy.swap( 1, 4 );
    test_suite.test_equality( dummy[1], 4, "Mapping() 05" );
    test_suite.test_equality( dummy[4], 1, "Mapping() 06" );
    }
    {
    std::vector< size_t > mapping;
    mapping.push_back( 1 );
    mapping.push_back( 0 );
    mapping.push_back( 3 );
    mapping.push_back( 4 );
    mapping.push_back( 2 );
    Mapping dummy( mapping );
    test_suite.test_equality( dummy[1], 0, "Mapping() 07" );
    test_suite.test_equality( dummy[4], 2, "Mapping() 08" );
    dummy.invert();
    test_suite.test_equality( dummy[0], 1, "Mapping() 09" );
    test_suite.test_equality( dummy[2], 4, "Mapping() 10" );
    }

}

