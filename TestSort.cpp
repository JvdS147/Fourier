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

#include "Sort.h"

#include "TestSuite.h"

#include <iostream>

void test_sort( TestSuite & test_suite )
{
    std::cout << "Now running tests for sort()." << std::endl;
    std::vector< int > objects;
    objects.push_back( -1 );
    objects.push_back( -5 );
    objects.push_back(  0 );
    objects.push_back(  9 );
    objects.push_back(  5 );
    objects.push_back(  8 );
    objects.push_back(  7 );
    objects.push_back(  8 );
    objects.push_back(  9 );
    // -1,-5, 0, 9, 8, 7, 8, 9
    // -5,-1, 0, 7, 8, 8, 9, 9
    {
    Mapping sorted_map = sort( objects );
    test_suite.test_equality( objects[ sorted_map.map( 0 ) ], -5,  "Sort() 0" );
    test_suite.test_equality( objects[ sorted_map.map( 1 ) ], -1,  "Sort() 1" );
    test_suite.test_equality( objects[ sorted_map.map( 2 ) ],  0,  "Sort() 2" );
    test_suite.test_equality( objects[ sorted_map.map( 3 ) ],  5,  "Sort() 3" );
    test_suite.test_equality( objects[ sorted_map.map( 4 ) ],  7,  "Sort() 4" );
    test_suite.test_equality( objects[ sorted_map.map( 5 ) ],  8,  "Sort() 5" );
    test_suite.test_equality( objects[ sorted_map.map( 6 ) ],  8,  "Sort() 6" );
    test_suite.test_equality( objects[ sorted_map.map( 7 ) ],  9,  "Sort() 7" );
    test_suite.test_equality( objects[ sorted_map.map( 8 ) ],  9,  "Sort() 8" );
    }
    {
    Mapping sorted_map = sort( objects, true );
    test_suite.test_equality( objects[ sorted_map.map( 0 ) ],  9,  "Sort() 10" );
    test_suite.test_equality( objects[ sorted_map.map( 1 ) ],  9,  "Sort() 11" );
    test_suite.test_equality( objects[ sorted_map.map( 2 ) ],  8,  "Sort() 12" );
    test_suite.test_equality( objects[ sorted_map.map( 3 ) ],  8,  "Sort() 13" );
    test_suite.test_equality( objects[ sorted_map.map( 4 ) ],  7,  "Sort() 14" );
    test_suite.test_equality( objects[ sorted_map.map( 5 ) ],  5,  "Sort() 15" );
    test_suite.test_equality( objects[ sorted_map.map( 6 ) ],  0,  "Sort() 16" );
    test_suite.test_equality( objects[ sorted_map.map( 7 ) ], -1,  "Sort() 17" );
    test_suite.test_equality( objects[ sorted_map.map( 8 ) ], -5,  "Sort() 18" );
    }

}

