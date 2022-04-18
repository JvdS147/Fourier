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

#include "SetOfNumbers.h"

#include "TestSuite.h"

#include <iostream>

void test_SetOfNumbers( TestSuite & test_suite )
{
    std::cout << "Now running tests for SetOfNumbers." << std::endl;
    std::vector< size_t > values;
    values.push_back( 1 );
    values.push_back( 2 );
    values.push_back( 4 );
    values.push_back( 1 );
    {
    SetOfNumbers dummy( values, SetOfNumbers::ALLOWED );
    test_suite.test_equality( dummy.size(), static_cast<size_t>(4), std::string( "SetOfNumbers()" ) );
    test_suite.test_equality( dummy.value( 0 ), static_cast<size_t>(1), "SetOfNumbers()" );
    test_suite.test_equality( dummy.value( 1 ), static_cast<size_t>(1), "SetOfNumbers()" );
    test_suite.test_equality( dummy.value( 2 ), static_cast<size_t>(2), "SetOfNumbers()" );
    test_suite.test_equality( dummy.value( 3 ), static_cast<size_t>(4), "SetOfNumbers()" );
    dummy.add( 6 );
    test_suite.test_equality( dummy.size(), static_cast<size_t>(5), "SetOfNumbers()" );
    dummy.add( 4 );
    test_suite.test_equality( dummy.size(), static_cast<size_t>(6), "SetOfNumbers()" );
    test_suite.test_equality( dummy.value( 5 ), static_cast<size_t>(6), "SetOfNumbers()" );
    dummy.remove( 1 );
    test_suite.test_equality( dummy.size(), static_cast<size_t>(5), "SetOfNumbers()" );
    test_suite.test_equality( dummy.value( 1 ), static_cast<size_t>(2), "SetOfNumbers()" );
    test_suite.test_equality( dummy.contains( 1 ), true, "SetOfNumbers()" );
    test_suite.test_equality( dummy.contains( 7 ), false, "SetOfNumbers()" );
    }
    {
    SetOfNumbers dummy( values, SetOfNumbers::AUTO_REMOVE );
    test_suite.test_equality( dummy.size(), static_cast<size_t>(3), "SetOfNumbers()" );
    dummy.add( 6 );
    test_suite.test_equality( dummy.size(), static_cast<size_t>(4), "SetOfNumbers()" );
    dummy.add( 4 );
    test_suite.test_equality( dummy.size(), static_cast<size_t>(4), "SetOfNumbers()" );
    }
    try
    {
    SetOfNumbers dummy( values, SetOfNumbers::THROW );
    test_suite.log_error( "SetOfNumbers::SetOfNumbers() should have thrown " );
    }
    catch ( std::exception & e )
    {
    }
    try
    {
    SetOfNumbers dummy( values, SetOfNumbers::ALLOWED );
    dummy.value( 4 );
    test_suite.log_error( "SetOfNumbers::value() should have thrown " );
    }
    catch ( std::exception & e )
    {
    }
}

