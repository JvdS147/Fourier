/* *********************************************
Copyright (c) 2013-2023, Cornelis Jan (Jacco) van de Streek
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

#include "Utilities.h"
#include "TestSuite.h"

#include <iostream>
#include <string>

void test_utilities( TestSuite & test_suite )
{
    std::cout << "Now running tests for Utilities." << std::endl;

    {
        std::vector< size_t > result = initialise_with_sequential_values( static_cast< size_t >( 5 ) );
        test_suite.test_equality( result.size(), static_cast< size_t >( 5 ), "initialise_with_sequential_values() 01" );
        test_suite.test_equality( result[0], static_cast< size_t >( 0 ), "initialise_with_sequential_values() 02" );
        test_suite.test_equality( result[4], static_cast< size_t >( 4 ), "initialise_with_sequential_values() 03" );
    }

// Recognises scientific notation with "E" or "e" such as -.234e-45
// double string2double( std::string input );

// int string2integer( const std::string & input );

// std::string double2string( const double input );

// Pads the string to e.g. "0001"
// If the length of the input value is longer than the padded length, a string with
// the length of the input value is returned (so the value is never corrupted).
// std::string size_t2string( const size_t input, const size_t padded_length = 0, const char padding_character = '0' );

}

