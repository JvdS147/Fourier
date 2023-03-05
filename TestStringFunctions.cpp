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

#include "StringFunctions.h"
#include "TestSuite.h"

#include <iostream>
#include <string>

void test_StringFunctions( TestSuite & test_suite )
{
    std::cout << "Now running tests for StringFunctions." << std::endl;

// char to_upper( const char argument );

// char to_lower( const char argument );

    {
        std::string result = to_upper( "" );
        test_suite.test_equality( result, std::string( "" ), "to_upper() 01" );
    }
    {
        std::string result = to_lower( "" );
        test_suite.test_equality( result, std::string( "" ), "to_lower() 01" );
    }
    {
        std::string result = remove( "", ' ' );
        test_suite.test_equality( result, std::string( "" ), "remove() 01" );
    }
    {
        std::string result = remove( "P_2_1_2_1_2_1", '_' );
        test_suite.test_equality( result, std::string( "P212121" ), "remove() 02" );
    }
    {
        std::string result = interlace( "", '_' );
        test_suite.test_equality( result, std::string( "" ), "interlace() 01" );
    }
    {
        std::string result = interlace( "P212121", '_' );
        test_suite.test_equality( result, std::string( "P_2_1_2_1_2_1" ), "interlace() 02" );
    }
    
    {
        std::string result = remove_delimiters( "", "", "" );
        test_suite.test_equality( result, std::string( "" ), "remove_delimiters() 01" );
    }
    {
        std::string result = remove_delimiters( "123", "12", "23" );
        test_suite.test_equality( result, std::string( "123" ), "remove_delimiters() 02" );
    }
    {
        std::string result = remove_delimiters( "1223", "12", "23" );
        test_suite.test_equality( result, std::string( "" ), "remove_delimiters() 03" );
    }
    {
        std::string result = remove_delimiters( "12r23", "12", "23" );
        test_suite.test_equality( result, std::string( "r" ), "remove_delimiters() 04" );
    }
    {
        std::string result = remove_delimiters( "12r23", "jj", "23" );
        test_suite.test_equality( result, std::string( "12r23" ), "remove_delimiters() 05" );
    }
    {
        std::string result = remove_delimiters( "12r23", "12", "jj" );
        test_suite.test_equality( result, std::string( "12r23" ), "remove_delimiters() 06" );
    }

    {
        size_t result = count_characters( "", '_' );
        test_suite.test_equality( result, static_cast<size_t>( 0 ), "count_characters() 01" );
    }
    {
        size_t result = count_characters( "_z__Tt_", '_' );
        test_suite.test_equality( result, static_cast<size_t>( 4 ), "count_characters() 02" );
    }
    {
        std::string result = centre( "1", 5, '=' );
        test_suite.test_equality( result, std::string( "==1==" ), "centre() 01" );
    }

// Splits a line into individual words, currently the separator is hard-coded to be a space or a tab
// "one word" and 'one word' are recognised as one word, the quotes are stripped.
// Empty words (either multiple spaces or e.g. "") are not retained.
// The following examples are not allowed:
//     one two th"ree
//     one two th'ree
//     one "two three
// @@ We need to allow escaped quotes, e.g. "He said \"yes\"."
// This kind of configurability suggests that this could be a class.
// std::vector< std::string > split( const std::string & input );

}

