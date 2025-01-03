#ifndef TESTSUITE_H
#define TESTSUITE_H

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

#include "BasicMathsFunctions.h"
#include "Complex.h"

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

class TestSuite
{
public:

    TestSuite() {}

    template< class T >
    void test_equality( const T & lhs, const T & rhs, const std::string & error_message )
    {
        if ( ! ( lhs == rhs ) )
            error_messages_.push_back( error_message );
    }

    void test_equality( const size_t lhs, const int rhs, const std::string & error_message )
    {
        if ( rhs < 0 )
            throw std::runtime_error( "test_equality( size_t, int ): int is negative in " + error_message );
        if ( ! ( lhs == static_cast< size_t >( rhs ) ) )
            error_messages_.push_back( error_message );
    }

    void test_equality_double( const double lhs, const double rhs, const std::string & error_message, double tolerance = TOLERANCE )
    {
        if ( absolute( lhs - rhs ) > tolerance )
            error_messages_.push_back( error_message );
    }

    void test_equality_Complex( const Complex lhs, const Complex rhs, const std::string & error_message, double tolerance = TOLERANCE )
    {
        if ( ( absolute( lhs.real() - rhs.real() ) > tolerance ) || ( absolute( lhs.imaginary() - rhs.imaginary() ) > tolerance ) )
            error_messages_.push_back( error_message );
    }

    void log_error( const std::string & error_message )
    {
        error_messages_.push_back( error_message );
    }

    void report() const;

private:
    std::vector< std::string > error_messages_;
};

#endif // TESTSUITE_H

