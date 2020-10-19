/* *********************************************
Copyright (c) 2013-2020, Cornelis Jan (Jacco) van de Streek
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

#include "Fraction.h"
#include "TestSuite.h"

#include <cmath>
#include <iostream>

namespace
{

    void test_one_fraction( TestSuite & test_suite,
                            const Fraction & fraction,
                            const int integer_part,
                            const int numerator,
                            const int denominator,
                            const std::string & error_message )
    {
//        std::cout << fraction.to_string() << std::endl;
        test_suite.test_equality( fraction.integer_part(), integer_part, error_message + " (integer_part)" );
        test_suite.test_equality( fraction.numerator(),    numerator,    error_message + " (numerator)" );
        test_suite.test_equality( fraction.denominator(),  denominator,  error_message + " (denominator)" );
    }

    int round_to_int( double x )
    {
        return ( x < 0 ) ? static_cast<int>( x - 0.5 ) : static_cast<int>( x + 0.5 );
    }

}

void test_fraction( TestSuite & test_suite )
{
    std::cout << "Now running tests for Fraction." << std::endl;
// Constructors
{
    Fraction fraction;
    test_one_fraction( test_suite, fraction, 0, 0, 1, "Fraction::Fraction() 01: error" );
}
{
    Fraction fraction( 0 );
    test_one_fraction( test_suite, fraction, 0, 0, 1, "Fraction::Fraction(int) 01: error" );
}
{
    Fraction fraction( 8 );
    test_one_fraction( test_suite, fraction, 8, 0, 1, "Fraction::Fraction(int) 02: error" );
}
{
    Fraction fraction( -8 );
    test_one_fraction( test_suite, fraction, -8, 0, 1, "Fraction::Fraction(int) 03: error" );
}
{
    try
    {
    Fraction fraction( 3, 0 );
    test_suite.log_error( "Fraction::Fraction() should have thrown 01" );
    }
    catch ( std::exception & e )
    {
    }
}
{
    try
    {
    Fraction fraction( 1, 0 );
    test_suite.log_error( "Fraction::Fraction() should have thrown 02" );
    }
    catch ( std::exception & e )
    {
    }
}
{
    Fraction fraction( 1, 4 );
    test_one_fraction( test_suite, fraction, 0, 1, 4, "Fraction::Fraction(int,int) 01: error" );
}
{
    Fraction fraction( 4, 6 );
    test_one_fraction( test_suite, fraction, 0, 2, 3, "Fraction::Fraction(int,int) 02: error" );
}
{
    Fraction fraction( 6, 4 );
    test_one_fraction( test_suite, fraction, 1, 1, 2, "Fraction::Fraction(int,int) 03: error" );
}
{
    Fraction fraction( 7, 7 );
    test_one_fraction( test_suite, fraction, 1, 0, 1, "Fraction::Fraction(int,int) 04: error" );
}
{
    Fraction fraction( -1, 4 );
    test_one_fraction( test_suite, fraction, 0, -1, 4, "Fraction::Fraction(int,int) 05: error" );
}
{
    Fraction fraction( -1, -4 );
    test_one_fraction( test_suite, fraction, 0, 1, 4, "Fraction::Fraction(int,int) 06: error" );
}
{
    Fraction fraction( 1, -4 );
    test_one_fraction( test_suite, fraction, 0, -1, 4, "Fraction::Fraction(int,int) 07: error" );
}
{
    Fraction fraction( -4, 6 );
    test_one_fraction( test_suite, fraction, 0, -2, 3, "Fraction::Fraction(int,int) 08: error" );
}
{
    Fraction fraction( 4, -6 );
    test_one_fraction( test_suite, fraction, 0, -2, 3, "Fraction::Fraction(int,int) 09: error" );
}
{
    Fraction fraction( -4, -6 );
    test_one_fraction( test_suite, fraction, 0, 2, 3, "Fraction::Fraction(int,int) 10: error" );
}
{
    Fraction fraction( -6, 4 );
    test_one_fraction( test_suite, fraction, -1, -1, 2, "Fraction::Fraction(int,int) 11: error" );
}
{
    Fraction fraction( 6, -4 );
    test_one_fraction( test_suite, fraction, -1, -1, 2, "Fraction::Fraction(int,int) 12: error" );
}
{
    Fraction fraction( -6, -4 );
    test_one_fraction( test_suite, fraction, 1, 1, 2, "Fraction::Fraction(int,int) 13: error" );
}
{
    Fraction fraction( -7, 7 );
    test_one_fraction( test_suite, fraction, -1, 0, 1, "Fraction::Fraction(int,int) 14: error" );
}
{
    Fraction fraction( 7, -7 );
    test_one_fraction( test_suite, fraction, -1, 0, 1, "Fraction::Fraction(int,int) 15: error" );
}
{
    Fraction fraction( -7, -7 );
    test_one_fraction( test_suite, fraction, 1, 0, 1, "Fraction::Fraction(int,int) 16: error" );
}
{
    Fraction fraction( 8, 2, 3 );
    test_one_fraction( test_suite, fraction, 8, 2, 3, "Fraction::Fraction(int,int,int) 01: error" );
}
{
    Fraction fraction( 8, 6, 4 );
    test_one_fraction( test_suite, fraction, 9, 1, 2, "Fraction::Fraction(int,int,int) 02: error" );
}
{
    Fraction fraction( 8, 7, 7 );
    test_one_fraction( test_suite, fraction, 9, 0, 1, "Fraction::Fraction(int,int,int) 03: error" );
}
{
    Fraction fraction( -8, 2, 3 );
    test_one_fraction( test_suite, fraction, -7, -1, 3, "Fraction::Fraction(int,int,int) 04: error" );
}
{
    Fraction fraction( -8, -2, 3 );
    test_one_fraction( test_suite, fraction, -8, -2, 3, "Fraction::Fraction(int,int,int) 05: error" );
}
{
    Fraction fraction( -8, 2, -3 );
    test_one_fraction( test_suite, fraction, -8, -2, 3, "Fraction::Fraction(int,int,int) 06: error" );
}
{
    Fraction fraction( -8, -2, -3 );
    test_one_fraction( test_suite, fraction, -7, -1, 3, "Fraction::Fraction(int,int,int) 07: error" );
}
{
    Fraction fraction( -8, -10, 8 );
    test_one_fraction( test_suite, fraction, -9, -1, 4, "Fraction::Fraction(int,int,int) 08: error" );
}
{
    Fraction fraction( -8, 10, -8 );
    test_one_fraction( test_suite, fraction, -9, -1, 4, "Fraction::Fraction(int,int,int) 09: error" );
}
{
    Fraction fraction( -8, 10, 8 );
    test_one_fraction( test_suite, fraction, -6, -3, 4, "Fraction::Fraction(int,int,int) 10: error" );
}
{
    Fraction fraction( -8, -10, -8 );
    test_one_fraction( test_suite, fraction, -6, -3, 4, "Fraction::Fraction(int,int,int) 11: error" );
}
{
    Fraction fraction( 8, -10, 8 );
    test_one_fraction( test_suite, fraction, 6, 3, 4, "Fraction::Fraction(int,int,int) 12: error" );
}
{
    Fraction fraction( 8, 10, -8 );
    test_one_fraction( test_suite, fraction, 6, 3, 4, "Fraction::Fraction(int,int,int) 13: error" );
}
{
    Fraction fraction( 8, 6, 4 );
    test_one_fraction( test_suite, fraction, 9, 1, 2, "Fraction::Fraction(int,int,int) 14: error" );
}
{
    Fraction fraction( 8, -6, -4 );
    test_one_fraction( test_suite, fraction, 9, 1, 2, "Fraction::Fraction(int,int,int) 15: error" );
}
{
    Fraction fraction( 8, 7, 7 );
    test_one_fraction( test_suite, fraction, 9, 0, 1, "Fraction::Fraction(int,int,int) 16: error" );
}
// to_string()
{
    Fraction fraction( 0, 0, 1 );
    test_suite.test_equality( fraction.to_string(), std::string("0"), "Fraction::to_string(): test 01" );
}
{
    Fraction fraction( 1, 3 );
    test_suite.test_equality( fraction.to_string(), std::string("1/3"), "Fraction::to_string(): test 02" );
}
{
    Fraction fraction( -7, -1, 3 );
    test_suite.test_equality( fraction.to_string(), std::string("-7 + -1/3"), "Fraction::to_string(): test 03" );
}
// to_double();
{
    Fraction fraction( 7 );
    test_suite.test_equality( fraction.to_double(), 7.0, "Fraction::to_double(): test 01" );
}
{
    Fraction fraction( -1, 4 );
    test_suite.test_equality( fraction.to_double(), -0.25, "Fraction::to_double(): test 02" );
}
{
    Fraction fraction( 7, 1, 4 );
    test_suite.test_equality( fraction.to_double(), 7.25, "Fraction::to_double(): test 03" );
}
{
    Fraction fraction( 1, 3 );
    test_suite.test_equality( fraction.to_double(), 1.0/3.0, "Fraction::to_double(): test 04" );
}
// absolute()
{
    Fraction fraction( -7, -1, 3 );
    fraction.absolute();
    test_one_fraction( test_suite, fraction, 7, 1, 3, "Fraction::absolute() 01: error" );
}
{
    Fraction fraction( 7, 1, 3 );
    fraction.absolute();
    test_one_fraction( test_suite, fraction, 7, 1, 3, "Fraction::absolute() 02: error" );
}
{
    Fraction fraction( 0, -1, 3 );
    fraction.absolute();
    test_one_fraction( test_suite, fraction, 0, 1, 3, "Fraction::absolute() 03: error" );
}
{
    Fraction fraction( 0, 1, 3 );
    fraction.absolute();
    test_one_fraction( test_suite, fraction, 0, 1, 3, "Fraction::absolute() 04: error" );
}
{
    Fraction fraction( -7, 0, 3 );
    fraction.absolute();
    test_one_fraction( test_suite, fraction, 7, 0, 1, "Fraction::absolute() 05: error" );
}
{
    Fraction fraction( 7, 0, 3 );
    fraction.absolute();
    test_one_fraction( test_suite, fraction, 7, 0, 1, "Fraction::absolute() 06: error" );
}
//    void reciprocal();
{
    Fraction fraction( 7 );
    fraction.reciprocal();
    test_one_fraction( test_suite, fraction, 0, 1, 7, "Fraction::reciprocal() 01: error" );
}
{
    Fraction fraction( -7 );
    fraction.reciprocal();
    test_one_fraction( test_suite, fraction, 0, -1, 7, "Fraction::reciprocal() 02: error" );
}
{
    Fraction fraction( 2, 3 );
    fraction.reciprocal();
    test_one_fraction( test_suite, fraction, 1, 1, 2, "Fraction::reciprocal() 03: error" );
}
{
    Fraction fraction( -2, 3 );
    fraction.reciprocal();
    test_one_fraction( test_suite, fraction, -1, -1, 2, "Fraction::reciprocal() 04: error" );
}
{
    try
    {
    Fraction fraction( 0, 1 );
    fraction.reciprocal();
    test_suite.log_error( "Fraction::reciprocal() should have thrown" );
    }
    catch ( std::exception & e )
    {
    }
}
// square()
{
    Fraction fraction(0,0,1);
    fraction.square();
    test_one_fraction( test_suite, fraction, 0, 0, 1, "Fraction::square() 01: error" );
}
{
    Fraction fraction(1,0,1);
    fraction.square();
    test_one_fraction( test_suite, fraction, 1, 0, 1, "Fraction::square() 02: error" );
}
{
    Fraction fraction(-2,-1,7);
    fraction.square();
    test_one_fraction( test_suite, fraction, 4, 29, 49, "Fraction::square() 03: error" );
}
{
    Fraction fraction(5,0,1);
    fraction.square();
    test_one_fraction( test_suite, fraction, 25, 0, 1, "Fraction::square() 04: error" );
}
{
    Fraction fraction(0,1,2);
    fraction.square();
    test_one_fraction( test_suite, fraction, 0, 1, 4, "Fraction::square() 05: error" );
}
{
    Fraction fraction(-5,-1,7);
    fraction.square();
    test_one_fraction( test_suite, fraction, 26, 22, 49, "Fraction::square() 06: error" );
}
// power( const int n )
{
    Fraction fraction(0,0,1);
    fraction.power(0);
    test_one_fraction( test_suite, fraction, 1, 0, 1, "Fraction::power() 01: error" );
}
{
    Fraction fraction(0,0,1);
    fraction.power(1);
    test_one_fraction( test_suite, fraction, 0, 0, 1, "Fraction::power() 02: error" );
}
{
    Fraction fraction(-5,-1,7);
    fraction.power(0);
    test_one_fraction( test_suite, fraction, 1, 0, 1, "Fraction::power() 03: error" );
}
{
    Fraction fraction(-5,-1,7);
    fraction.power(1);
    test_one_fraction( test_suite, fraction, -5, -1, 7, "Fraction::power() 04: error" );
}
{
    Fraction fraction(-5,-1,7);
    fraction.power(2);
    test_one_fraction( test_suite, fraction, 26, 22, 49, "Fraction::power() 05: error" );
}
{
    Fraction fraction(-5,-1,7);
    fraction.power(3);
    test_one_fraction( test_suite, fraction, -136, -8, 343, "Fraction::power() 06: error" );
}
{
    Fraction fraction(-5,-1,7);
    fraction.power(-3);
    test_one_fraction( test_suite, fraction, 0, -343, 46656, "Fraction::power() 07: error" );
}
{
    try
    {
    Fraction fraction( 0, 1 );
    fraction.power(-2);
    test_suite.log_error( "Fraction::power() should have thrown" );
    }
    catch ( std::exception & e )
    {
    }
}
// Fraction operator+( const Fraction & rhs ) const
{
    Fraction lhs( 7, 2, 3 );
    Fraction rhs( 2, 1, 4 );
    Fraction result = lhs + rhs;
    test_one_fraction( test_suite, result, 9, 11, 12, "Fraction::operator+(const Fraction&) 01: error" );
}
{
    Fraction lhs( 7, 2, 3 );
    Fraction rhs( -7, -2, 3 );
    Fraction result = lhs + rhs;
    test_one_fraction( test_suite, result, 0, 0, 1, "Fraction::operator+(const Fraction&) 02: error" );
}
{
    Fraction lhs( -7, -1, 3 );
    Fraction rhs( -5, -1, 6 );
    Fraction result = lhs + rhs;
    test_one_fraction( test_suite, result, -12, -1, 2, "Fraction::operator+(const Fraction&) 03: error" );
}
{
    Fraction lhs( -7, -2, 3 );
    Fraction rhs( 2, 1, 4 );
    Fraction result = lhs + rhs;
    test_one_fraction( test_suite, result, -5, -5, 12, "Fraction::operator+(const Fraction&) 04: error" );
}
// Fraction operator-( const Fraction & rhs ) const
{
    Fraction lhs( 2, 2, 3 );
    Fraction rhs( 7, 1, 4 );
    Fraction result = lhs - rhs;
    test_one_fraction( test_suite, result, -4, -7, 12, "Fraction::operator-(const Fraction&) 01: error" );
}
{
    Fraction lhs( 7, 2, 3 );
    Fraction rhs( 7, 2, 3 );
    Fraction result = lhs - rhs;
    test_one_fraction( test_suite, result, 0, 0, 1, "Fraction::operator-(const Fraction&) 02: error" );
}
{
    Fraction lhs( -7, -1, 3 );
    Fraction rhs( 5, 1, 6 );
    Fraction result = lhs - rhs;
    test_one_fraction( test_suite, result, -12, -1, 2, "Fraction::operator-(const Fraction&) 03: error" );
}
{
    Fraction lhs( -7, -2, 3 );
    Fraction rhs( 2, 1, 4 );
    Fraction result = lhs - rhs;
    test_one_fraction( test_suite, result, -9, -11, 12, "Fraction::operator-(const Fraction&) 04: error" );
}
// Fraction operator*( const Fraction& rhs ) const
{
    Fraction lhs( -7, -2, 3 );
    Fraction rhs( 2, 1, 4 );
    Fraction result = lhs * rhs;
    test_one_fraction( test_suite, result, -17, -1, 4, "Fraction::operator*(const Fraction&) 01: error" );
}
{
    Fraction lhs( -7, -2, 3 );
    Fraction rhs( -2, -1, 4 );
    Fraction result = lhs * rhs;
    test_one_fraction( test_suite, result, 17, 1, 4, "Fraction::operator*(const Fraction&) 02: error" );
}
{
    Fraction lhs( -7, -2, 3 );
    Fraction rhs( 1, 0, 1 );
    Fraction result = lhs * rhs;
    test_one_fraction( test_suite, result, -7, -2, 3, "Fraction::operator*(const Fraction&) 03: error" );
}
{
    Fraction lhs( -7, -2, 3 );
    Fraction rhs( 0, 0, 1 );
    Fraction result = lhs * rhs;
    test_one_fraction( test_suite, result, 0, 0, 1, "Fraction::operator*(const Fraction&) 04: error" );
}
// Fraction operator/( const Fraction& rhs ) const
{
    Fraction lhs( -7, -2, 3 );
    Fraction rhs( 2, 1, 4 );
    Fraction result = lhs / rhs;
    test_one_fraction( test_suite, result, -3, -11, 27, "Fraction::operator/(const Fraction&) 01: error" );
}
{
    Fraction lhs( -7, -2, 3 );
    Fraction rhs( -2, -1, 4 );
    Fraction result = lhs / rhs;
    test_one_fraction( test_suite, result, 3, 11, 27, "Fraction::operator/(const Fraction&) 02: error" );
}
{
    Fraction lhs( -7, -2, 3 );
    Fraction rhs( 0 );
    try
    {
    /*Fraction result =*/ lhs / rhs;
    test_suite.log_error( "Fraction::operator/(const Fraction&) should have thrown" );
    }
    catch ( std::exception & e )
    {
    }
}
{
    Fraction lhs( -7, -2, 3 );
    Fraction rhs( 1, 0, 1 );
    Fraction result = lhs / rhs;
    test_one_fraction( test_suite, result, -7, -2, 3, "Fraction::operator/(const Fraction&) 03: error" );
}
// Fraction operator+() const
{
    Fraction fraction( 2, 3 );
    fraction = +fraction;
    test_one_fraction( test_suite, fraction, 0, 2, 3, "Fraction::operator+() 01: error" );
}
{
    Fraction fraction( -2, 3 );
    fraction = +fraction;
    test_one_fraction( test_suite, fraction, 0, -2, 3, "Fraction::operator+() 02: error" );
}
// Fraction operator-() const
{
    Fraction fraction( 2, 3 );
    fraction = -fraction;
    test_one_fraction( test_suite, fraction, 0, -2, 3, "Fraction::operator-() 01: error" );
}
{
    Fraction fraction( -2, 3 );
    fraction = -fraction;
    test_one_fraction( test_suite, fraction, 0, 2, 3, "Fraction::operator-() 02: error" );
}
// Fraction& operator+=( const Fraction & rhs);
{
    Fraction lhs( 7, 2, 3 );
    Fraction rhs( 2, 1, 4 );
    lhs += rhs;
    test_one_fraction( test_suite, lhs, 9, 11, 12, "Fraction::operator+=(const Fraction&) 01: error" );
}
{
    Fraction lhs( 7, 2, 3 );
    Fraction rhs( -7, -2, 3 );
    lhs += rhs;
    test_one_fraction( test_suite, lhs, 0, 0, 1, "Fraction::operator+=(const Fraction&) 02: error" );
}
{
    Fraction lhs( -7, -1, 3 );
    Fraction rhs( -5, -1, 6 );
    lhs += rhs;
    test_one_fraction( test_suite, lhs, -12, -1, 2, "Fraction::operator+=(const Fraction&) 03: error" );
}
{
    Fraction lhs( -7, -2, 3 );
    Fraction rhs( 2, 1, 4 );
    lhs += rhs;
    test_one_fraction( test_suite, lhs, -5, -5, 12, "Fraction::operator+=(const Fraction&) 04: error" );
}
// Fraction& operator-=( const Fraction & rhs);
{
    Fraction lhs( 2, 2, 3 );
    Fraction rhs( 7, 1, 4 );
    lhs -= rhs;
    test_one_fraction( test_suite, lhs, -4, -7, 12, "Fraction::operator-=(const Fraction&) 01: error" );
}
{
    Fraction lhs( 7, 2, 3 );
    Fraction rhs( 7, 2, 3 );
    lhs -= rhs;
    test_one_fraction( test_suite, lhs, 0, 0, 1, "Fraction::operator-=(const Fraction&) 02: error" );
}
{
    Fraction lhs( -7, -1, 3 );
    Fraction rhs( 5, 1, 6 );
    lhs -= rhs;
    test_one_fraction( test_suite, lhs, -12, -1, 2, "Fraction::operator-=(const Fraction&) 03: error" );
}
{
    Fraction lhs( -7, -2, 3 );
    Fraction rhs( 2, 1, 4 );
    lhs -= rhs;
    test_one_fraction( test_suite, lhs, -9, -11, 12, "Fraction::operator-=(const Fraction&) 04: error" );
}
// Fraction& operator*=( const Fraction & rhs);
{
    Fraction lhs( -7, -2, 3 );
    Fraction rhs( 2, 1, 4 );
    lhs *= rhs;
    test_one_fraction( test_suite, lhs, -17, -1, 4, "Fraction::operator*=(const Fraction&) 01: error" );
}
{
    Fraction lhs( -7, -2, 3 );
    Fraction rhs( -2, -1, 4 );
    lhs *= rhs;
    test_one_fraction( test_suite, lhs, 17, 1, 4, "Fraction::operator*=(const Fraction&) 02: error" );
}
// Fraction& operator/=( const Fraction & rhs);
{
    Fraction lhs( -7, -2, 3 );
    Fraction rhs( 2, 1, 4 );
    lhs /= rhs;
    test_one_fraction( test_suite, lhs, -3, -11, 27, "Fraction::operator/=(const Fraction&) 01: error" );
}
{
    Fraction lhs( -7, -2, 3 );
    Fraction rhs( -2, -1, 4 );
    lhs /= rhs;
    test_one_fraction( test_suite, lhs, 3, 11, 27, "Fraction::operator/=(const Fraction&) 02: error" );
}
{
    Fraction lhs( -7, -2, 3 );
    Fraction rhs( 0 );
    try
    {
    lhs /= rhs;
    test_suite.log_error( "Fraction::operator/=(const Fraction&) should have thrown" );
    }
    catch ( std::exception & e )
    {
    }
}
// Fraction& operator++();    // Prefix
{
    Fraction fraction( -2, 3 );
    Fraction fraction2 = ++fraction;
    test_one_fraction( test_suite, fraction, 0, 1, 3, "Fraction::operator++() 01a: error" );
    test_one_fraction( test_suite, fraction2, 0, 1, 3, "Fraction::operator++() 01b: error" );
}
{
    Fraction fraction( 7, 2, 3 );
    ++fraction;
    test_one_fraction( test_suite, fraction, 8, 2, 3, "Fraction::operator++() 02: error" );
}
{
    Fraction fraction( -4, -2, 3 );
    ++fraction;
    test_one_fraction( test_suite, fraction, -3, -2, 3, "Fraction::operator++() 03: error" );
}
{
    Fraction fraction( -1, -2, 3 );
    ++fraction;
    test_one_fraction( test_suite, fraction, 0, -2, 3, "Fraction::operator++() 04: error" );
}
// Fraction  operator++(int); // Postfix
{
    Fraction fraction( -2, 3 );
    Fraction fraction2 = fraction++;
    test_one_fraction( test_suite, fraction, 0, 1, 3, "Fraction::operator++(int) 01a: error" );
    test_one_fraction( test_suite, fraction2, 0, -2, 3, "Fraction::operator++(int) 01b: error" );
}
{
    Fraction fraction( 7, 2, 3 );
    fraction++;
    test_one_fraction( test_suite, fraction, 8, 2, 3, "Fraction::operator++(int) 02: error" );
}
{
    Fraction fraction( -4, -2, 3 );
    fraction++;
    test_one_fraction( test_suite, fraction, -3, -2, 3, "Fraction::operator++(int) 03: error" );
}
{
    Fraction fraction( -1, -2, 3 );
    fraction++;
    test_one_fraction( test_suite, fraction, 0, -2, 3, "Fraction::operator++(int) 04: error" );
}
// Fraction& operator--();    // Prefix
    {
    Fraction fraction( 2, 3 );
    Fraction fraction2 = --fraction;
    test_one_fraction( test_suite, fraction, 0, -1, 3, "Fraction::operator--() 01a: error" );
    test_one_fraction( test_suite, fraction2, 0, -1, 3, "Fraction::operator--() 01b: error" );
}
{
    Fraction fraction( 7, 2, 3 );
    --fraction;
    test_one_fraction( test_suite, fraction, 6, 2, 3, "Fraction::operator--() 02: error" );
}
{
    Fraction fraction( -4, -2, 3 );
    --fraction;
    test_one_fraction( test_suite, fraction, -5, -2, 3, "Fraction::operator--() 03: error" );
}
{
    Fraction fraction( 1, 2, 3 );
    --fraction;
    test_one_fraction( test_suite, fraction, 0, 2, 3, "Fraction::operator--() 04: error" );
}
// Fraction  operator--(int); // Postfix
{
    Fraction fraction( 2, 3 );
    Fraction fraction2 = fraction--;
    test_one_fraction( test_suite, fraction, 0, -1, 3, "Fraction::operator--(int) 01a: error" );
    test_one_fraction( test_suite, fraction2, 0, 2, 3, "Fraction::operator--(int) 01b: error" );
}
{
    Fraction fraction( 7, 2, 3 );
    fraction--;
    test_one_fraction( test_suite, fraction, 6, 2, 3, "Fraction::operator--(int) 02: error" );
}
{
    Fraction fraction( -4, -2, 3 );
    fraction--;
    test_one_fraction( test_suite, fraction, -5, -2, 3, "Fraction::operator--(int) 03: error" );
}
{
    Fraction fraction( 1, 2, 3 );
    fraction--;
    test_one_fraction( test_suite, fraction, 0, 2, 3, "Fraction::operator--(int) 04: error" );
}
// bool operator==( const Fraction & rhs ) const
{
    Fraction fraction1( 8, 6, -4 );
    Fraction fraction2( 6, 1, 2 );
    test_suite.test_equality( fraction1 == fraction2, true , "Fraction::operator== 01: error" );
}
{
    Fraction fraction1( 1, 100 );
    Fraction fraction2( 1, 3 );
    test_suite.test_equality( fraction1 == fraction2, false , "Fraction::operator== 02: error" );
}
// bool operator!=( const Fraction& rhs ) const
{
    Fraction fraction1( 8, 6, -4 );
    Fraction fraction2( 6, 1, 2 );
    test_suite.test_equality( fraction1 != fraction2, false , "Fraction::operator!= 01: error" );
}
{
    Fraction fraction1( 1, 100 );
    Fraction fraction2( 1, 3 );
    test_suite.test_equality( fraction1 != fraction2, true , "Fraction::operator!= 02: error" );
}
// bool operator< ( const Fraction& rhs ) const
{
    Fraction fraction1( 1, 100 );
    Fraction fraction2( 1, 3 );
    test_suite.test_equality( fraction1 < fraction2, true , "Fraction::operator< 01: error" );
}
{
    Fraction fraction1( -1, 3 );
    Fraction fraction2( 1, 100 );
    test_suite.test_equality( fraction1 < fraction2, true , "Fraction::operator< 02: error" );
}
{
    Fraction fraction1( -1, -1, 100 );
    Fraction fraction2( -1, 3 );
    test_suite.test_equality( fraction1 < fraction2, true , "Fraction::operator< 03: error" );
}
{
    Fraction fraction1( 1, 3 );
    Fraction fraction2( 1, 3 );
    test_suite.test_equality( fraction1 < fraction2, false , "Fraction::operator< 04: error" );
}
{
    Fraction fraction1( 1, 3 );
    Fraction fraction2( 1, 4 );
    test_suite.test_equality( fraction1 < fraction2, false , "Fraction::operator< 05: error" );
}
// bool operator> ( const Fraction& rhs ) const
{
    Fraction fraction1( 1, 100 );
    Fraction fraction2( 1, 3 );
    test_suite.test_equality( fraction1 > fraction2, false , "Fraction::operator> 01: error" );
}
{
    Fraction fraction1( -1, 3 );
    Fraction fraction2( 1, 100 );
    test_suite.test_equality( fraction1 > fraction2, false , "Fraction::operator> 02: error" );
}
{
    Fraction fraction1( -1, -1, 100 );
    Fraction fraction2( -1, 3 );
    test_suite.test_equality( fraction1 > fraction2, false , "Fraction::operator> 03: error" );
}
{
    Fraction fraction1( 1, 3 );
    Fraction fraction2( 1, 3 );
    test_suite.test_equality( fraction1 > fraction2, false , "Fraction::operator> 04: error" );
}
{
    Fraction fraction1( 1, 3 );
    Fraction fraction2( 1, 4 );
    test_suite.test_equality( fraction1 > fraction2, true , "Fraction::operator> 05: error" );
}
// bool operator>=( const Fraction& rhs ) const
{
    Fraction fraction1( 1, 100 );
    Fraction fraction2( 1, 3 );
    test_suite.test_equality( fraction1 >= fraction2, false , "Fraction::operator>= 01: error" );
}
{
    Fraction fraction1( -1, 3 );
    Fraction fraction2( 1, 100 );
    test_suite.test_equality( fraction1 >= fraction2, false , "Fraction::operator>= 02: error" );
}
{
    Fraction fraction1( -1, -1, 100 );
    Fraction fraction2( -1, 3 );
    test_suite.test_equality( fraction1 >= fraction2, false , "Fraction::operator>= 03: error" );
}
{
    Fraction fraction1( 1, 3 );
    Fraction fraction2( 1, 3 );
    test_suite.test_equality( fraction1 >= fraction2, true , "Fraction::operator>= 04: error" );
}
{
    Fraction fraction1( 1, 3 );
    Fraction fraction2( 1, 4 );
    test_suite.test_equality( fraction1 >= fraction2, true , "Fraction::operator>= 05: error" );
}
// bool operator<=( const Fraction& rhs ) const
{
    Fraction fraction1( 1, 100 );
    Fraction fraction2( 1, 3 );
    test_suite.test_equality( fraction1 <= fraction2, true , "Fraction::operator<= 01: error" );
}
{
    Fraction fraction1( -1, 3 );
    Fraction fraction2( 1, 100 );
    test_suite.test_equality( fraction1 <= fraction2, true , "Fraction::operator<= 02: error" );
}
{
    Fraction fraction1( -1, -1, 100 );
    Fraction fraction2( -1, 3 );
    test_suite.test_equality( fraction1 <= fraction2, true , "Fraction::operator<= 03: error" );
}
{
    Fraction fraction1( 1, 3 );
    Fraction fraction2( 1, 3 );
    test_suite.test_equality( fraction1 <= fraction2, true , "Fraction::operator<= 04: error" );
}
{
    Fraction fraction1( 1, 3 );
    Fraction fraction2( 1, 4 );
    test_suite.test_equality( fraction1 <= fraction2, false , "Fraction::operator<= 05: error" );
}
// double   operator*( const Fraction& lhs, const double    rhs )
{
    Fraction lhs( -1, -1, 2);
    double   rhs( 2.5 );
    double result = lhs * rhs;
    test_suite.test_equality( result, -3.75, "operator*(const Fraction&,const double) 01: error" );
}
// double   operator*( const double    lhs, const Fraction& rhs )
{
    Fraction rhs( -1, -1, 2);
    double   lhs( 2.5 );
    double result = lhs * rhs;
    test_suite.test_equality( result, -3.75, "operator*(const double,const Fraction&) 01: error" );
}
// Fraction operator*( const Fraction& lhs, const int       rhs )
{
    Fraction rhs( -1, -1, 2);
    int      lhs( 2 );
    Fraction result = lhs * rhs;
    test_one_fraction( test_suite, result, -3, 0, 1, "operator*(const Fraction&,const int) 01: error" );
}
{
    Fraction rhs( -1, -2, 3);
    int      lhs( 2 );
    Fraction result = lhs * rhs;
    test_one_fraction( test_suite, result, -3, -1, 3, "operator*(const Fraction&,const int) 02: error" );
}
// Fraction operator*( const int       lhs, const Fraction& rhs )
{
    Fraction lhs( -1, -1, 2);
    int      rhs( 2 );
    Fraction result = lhs * rhs;
    test_one_fraction( test_suite, result, -3, 0, 1, "operator*(const int,const Fraction&) 01: error" );
}
{
    Fraction lhs( -1, -2, 3);
    int      rhs( 2 );
    Fraction result = lhs * rhs;
    test_one_fraction( test_suite, result, -3, -1, 3, "operator*(const int,const Fraction&) 02: error" );
}
// double   operator+( const Fraction& lhs, const double    rhs )
{
    Fraction lhs( -1, -1, 4);
    double   rhs( 2.75 );
    double result = lhs + rhs;
    test_suite.test_equality( result, 1.5, "operator+(const Fraction&,const double) 01: error" );
}
// double   operator+( const double    lhs, const Fraction& rhs )
{
    Fraction rhs( -1, -1, 4);
    double   lhs( 2.75 );
    double result = lhs + rhs;
    test_suite.test_equality( result, 1.5, "operator+(const double,const Fraction&) 01: error" );
}
// Fraction operator+( const Fraction& lhs, const int       rhs )
{
    Fraction lhs(-1,-1,3);
    int      rhs(3);
    Fraction result = lhs + rhs;
    test_one_fraction( test_suite, result, 1, 2, 3, "operator+(const Fraction&,const int) 01: error" );
}
{
    Fraction lhs(-1,-1,3);
    int      rhs(-3);
    Fraction result = lhs + rhs;
    test_one_fraction( test_suite, result, -4, -1, 3, "operator+(const Fraction&,const int) 02: error" );
}
{
    Fraction lhs(-1,-1,3);
    int      rhs(0);
    Fraction result = lhs + rhs;
    test_one_fraction( test_suite, result, -1, -1, 3, "operator+(const Fraction&,const int) 03: error" );
}
{
    Fraction lhs(0,0,1);
    int      rhs(-3);
    Fraction result = lhs + rhs;
    test_one_fraction( test_suite, result, -3, 0, 1, "operator+(const Fraction&,const int) 04: error" );
}
// Fraction operator+( const int       lhs, const Fraction& rhs )
{
    int      lhs( 3 );
    Fraction rhs( -1, -1, 3 );
    Fraction result = lhs + rhs;
    test_one_fraction( test_suite, result, 1, 2, 3, "operator+(const int,const Fraction&) 01: error" );
}
{
    int      lhs( -3 );
    Fraction rhs( -1, -1, 3 );
    Fraction result = lhs + rhs;
    test_one_fraction( test_suite, result, -4, -1, 3, "operator+(const int,const Fraction&) 02: error" );
}
{
    int      lhs( 0 );
    Fraction rhs( -1, -1, 3 );
    Fraction result = lhs + rhs;
    test_one_fraction( test_suite, result, -1, -1, 3, "operator+(const int,const Fraction&) 03: error" );
}
{
    int      lhs( -3 );
    Fraction rhs( 0, 0, 1 );
    Fraction result = lhs + rhs;
    test_one_fraction( test_suite, result, -3, 0, 1, "operator+(const int,const Fraction&) 04: error" );
}
// double   operator/( const Fraction& lhs, const double    rhs )
{
    Fraction lhs( -1, -1, 2 );
    double   rhs( 2.3 );
    double result = lhs / rhs;
    test_suite.test_equality( result, -1.5 / 2.3, "operator/(const Fraction&,const double) 01: error" );
}
// double   operator/( const double    lhs, const Fraction& rhs )
{
    double   lhs( -1.5 );
    Fraction rhs( 2, 3, 10 );
    double result = lhs / rhs;
    test_suite.test_equality( result, -1.5 / 2.3, "operator/(const double,const Fraction&) 01: error" );
}
// Fraction operator/( const Fraction& lhs, const int       rhs )
{
    Fraction lhs( -2, -1, 3) ;
    int      rhs( 6 );
    Fraction result = lhs / rhs;
    test_one_fraction( test_suite, result, 0, -7, 18, "operator/(const Fraction&,const int) 01: error" );
}
{
    Fraction lhs( -7, -2, 3 );
    int rhs( 0 );
    try
    {
    /*Fraction result =*/ lhs / rhs;
    test_suite.log_error( "operator/(const Fraction&,const int) should have thrown" );
    }
    catch ( std::exception& e )
    {
    }
}
// Fraction operator/( const int       lhs, const Fraction& rhs )
{
    int      lhs( 6 );
    Fraction rhs( -2, -1, 3 );
    Fraction result = lhs / rhs;
    test_one_fraction( test_suite, result, -2, -4, 7, "operator/(const int,const Fraction&) 01: error" );
}
{
    int      lhs(6);
    Fraction rhs( 0, 0, 1 );
    try
    {
        /*Fraction result =*/ lhs / rhs;
        test_suite.log_error( "operator/(const int,const Fraction&) should have thrown" );
    }
    catch ( std::exception& e )
    {

    }
}
// double   operator-( const Fraction& lhs, const double    rhs )
{
    Fraction lhs( -1, -1, 2 );
    double   rhs(2.3);
    test_suite.test_equality( lhs - rhs, -1.5 - 2.3, "operator-(const Fraction&,const double) 01: error" );
}
// double   operator-( const double    lhs, const Fraction& rhs )
{
    double   lhs( -1.5 );
    Fraction rhs( 2, 3, 10 );
    test_suite.test_equality( lhs - rhs, -1.5 - 2.3, "operator-(const double,const Fraction&) 01: error" );
}
// Fraction operator-( const Fraction& lhs, const int       rhs )
{
    Fraction lhs(-7,-1,4);
    int      rhs(5);
    test_one_fraction( test_suite, lhs - rhs, -12, -1, 4, "operator-(const Fraction&,const int) 01: error" );
}
// Fraction operator-( const int       lhs, const Fraction& rhs )
{
    int      lhs(-5);
    Fraction rhs(-3,-1,4);
    test_one_fraction( test_suite, lhs - rhs, -1, -3, 4, "operator-(const int,const Fraction&) 01: error" );
}
// Fraction from_floating_point( const double theDouble, const Fraction& smallestUnit )
{
    double original_value( 0.333 );
    Fraction fraction( double2fraction( original_value, Fraction( 1, 24 ) ) );
    test_one_fraction( test_suite, fraction, 0, 1, 3, "Fraction from_floating_point() 01: error" );
}
{
    double original_value( 0.333 );
    Fraction fraction( double2fraction( original_value, Fraction( 1, 10 ) ) );
    test_one_fraction( test_suite, fraction, 0, 3, 10, "Fraction from_floating_point() 02a: error" );
    test_suite.test_equality_double( fraction.to_double() - original_value, -0.033, "Fraction from_floating_point() 02b: error" );
}
{
    double maximum_error( 0.0 );
    int worst_i = -1;
    for ( int i( 0 ); i <= 24; ++i )
    {
        // We generate 0/24, 1/24 etc. as doubles rounded to three decimal places:
        // "0.000", "0.333", "0.250", "0.167", ...
        // This is what we can expect in e.g. .res files
        // Then we determine the maximum error.
        double original_value( static_cast<double>(i)/24.0 );
        original_value = static_cast<double>( round_to_int( original_value * 1000.0 ) ) / 1000.0;
        Fraction fraction = double2fraction( original_value, Fraction( 1, 24 ) );
        double error_i = fabs( fraction.to_double() - original_value );
        if ( maximum_error < error_i )
        {
            worst_i = i;
            maximum_error = error_i;
        }
    }
    test_suite.test_equality_double( maximum_error, 1.0/3000.0, "Fraction from_floating_point() : maximum error has changed" );
}
{
    double initial_value( 0.333 );
    int best_i = 1;
    double smallest_error = fabs( initial_value );
    for ( int i( 1 ); i <= 1000; ++i )
    {
        Fraction trial_i = double2fraction( initial_value, Fraction( 1, i ) );
        double error_i = fabs( trial_i.to_double() - initial_value );
 //       if ( error_i <= 0.0005 )
 //           std::cout << "i = " << i << ", Fraction = " << trial_i.to_string() << ", Current error is :" << error_i << std::endl;
        if ( error_i < smallest_error )
        {
            best_i = i;
            smallest_error = error_i;
        }
    }
  //  std::cout << "Smallest error = " << smallest_error << std::endl;
  //  std::cout << "Best i = " << best_i << std::endl;
    Fraction best_fraction( 1, best_i );
}
{
    Fraction fraction = Farey( 0.605551, 30 ); // 3/5 is wrong, it should be 17/28
    test_one_fraction( test_suite, fraction, 0, 17, 28, "Fraction Farey() 01: error" );
}
{
    Fraction fraction = Farey( -23.605551, 30 ); // 3/5 is wrong, it should be 17/28
    test_one_fraction( test_suite, fraction, -23, -17, 28, "Fraction Farey() -01: error" );
}
{
    Fraction fraction = Farey( 0.5, 100 );
    test_one_fraction( test_suite, fraction, 0, 1, 2, "Fraction Farey() 02: error" );
}
{
    Fraction fraction = Farey( 1.0/3.0, 1000 );
    test_one_fraction( test_suite, fraction, 0, 1, 3, "Fraction Farey() 03: error" );
}
{
    try
    {
    /* Fraction fraction = */ Farey( 0.5, 0 );
    test_suite.log_error( "Farey() should have thrown 01" );
    }
    catch ( std::exception& e )
    {
    }
}

}

