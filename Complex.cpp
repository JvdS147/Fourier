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

#include "Complex.h"
#include "BasicMathsFunctions.h"

#include <stdexcept>
#include <iostream>

// ********************************************************************************

Complex::Complex( const double real, const double imaginary ):
real_(real),
imaginary_(imaginary)
{
}

// ********************************************************************************

bool Complex::nearly_zero( const double tolerance ) const
{
    return ( ::nearly_zero( real_, tolerance ) && ::nearly_zero( imaginary_, tolerance ) );
}

// ********************************************************************************

void Complex::reciprocal()
{
    if ( nearly_zero() )
        throw std::runtime_error( "Complex::reciprocal(): divide by zero." );
    *this = Complex( real_ / norm2(), -imaginary_ / norm2() );    
}

// ********************************************************************************

void Complex::square()
{
    *this = Complex( ::square( real_ ) - ::square( imaginary_ ), 2.0 * real_ * imaginary_ );
}

// ********************************************************************************

void Complex::power( const int n )
{
    // Note: 0^0 = 1
    if ( n == 0 )
    {
        real_ = 1.0;
        imaginary_ = 0.0;
        return;
    }
    if ( n < 0 )
        reciprocal();
    size_t abs_n = abs( n );
    if ( ( abs_n == 1 ) || nearly_zero() )
        return;
    if ( abs_n == 2 )
    {
        square();
        return;
    }
    // Many clever things can be done here, e.g. using the binary power method or recursion.
    // Note that we gain something from using the binary method because our square() function is optimised
    // for Complex objects and faster than simply x * x (with x a Complex).
    // The current implementation is extremely simple.
    Complex z( *this );
    for ( size_t i(1); i != abs_n; ++i )
        *this = *this * z;
}

// ********************************************************************************

double Complex::norm() const
{
    return sqrt( norm2() );
}

// ********************************************************************************

double Complex::norm2() const
{
    return ::square( real_ ) + ::square( imaginary_ );
}

// ********************************************************************************

void Complex::show() const
{
    std::cout << "Re = " << real_ << ", Im = " << imaginary_<< std::endl;
}

// ********************************************************************************

Complex & Complex::operator+=( const Complex & rhs )
{
    real_      += rhs.real_;
    imaginary_ += rhs.imaginary_;
    return *this;
}

// ********************************************************************************

Complex & Complex::operator-=( const Complex & rhs )
{
    real_      -= rhs.real_;
    imaginary_ -= rhs.imaginary_;
    return *this;
}

// ********************************************************************************

Complex & Complex::operator*=( const Complex & rhs )
{
  //  numerator_    = integer_part_*rhs.numerator_*denominator_ + rhs.integer_part_*numerator_*rhs.denominator_ + numerator_*rhs.numerator_;
  //  denominator_ *= rhs.denominator_;
  //  integer_part_ *= rhs.integer_part_;
    return *this;
}

// *************************Complex*******************************************************

Complex & Complex::operator/=( const Complex & rhs )
{
 //   numerator_   = integer_part_*denominator_*rhs.denominator_     + numerator_*rhs.denominator_;
 //   denominator_ = rhs.integer_part_*denominator_*rhs.denominator_ + rhs.numerator_*denominator_;
 //   integer_part_ = 0;
    return *this;
}

// ********************************************************************************

Complex & Complex::operator++()    // Prefix
{
    ++real_;
    return *this;
}

// ********************************************************************************

Complex  Complex::operator++(int) // Postfix
{
    Complex old( *this );
    ++(*this);
    return old;
}

// ********************************************************************************

Complex & Complex::operator--()    // Prefix
{
    --real_;
    return *this;
}

// ********************************************************************************

Complex  Complex::operator--(int) // Postfix
{
    Complex old( *this );
    --(*this);
    return old;
}

// ********************************************************************************

bool nearly_equal( const Complex & lhs, const Complex & rhs, const double tolerance )
{
    return ( nearly_equal( lhs.real(), rhs.real(), tolerance ) && nearly_equal( lhs.imaginary(), rhs.imaginary(), tolerance ) );
}

// ********************************************************************************

Complex operator+( const Complex & lhs, const Complex & rhs )
{
    return Complex( lhs.real() + rhs.real(), lhs.imaginary() + rhs.imaginary() );
}

Complex operator-( const Complex & lhs, const Complex & rhs )
{
    return Complex( lhs.real() - rhs.real(), lhs.imaginary() - rhs.imaginary() );
}

Complex operator*( const Complex & lhs, const Complex & rhs )
{
    return Complex( lhs.real() * rhs.real() - lhs.imaginary() * rhs.imaginary(), lhs.real() * rhs.imaginary() + rhs.real() * lhs.imaginary() );
}

Complex operator/( const Complex & lhs, Complex rhs )
{
    rhs.reciprocal();
    return lhs * rhs;
}

Complex operator+( const double lhs, const Complex & rhs )
{
    return Complex( lhs + rhs.real(), rhs.imaginary() );
}

Complex operator-( const double lhs, const Complex & rhs )
{
    return Complex( lhs - rhs.real(), -rhs.imaginary() );
}

Complex operator*( const double lhs, const Complex & rhs )
{
    return Complex( lhs * rhs.real(), lhs * rhs.imaginary() );
}

Complex operator/( const double lhs, Complex rhs )
{
    rhs.reciprocal();
    return lhs * rhs;
}

Complex operator+( const Complex & lhs, const double rhs )
{
    return Complex( lhs.real() + rhs, lhs.imaginary() );
}

Complex operator-( const Complex & lhs, const double rhs )
{
    return Complex( lhs.real() - rhs, lhs.imaginary() );
}

Complex operator*( const Complex & lhs, const double rhs )
{
    return Complex( rhs * lhs.real(), rhs * lhs.imaginary() );
}

Complex operator/( const Complex & lhs, const double rhs )
{
    if ( nearly_zero( rhs ) )
        throw std::runtime_error( "operator/( Complex, double ): divide by zero." );
    return Complex( lhs.real() / rhs, lhs.imaginary() / rhs );
}

// ********************************************************************************

Complex square( const Complex value )
{
    Complex result( value );
    result.square();
    return result;
}

// ********************************************************************************

Complex exponential( const Complex z )
{
    return Complex( exp( z.real() ) * cos( z.imaginary() ), exp( z.real() ) * sin( z.imaginary() ) );
}

// ********************************************************************************

