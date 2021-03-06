#ifndef COMPLEX_H
#define COMPLEX_H

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

/*
    Mainly a wrapper for std::complex< double > .
    
    A quick Google search suggests that std::complex< double > * double etc. (mixing built-ins with std::complex) does not work.
    Also, apparently std::complex< double >::norm() returns the square of the norm, not the norm.
*/
class Complex
{
public:

    explicit Complex( const double real = 0.0, const double imaginary = 0.0 );

    static Complex i() { return Complex( 0.0, 1.0 ); }

    double real() const { return real_; }

    double imaginary() const { return imaginary_; }

    void conjugate() { imaginary_ = -imaginary_; }

    // This could either test if both the real and the imaginary part are close to 0.0
    // or it could test if norm() or norm2() is close to 0.0. I have decided to do the
    // first.
    bool nearly_zero( const double tolerance = 0.000001 ) const;
    
    void reciprocal();

    void square();

    void power( const int n );

    double norm() const;

    double norm2() const;

    void show() const;

    Complex operator+() const { return Complex( *this ); }

    Complex operator-() const
    {
        Complex result( *this );
        result.negate();
        return result;
    }

    Complex & operator+=( const Complex & rhs );
    Complex & operator-=( const Complex & rhs );
    Complex & operator*=( const Complex & rhs );
    Complex & operator/=( const Complex & rhs );

    Complex & operator+=( const double & rhs );
    Complex & operator-=( const double & rhs );
    Complex & operator*=( const double & rhs );
    Complex & operator/=( const double & rhs );

    Complex & operator++();    // Prefix
    Complex   operator++(int); // Postfix
    Complex & operator--();    // Prefix
    Complex   operator--(int); // Postfix

//    bool operator==( const Complex & rhs ) const;
//    bool operator!=( const Complex & rhs ) const { return ! ( *this == rhs ); }
//    bool operator< ( const Complex & rhs ) const;
//    bool operator> ( const Complex & rhs ) const { return ( rhs < *this ); }
//    bool operator>=( const Complex & rhs ) const { return ! ( *this < rhs ); }
//    bool operator<=( const Complex & rhs ) const { return ! ( rhs < *this ); }

private:

    double real_;
    double imaginary_;

    inline void negate() { real_ = -real_; imaginary_ = -imaginary_; }

};

bool nearly_equal( const Complex & lhs, const Complex & rhs, const double tolerance = 0.000001 );

Complex operator+( const Complex & lhs, const Complex & rhs );
Complex operator-( const Complex & lhs, const Complex & rhs );
Complex operator*( const Complex & lhs, const Complex & rhs );
Complex operator/( const Complex & lhs, Complex rhs );

Complex operator+( const double lhs, const Complex & rhs );
Complex operator-( const double lhs, const Complex & rhs );
Complex operator*( const double lhs, const Complex & rhs );
Complex operator/( const double lhs, Complex rhs );

Complex operator+( const Complex & lhs, const double rhs );
Complex operator-( const Complex & lhs, const double rhs );
Complex operator*( const Complex & lhs, const double rhs );
Complex operator/( const Complex & lhs, const double rhs );

Complex square( const Complex value );

Complex exponential( const Complex z );

#endif // COMPLEX_H

