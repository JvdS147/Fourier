#ifndef FRACTION_H
#define FRACTION_H

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

#include <iosfwd>
#include <string>

/*
 * Fraction class.
 *
 * Fractions are automatically reduced:
 *
 * Fraction fraction( 3, 6 );
 * fraction.numerator() now returns 1,
 * fraction.denominator() now returns 2.
 *
 * Fractions have an integer part:
 *
 * Fraction fraction( 6, 6 );
 * fraction.integerPart() now returns 1,
 * fraction.numerator() now returns 0,
 * fraction.denominator() now returns 1.
 *
 * The denominator is always made positive.
 *
 * The class expects that its numerator, its denumerator and its integer part can be stored in an int, the
 * class does not expect to store large integers.
 *
 */
class Fraction
{
public:

    explicit Fraction( const int integer_part = 0 ) : integer_part_(integer_part), numerator_(0), denominator_(1) {}

    // Must be of the type "1+4/5", "-1-4/5". "-1-4/-5" or "-1-4/+5" works but is discouraged.
    // Adding this constructor would make it necessary to add Utilities.h to Fraction.cpp
//    explicit Fraction( const std::string fraction_str );

    Fraction( const int numerator, const int denominator );

    // The Fraction is constructed to have the value "integerPart + ( numerator / denominator )",
    // but fractions are reduced, integer parts are factored out and the integer part and the fraction part are
    // given the same sign. E.g. Fraction( 8, 6, -4 ) is constructed as
    // 8 + ( 6 / -4 ), which is
    // 8 + -1 + ( 2 / -4 ), which is
    // 7 + ( -1 / 2 ), which is
    // 6 + ( 1 / 2 )
    // Fraction(-7,4,4) is constructed as "-6 + (0/1)"
    // Fraction(-7,2,-4) is constructed as "-7 + (-1/2)"
    Fraction( const int integer_part, const int numerator, const int denominator );

    static Fraction zero() { return Fraction( 0 ); }
    static Fraction one() { return Fraction( 1 ); }

    inline int integer_part() const { return integer_part_; }
    inline Fraction fractional_part() const { Fraction result; result.numerator_ = this->numerator_; result.denominator_ = this->denominator_; return result; }
    inline int numerator()    const { return numerator_; }
    inline int denominator()  const { return denominator_; }

    // Could make this such that it allows implicit conversions
    double to_double() const { return static_cast<double>( integer_part_ ) + (static_cast<double>( numerator_ )/static_cast<double>( denominator_ )) ; }

    // Returns true if the fraction part is 0
    bool is_integer() const { return ( numerator_ == 0 ); }

    // Returns true if the integer part is 0 and the fraction part is non-zero
    bool is_pure_fraction() const { return ( ( integer_part_ == 0 ) && ( numerator_ != 0 ) ); }

    // Returns true if the Fraction equals 0
    bool is_zero() const { return ( ( integer_part_ == 0 ) && ( numerator_ == 0 ) ); }

    // The absolute value is taken.
    void absolute();

    // The reciprocal value is taken. Throws if the fraction is 0.
    void reciprocal();

    // Raises the Fraction to the power 2
    void square();

    // Raises the Fraction to n. Fraction( 0 ).power( 0 ) is 1. Throws if the fraction is 0 and n is negative.
    void power( const int n );

    Fraction operator+( const Fraction & rhs ) const { return Fraction(*this) += rhs; }
    Fraction operator-( const Fraction & rhs ) const { return Fraction(*this) -= rhs; }
    Fraction operator*( const Fraction & rhs ) const { return Fraction(*this) *= rhs; }
    Fraction operator/( const Fraction & rhs ) const { return Fraction(*this) /= rhs; }

    Fraction operator+() const { return Fraction( *this ); }

    Fraction operator-() const
    {
        Fraction result( *this );
        result.negate();
        return result;
    }

    Fraction & operator+=( const Fraction & rhs );
    Fraction & operator-=( const Fraction & rhs );
    Fraction & operator*=( const Fraction & rhs );
    Fraction & operator/=( const Fraction & rhs );

    Fraction & operator++();    // Prefix
    Fraction   operator++(int); // Postfix
    Fraction & operator--();    // Prefix
    Fraction   operator--(int); // Postfix

    bool operator==( const Fraction & rhs ) const;
    bool operator!=( const Fraction & rhs ) const { return ! ( *this == rhs ); }
    bool operator< ( const Fraction & rhs ) const;
    bool operator> ( const Fraction & rhs ) const { return ( rhs < *this ); }
    bool operator>=( const Fraction & rhs ) const { return ! ( *this < rhs ); }
    bool operator<=( const Fraction & rhs ) const { return ! ( rhs < *this ); }

    // Returns the Fraction in string form, e.g. "0", "3", "1/2", "2 + 2/3", "-9 + -1/4"
    std::string to_string() const;

private:

    int integer_part_; // Always same sign as numerator_.
    int numerator_;   // Always same sign as integerPart_. Always smaller than denominator_.
    int denominator_; // Always positive. Always greater than numerator_. Always 1 when numerator_ is 0.

    void clean_up();

    inline void negate() { integer_part_ = -integer_part_; numerator_ = -numerator_; }

};

Fraction absolute( const Fraction & fraction );

std::ostream & operator<<( std::ostream & os, const Fraction fraction );

// We could remove all the overloads taking one double argument and either make it the
// responsibility of the programmer to convert the Fraction to a double first, or we could
// add a conversion operator to allow implicit conversion from a Fraction to a double.
// These should probably be inline, most of them could be replaced by member functions
// Half of them could be expressed in terms of the other half.
// Several of them would not be necessary if there were implicit conversions to double or from int.
double   operator*( const Fraction & lhs, const double     rhs );
double   operator*( const double     lhs, const Fraction & rhs );
Fraction operator*( const Fraction & lhs, const int        rhs );
Fraction operator*( const int        lhs, const Fraction & rhs );
Fraction operator*( const Fraction & lhs, const size_t     rhs );
Fraction operator*( const size_t     lhs, const Fraction & rhs );

double   operator+( const Fraction & lhs, const double     rhs );
double   operator+( const double     lhs, const Fraction & rhs );
Fraction operator+( const Fraction & lhs, const int        rhs );
Fraction operator+( const int        lhs, const Fraction & rhs );
Fraction operator+( const Fraction & lhs, const size_t     rhs );
Fraction operator+( const size_t     lhs, const Fraction & rhs );

double   operator/( const Fraction & lhs, const double     rhs );
double   operator/( const double     lhs, const Fraction & rhs );
Fraction operator/( const Fraction & lhs, const int        rhs );
Fraction operator/( const int        lhs, const Fraction & rhs );
Fraction operator/( const Fraction & lhs, const size_t     rhs );
Fraction operator/( const size_t     lhs, const Fraction & rhs );

double   operator-( const Fraction & lhs, const double     rhs );
double   operator-( const double     lhs, const Fraction & rhs );
Fraction operator-( const Fraction & lhs, const int        rhs );
Fraction operator-( const int        lhs, const Fraction & rhs );
Fraction operator-( const Fraction & lhs, const size_t     rhs );
Fraction operator-( const size_t     lhs, const Fraction & rhs );

// @@ We need bool operator<( const Fraction & lhs, const double rhs ); etc.

Fraction double2fraction( const double target, const Fraction & smallest_unit );

// Finds the best approximation to a floating point value for which the
// denominator is smaller than or equal to maximum_denominator
Fraction Farey( double target, const int maximum_denominator );

Fraction Farey_terminate_on_error( const double target, const double epsilon );

// To calculate the average of four values:
// Fraction average = average( value_1, value_2 );
// average = average( value_3, average, Fraction( 2 ) );
// average = average( value_4, average, Fraction( 3 ) );
// Even this works:
//    Fraction prev_estimate; // No need to initialise...
//    Fraction next_estimate; // No need to initialise...
//    size_t iStep( 0 );
//    while ( iStep < 1000000 )
//    {
//        prev_estimate = next_estimate;
//        Fraction current_value = some_function();
//        next_estimate = average( current_value, prev_estimate, Fraction( iStep ) );
//        ++iStep;
//    }
Fraction average( const Fraction & lhs, const Fraction & rhs, const Fraction & weight = Fraction( 1 ) );

bool triquality( const Fraction & x1, const Fraction & x2, const Fraction & x3 );

#endif // FRACTION_H

