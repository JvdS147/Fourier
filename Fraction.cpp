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

#include "Fraction.h"
#include "MathsFunctions.h"

#include <stdexcept>
#include <sstream>
#include <cmath>
#include <iostream> // For debugging

// ********************************************************************************

Fraction::Fraction( const int numerator, const int denominator )
 : integer_part_(0),
   numerator_(numerator),
   denominator_(denominator)
{
    clean_up();
}

// ********************************************************************************

Fraction::Fraction( const int integer_part, const int numerator, const int denominator )
 : integer_part_(integer_part),
   numerator_(numerator),
   denominator_(denominator)
{
    clean_up();
}

// ********************************************************************************

void Fraction::absolute()
{
    if ( ( integer_part_ < 0 ) ||
         ( numerator_    < 0 ) )
        negate();
}

// ********************************************************************************

void Fraction::reciprocal()
{
    int old_numerator = numerator_;
    numerator_   = denominator_;
    denominator_ = denominator_*integer_part_ + old_numerator;
    integer_part_ = 0;
    clean_up();
}

// ********************************************************************************

void Fraction::square()
{
    // Could potentially lead to overflows for temporary results.
    numerator_    = 2*integer_part_*numerator_*denominator_ + numerator_*numerator_;
    denominator_  *= denominator_;
    integer_part_ *= integer_part_;
    clean_up();
}

// ********************************************************************************

void Fraction::power( const int n )
{
    // Note: 0^0 = 1
    if ( n == 0 )
    {
        integer_part_ = 1;
        numerator_    = 0;
        denominator_  = 1;
        return;
    }
    if ( n < 0 )
        reciprocal();
    size_t abs_n = abs( n );
    if ( ( abs_n == 1 ) || is_zero() )
        return;
    if ( abs_n == 2 )
    {
        square();
        return;
    }
    // Many clever things can be done here, e.g. using the binary power method or recursion.
    // Note that we gain something from using the binary method because our square() function is optimised
    // for Fraction objects and faster than simply x * x (with x a Fraction).
    // The current implementation is extremely simple.
    Fraction fraction( *this );
    for ( size_t i(1); i != abs_n; ++i )
        *this *= fraction;
}

// ********************************************************************************

// We could define a private static function safe_add() that tries to add two fractions in such a way that
// GCDs for temporary results are factored out to prevent overflows for temporary results:
// safe_add( lhs_numerator, lhs_denominator, rhs_numerator, rhs_denominator )

Fraction & Fraction::operator+=( const Fraction & rhs )
{
    // Could potentially lead to overflows for temporary results.
    integer_part_ += rhs.integer_part_;
    numerator_     = rhs.denominator_*numerator_ + denominator_*rhs.numerator_;
    denominator_  *= rhs.denominator_;
    clean_up();
    return *this;
}

// ********************************************************************************

Fraction & Fraction::operator-=( const Fraction & rhs )
{
    // Could potentially lead to overflows for temporary results.
    integer_part_ -= rhs.integer_part_;
    numerator_     = rhs.denominator_*numerator_ - denominator_*rhs.numerator_;
    denominator_  *= rhs.denominator_;
    clean_up();
    return *this;
}

// ********************************************************************************

Fraction & Fraction::operator*=( const Fraction & rhs )
{
    // Could potentially lead to overflows for temporary results.
    numerator_     = integer_part_*rhs.numerator_*denominator_ + rhs.integer_part_*numerator_*rhs.denominator_ + numerator_*rhs.numerator_;
    denominator_  *= rhs.denominator_;
    integer_part_ *= rhs.integer_part_;
    clean_up();
    return *this;
}

// ********************************************************************************

Fraction & Fraction::operator/=( const Fraction & rhs )
{
    // Could potentially lead to overflows for temporary results.
    numerator_    = integer_part_*denominator_*rhs.denominator_     + numerator_*rhs.denominator_;
    denominator_  = rhs.integer_part_*denominator_*rhs.denominator_ + rhs.numerator_*denominator_;
    integer_part_ = 0;
    clean_up();
    return *this;
}

// ********************************************************************************

Fraction & Fraction::operator++()    // Prefix
{
    ++integer_part_;
    clean_up();
    return *this;
}

// ********************************************************************************

Fraction  Fraction::operator++(int) // Postfix
{
    Fraction old( *this );
    ++(*this);
    return old;
}

// ********************************************************************************

Fraction & Fraction::operator--()    // Prefix
{
    --integer_part_;
    clean_up();
    return *this;
}

// ********************************************************************************

Fraction  Fraction::operator--(int) // Postfix
{
    Fraction old( *this );
    --(*this);
    return old;
}

// ********************************************************************************

bool Fraction::operator==( const Fraction & rhs ) const
{
    return ( ( integer_part_ == rhs.integer_part_ ) &&
             ( numerator_    == rhs.numerator_    ) &&
             ( denominator_  == rhs.denominator_  ) );
}

// ********************************************************************************

bool Fraction::operator<( const Fraction & rhs ) const
{
    if ( integer_part_ < rhs.integer_part_ )
        return true;
    if ( rhs.integer_part_ < integer_part_ )
        return false;
    // When we are here, the integer parts are equal. We now want to test:
    //
    //         a    ?    c
    //        ---   <   ---
    //         b         d
    //
    // Because denominators b and d are always positive, we are allowed to multiply by them without swapping lhs / rhs:
    //
    //                ?
    //        d * a   <   b * c
    //
    // (This is a relatively expensive way of doing things, because it requires two multiplications, and
    // in principle the multiplications can even lead to overflows. The alternative would be a series of if-statements.)
    return ( (rhs.denominator_ * numerator_) < (denominator_ * rhs.numerator_) );
}

// ********************************************************************************

std::string Fraction::to_string() const
{
    std::stringstream o;
    if ( ! is_pure_fraction() )
        o << integer_part_;
    if ( ( ! is_pure_fraction() ) &&
         ( ! is_integer()       ) )
        o << " + ";
    if ( ! is_integer() )
        o << numerator_ << "/" << denominator_;
    return o.str();
}

// ********************************************************************************

void Fraction::clean_up()
{
    if ( denominator_ == 0 )
        throw std::runtime_error( "Fraction::clean_up(): divide by zero." );
    // The denominator is always positive.
    if ( denominator_ < 0 )
    {
        denominator_ = -denominator_;
        numerator_ = -numerator_;
    }
    // Factor out the integer part from the fraction part
    int integer_part = numerator_ / denominator_;
    integer_part_ += integer_part;
    numerator_ -= integer_part * denominator_;
    // Reduce the fraction part
    // 2/4 ==> 1/2
    // 0/7 ==> 0/1
    int common_factor = ( numerator_ < 0 ) ? greatest_common_divisor( -numerator_, denominator_ ) : greatest_common_divisor( numerator_, denominator_ );
    numerator_   /= common_factor;
    denominator_ /= common_factor;
    // When we are here, | numerator_ | < denominator_ and denominator_ always positive.
    // But integerPart_ and numerator_ can still have opposite signs ( "8 + -2/3" ).
    // Because the integer part is always numerically greater than the fraction part when we are here,
    // we must always adjust the fraction part to have the same sign as the integer part.
    if ( ( integer_part_ < 0 ) &&
         ( numerator_    > 0 ) )
    {
        numerator_ -= denominator_;
        ++integer_part_;
        clean_up();
    }
    if ( ( integer_part_ > 0 ) &&
         ( numerator_    < 0 ) )
    {
        numerator_ += denominator_;
        --integer_part_;
        clean_up();
    }
}

// ********************************************************************************

Fraction absolute( const Fraction & fraction )
{
    Fraction result( fraction );
    result.absolute();
    return result;
}

// ********************************************************************************

std::ostream & operator<<( std::ostream & os, const Fraction fraction )
{
    os << fraction.to_string();
    return os;
}

// ********************************************************************************

double   operator*( const Fraction & lhs, const double    rhs )
{
    return lhs.to_double() * rhs;
}

// ********************************************************************************

double   operator*( const double    lhs, const Fraction & rhs )
{
    return lhs * rhs.to_double();
}

// ********************************************************************************

Fraction operator*( const Fraction & lhs, const int       rhs )
{
    return Fraction( rhs * lhs.integer_part(),
                     rhs * lhs.numerator(),
                     lhs.denominator() );
}

// ********************************************************************************

Fraction operator*( const int       lhs, const Fraction & rhs )
{
    return Fraction( lhs * rhs.integer_part(),
                     lhs * rhs.numerator(),
                     rhs.denominator() );
}

// ********************************************************************************

Fraction operator*( const Fraction & lhs, const size_t    rhs )
{
    return Fraction( static_cast<int>(rhs) * lhs.integer_part(),
                     static_cast<int>(rhs) * lhs.numerator(),
                     lhs.denominator() );
}

// ********************************************************************************

Fraction operator*( const size_t    lhs, const Fraction & rhs )
{
    return Fraction( static_cast<int>(lhs) * rhs.integer_part(),
                     static_cast<int>(lhs) * rhs.numerator(),
                     rhs.denominator() );
}

// ********************************************************************************

double   operator+( const Fraction & lhs, const double    rhs )
{
    return lhs.to_double() + rhs;
}

// ********************************************************************************

double   operator+( const double    lhs, const Fraction & rhs )
{
    return lhs + rhs.to_double();
}

// ********************************************************************************

Fraction operator+( const Fraction & lhs, const int       rhs )
{
    return Fraction( lhs.integer_part() + rhs, lhs.numerator(), lhs.denominator() );
}

// ********************************************************************************

Fraction operator+( const int       lhs, const Fraction & rhs )
{
    return rhs + lhs;
}

// ********************************************************************************

Fraction operator+( const Fraction & lhs, const size_t    rhs )
{
    return Fraction( lhs.integer_part() + rhs, lhs.numerator(), lhs.denominator() );
}

// ********************************************************************************

Fraction operator+( const size_t    lhs, const Fraction & rhs )
{
    return rhs + lhs;
}

// ********************************************************************************

double   operator/( const Fraction & lhs, const double    rhs )
{
    return lhs.to_double() / rhs;
}

// ********************************************************************************

double   operator/( const double    lhs, const Fraction & rhs )
{
    return lhs / rhs.to_double();
}

// ********************************************************************************

Fraction operator/( const Fraction & lhs, const int       rhs )
{
    return Fraction( lhs.integer_part()*lhs.denominator() + lhs.numerator(),
                     lhs.denominator() * rhs );
}

// ********************************************************************************

Fraction operator/( const int       lhs, const Fraction & rhs )
{
    return Fraction( lhs * rhs.denominator(),
                     rhs.integer_part()*rhs.denominator() + rhs.numerator() );
}

// ********************************************************************************

Fraction operator/( const Fraction & lhs, const size_t    rhs )
{
    return Fraction( lhs.integer_part()*lhs.denominator() + lhs.numerator(),
                     lhs.denominator() * static_cast<int>(rhs) );
}

// ********************************************************************************

Fraction operator/( const size_t    lhs, const Fraction & rhs )
{
    return Fraction( static_cast<int>(lhs) * rhs.denominator(),
                     rhs.integer_part()*rhs.denominator() + rhs.numerator() );
}

// ********************************************************************************

double   operator-( const Fraction & lhs, const double    rhs )
{
    return lhs.to_double() - rhs;
}

// ********************************************************************************

double   operator-( const double    lhs, const Fraction & rhs )
{
    return lhs - rhs.to_double();
}

// ********************************************************************************

Fraction operator-( const Fraction & lhs, const int       rhs )
{
    return Fraction( lhs.integer_part() - rhs, lhs.numerator(), lhs.denominator() );
}

// ********************************************************************************

Fraction operator-( const int       lhs, const Fraction & rhs )
{
    return lhs + ( -rhs );
}

// ********************************************************************************

Fraction operator-( const Fraction & lhs, const size_t    rhs )
{
    return Fraction( lhs.integer_part() - static_cast<int>(rhs), lhs.numerator(), lhs.denominator() );
}

// ********************************************************************************

Fraction operator-( const size_t    lhs, const Fraction & rhs )
{
    return Fraction(lhs) + ( -rhs );
}

// ********************************************************************************

// Returns a Fraction from a double. smallest_unit is used for rounding.
// This is useful e.g. for translation components in space-group symmetry operators,
// which are multiples of 1/24. Fraction can be anything, e.g. "-1 + -3/4",
// but in practice only values like "1/24" or "1/100" are likely to make sense.
// (Note: translation components in space-group symmetry operators are really multiples
// of 1/24, not of 1/12, because there are cubic space groups like F 41 3 2 that have
// symmetry operators with translation components like 3/8, 5/8, 7/8.)
//
// Examples:
// Fraction fraction = double2fraction( 0.333, Fraction( 1, 24 ) ); // 1/3
// Fraction fraction = double2fraction( 0.2, Fraction( 1, 24 ) ); // 5/24 !! = 0.208333
// Fraction fraction = double2fraction( 0.25, Fraction( 1, 24 ) ); // 1/4
// Fraction fraction = double2fraction( 0.167, Fraction( 1, 24 ) ); // 1/6
// Fraction fraction = double2fraction( -0.333, Fraction( 1, 24 ) ); // -1/3
//
// We could make this a constructor, but a non-friend non-member function is perhaps better.
// The rounding error that is introduced can be calculated as follows:
// double original_value( 0.333 );
// Fraction fraction_from_double( double2fraction( original_value, Fraction( 1, 10 ) ) );
// // fraction_from_double is now 0.3, i.e. 3/10
// double error = fraction_from_double.to_double() - original_value;
// std::abs(error) is now 0.033 (within the numerical accuracy of a floating point number).
// Being able to calculate the error suggests that a natural way to generate general
// human-readable Fractions is:
//
//double initial_value( 0.2 );
//int best_i = 1;
//double smallest_error = std::abs( initial_value );
//for ( int i( 1 ); i <= 10; ++i )
//{
//    // Following is wrong if initial_value == 0.0?
//    Fraction trial_i = double2fraction( initial_value, Fraction( 1, i ) );
//    double error_i = std::abs( trial_i.to_double() - initial_value );
//    if ( error_i < smallest_error )
//    {
//        best_i = i;
//        smallest_error = error_i;
//    }
//}
// Fraction best_fraction( 1, best_i );
// Note that even though 1/24 is smaller than 1/10, and one would therefore expect that a granularity of 1/24
// can therefore always find a better approximation to a fraction than 1/10, this is not so:
// 0.2 (= 1/5) can be represented exactly using a granularity of 1/10 or even 1/5, whereas it cannot be represented
// exactly with 1/24.
// Actually, there may be an even better way: in e.g. "0.333", it is clear that the rounding error is <= 0.0005.
// So we are looking for the fraction n/m with the smallest values for n and m with a residual error smaller than 0.0005.
// So even though the residual error is 0.0 if we use 333/1000, the residual error is 0.000333..., which is < 0.0005, already for 1/3.
Fraction double2fraction( const double target, const Fraction & smallest_unit )
{
    double multiple( target / smallest_unit.to_double() );
    int multiple_as_int = ( multiple < 0 ) ? static_cast<int>( multiple - 0.5 ) : static_cast<int>( multiple + 0.5 );
    return ( multiple_as_int * smallest_unit );
}

// ********************************************************************************

Fraction Farey( double target, const int maximum_denominator )
{
    if ( maximum_denominator < 1 )
        throw std::runtime_error( "Farey(): maximum_denominator must be greater than 0." );

    bool is_negative = false;
    if ( target < 0.0 )
    {
        target = -target;
        is_negative = true;
    }
    double integer_part;
    target = modf( target, &integer_part );
    
    Fraction lower_limit( 0, 1 );
    Fraction upper_limit( 1, 1 ); // Note that this is stored as Fraction( 1, 0, 1 )
    Fraction mediant( 1, 2 );
    while ( ( lower_limit.denominator() + upper_limit.denominator() ) <= maximum_denominator )
    {
        mediant = Fraction( lower_limit.numerator() + upper_limit.numerator() + upper_limit.integer_part(), lower_limit.denominator() + upper_limit.denominator() );
        if ( mediant.to_double() < target )
            lower_limit = mediant;
        else
            upper_limit = mediant;
    }
    Fraction result = ( std::abs( lower_limit.to_double() - target ) < std::abs( upper_limit.to_double() - target ) ) ? lower_limit : upper_limit;

    result += Fraction( round_to_int( integer_part ) );
    if ( is_negative )
        result = -result;
        
    return result;
}

// ********************************************************************************

// Implementation is elegant but inefficient: every call to Farey() extracts the integer part
// and tests if the target is negative, which later get added back in again. Especially the first
// call, Farey( target, 1 ), which can only return 0 or 1, is ludicrously inefficient.
Fraction Farey_terminate_on_error( const double target, const double epsilon )
{
    size_t i( 1 );
    Fraction result = Farey( target, i );
    while ( ! ( std::abs( result.to_double() - target ) < epsilon ) )
    {
        ++i;
        result = Farey( target, i );
    }
    return result;
}

// ********************************************************************************

