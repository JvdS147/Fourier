#ifndef UTILITIES_H
#define UTILITIES_H

/* *********************************************
Copyright (c) 2013-2024, Cornelis Jan (Jacco) van de Streek
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

#include <iostream> // @@ Only necessary for show(), this is bad
#include <string>
#include <vector>

// E.g. initialise_with_sequential_numbers( 5 ) returns { 0, 1, 2, 3, 4 }.
template <class T>
std::vector< T > initialise_with_sequential_values( const T & end_value )
{
    std::vector< T > result;
    result.reserve( end_value );
    T t(0);
    while ( t != end_value )
    {
        result.push_back( t );
        ++t;
    }
    return result;
}

template <class T>
void show( const std::vector< T > & values )
{
    for ( size_t i( 0 ); i != values.size(); ++i )
        std::cout << values[i] << std::endl;
}

template <class T>
std::vector< T > add( const std::vector< T > & lhs, const std::vector< T > & rhs )
{
    std::vector< T > result;
    result.reserve( lhs.size() + rhs.size() );
    for ( size_t i( 0 ); i != lhs.size(); ++i )
        result.push_back( lhs[i] );
    for ( size_t i( 0 ); i != rhs.size(); ++i )
        result.push_back( rhs[i] );
    return result;
}

template <class T>
bool contains( const std::vector< T > & values, const T & value )
{
    for ( size_t i( 0 ); i != values.size(); ++i )
    {
        if ( values[i] == value )
            return true;
    }
    return false;
}

template <class T>
bool nearly_contains( const std::vector< T > & values, const T & value )
{
    for ( size_t i( 0 ); i != values.size(); ++i )
    {
        if ( nearly_equal( values[i], value ) )
            return true;
    }
    return false;
}

// Outputs between 0 and max_size (inclusive) number of c characters, proportional to value.
// Basically generates a quick and dirty ASCII histogram (rotated by 90 degrees).
std::string ASCII_histogram( const double min, const double max, const double value, const size_t max_size, const char c );

// 1/0, true/false, t/f, yes/no, y/n, not case sensitive
bool string2bool( const std::string & input );

// For internal use only.
double string2double_2( const std::string & input, const bool float_allowed );

// Recognises scientific notation with "E" or "e" such as -.234e-45.
double string2double( std::string input );

int string2integer( const std::string & input );

size_t string2size_t( const std::string & input );

// The configurability of double2string() suggests that a class is called for.
// The class could be fed all numbers to be printed so that it can determine optimum values for length and precision
// E.g. explicit "+"
// E.g. allow_scientific_notation
class Double2String
{
public:

    // Default constructor.
    Double2String();

    bool explicit_plus() { return explicit_plus_; }
    void set_explicit_plus( const bool explicit_plus ) { explicit_plus_ = explicit_plus; }

    std::string convert( const double ) const;

private:
    bool explicit_plus_;
    size_t length_;
    size_t precision_;
    bool pad_plus_;

};

std::string double2string( const double input );

std::string double2string_2( const double input, const size_t ndecimals );

// Pads the string to e.g. "   1.000"
// If the length of the input value is longer than the padded length, a string with
// the length of the input value is returned (so the value is never corrupted).
// Precision takes precedence over maximum length, again the value (including its precision) is never corrupted
std::string double2string( const double input, const size_t precision, const size_t padded_length = 0, const char padding_character = ' ' );

// Aligns positive and negative numbers by prepending positive numbers with the padding_character.
// Output would look like this:
// "+0.1234"        or        " 0.1234"
// "-0.1234"        or        "-0.1234"
// "+0.1234"        or        " 0.1234"
// "-0.1234"        or        "-0.1234"
std::string double2string_pad_plus( const double input, const size_t precision, const char padding_character = ' ' );

// Pads the string to e.g. "0001"
// If the length of the input value is longer than the padded length, a string with
// the length of the input value is returned (so the value is never corrupted).
std::string size_t2string( const size_t input, const size_t padded_length = 0, const char padding_character = '0' );

// Pads the string to e.g. "0001"
// If the length of the input value is longer than the padded length, a string with
// the length of the input value is returned (so the value is never corrupted).
std::string int2string( const int input, const size_t padded_length = 0, const char padding_character = '0' );

std::string bool2string( const bool input );

#endif // UTILITIES_H

