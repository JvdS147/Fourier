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

#include "Utilities.h"
#include "BasicMathsFunctions.h" // For round_to_int(), but this has got to lead to circular references sooner or later
#include "StringFunctions.h"

#include <cmath>
#include <sstream>
#include <stdexcept>

// ********************************************************************************

std::string ASCII_histogram( const double min, const double max, const double value, const size_t max_size, const char c )
{
    size_t n;
    if ( value < min )
        n = 0;
    else if ( value > max )
        n = max_size;
    else // This function still suffers from rounding errors. What happens if the data point happens to be start_ or finish_?
        n = round_to_size_t( (max_size+1) * ( ( value - min ) / ( max - min ) ) );
    return std::string( n, c );
}

// ********************************************************************************

// 1/0, true/false, t/f, yes/no, y/n, not case sensitive
bool string2bool( const std::string & input )
{
    if ( input == "1" )
        return true;
    if ( input == "0" )
        return false;
    std::string input_2 = to_upper( input );
    if ( input_2 == "Y" )
        return true;
    if ( input_2 == "N" )
        return false;
    if ( input_2 == "YES" )
        return true;
    if ( input_2 == "NO" )
        return false;
    if ( input_2 == "T" )
        return true;
    if ( input_2 == "F" )
        return false;
    if ( input_2 == "TRUE" )
        return true;
    if ( input_2 == "FALSE" )
        return false;
    throw std::runtime_error( "string2bool(): input string >" + input + "< is not a valid boolean value." );
}

// ********************************************************************************

double string2double_2( const std::string & input, const bool float_allowed )
{
    if ( input.empty() )
        throw std::runtime_error( "string2double_2(): input string is empty." );
    double result( 0.0 );
    bool is_negative( false );
    size_t iPos( 0 );
    if ( input[iPos] == '+' )
        ++iPos;
    else if ( input[iPos] == '-' )
    {
        is_negative = true;
        ++iPos;
    }
    bool after_period( false );
    bool at_least_one_digit_found( false );
    double power( 1.0 );
    double exponent( 0.0 );
    bool exponent_found( false );
    while ( iPos < input.length() )
    {
        if ( input[iPos] == '.' )
        {
            if ( ! float_allowed )
                throw std::runtime_error( "string2double_2(): period found in integer value : >" + input + "<" );
            if ( after_period )
                throw std::runtime_error( "string2double_2(): second period found : >" + input + "<" );
            after_period = true;
            ++iPos;
            continue;
        }
        if ( ( input[iPos] == 'E' ) || ( input[iPos] == 'e' ) )
//        if ( ( input[iPos] == 'E' ) || ( input[iPos] == 'e' ) || ( input[iPos] == 'D' ) || ( input[iPos] == 'd' ) )
        {
            if ( ! float_allowed )
                throw std::runtime_error( "string2double_2(): exponent found in integer value : >" + input + "<" );
            if ( ! at_least_one_digit_found )
                throw std::runtime_error( "string2double(): no digits before exponent : >" + input + "<" );
            ++iPos;
            if ( iPos == input.length() )
                throw std::runtime_error( "string2double(): no digits after exponent : >" + input + "<" );
            exponent_found = true;
            exponent = string2double_2( input.substr( iPos ), false );
            break;
        }
        if ( is_digit( input[iPos] ) )
        {
            at_least_one_digit_found = true;
            double value = double( int( input[iPos] ) - int( '0' ) );
            if ( after_period )
                power /= 10.0;
            else
                result *= 10.0;
            result += value * power;
            ++iPos;
            continue;
        }
        throw std::runtime_error( "string2double_2(): invalid character found : >" + input + "<" );
    }
    if ( ! at_least_one_digit_found )
        throw std::runtime_error( "string2double_2(): no digits found : >" + input + "<" );
    if ( exponent_found )
        result *= std::pow( 10.0, exponent );
    if ( is_negative )
        result = -result;
    return result;
}

// ********************************************************************************

double string2double( std::string input )
{
    size_t iPos1 = input.find_first_of( "(" );
    size_t iPos2 = input.find_first_of( ")" );
    if ( iPos1 == std::string::npos )
    {
        if ( iPos2 != std::string::npos )
            throw std::runtime_error( "string2double(): parentheses not closed properly :  >" + input + "<" );
    }
    else
    {
        if ( iPos2 != input.length()-1 )
            throw std::runtime_error( "string2double(): parentheses not closed properly :  >" + input + "<" );
        input.erase( iPos1 );
    }
    return string2double_2( input, true );
}

// ********************************************************************************

int string2integer( const std::string & input )
{
    return round_to_int( string2double_2( input, false ) );
}

// ********************************************************************************

std::string double2string( const double input )
{
    std::stringstream string_stream;
    string_stream << input;
    return string_stream.str();
}

// ********************************************************************************

std::string double2string_2( const double input, const size_t ndecimals )
{
    std::stringstream string_stream;
    string_stream << std::fixed;
    string_stream.precision( ndecimals );
    string_stream << input;
    return string_stream.str();
}

// ********************************************************************************

std::string double2string( const double input, const size_t precision, const size_t padded_length, const char padding_character )
{
    std::stringstream string_stream;
    string_stream << std::fixed;
    string_stream.precision( precision );
    string_stream << input;
    std::string result = string_stream.str();
    while ( result.size() < padded_length )
        result = padding_character + result;
    return result;
}

// ********************************************************************************

std::string double2string_pad_plus( const double input, const size_t precision, const char padding_character )
{
    std::stringstream string_stream;
    string_stream << std::fixed;
    string_stream.precision( precision );
    string_stream << input;
    std::string result = string_stream.str();
    if ( input >= 0.0 )
        result = padding_character + result;
    return result;
}

// ********************************************************************************

std::string size_t2string( const size_t input, const size_t padded_length, const char padding_character )
{
    std::stringstream string_stream;
    string_stream << input;
    std::string result = string_stream.str();
    while ( result.size() < padded_length )
        result = padding_character + result;
    return result;
}

// ********************************************************************************

std::string int2string( const int input, const size_t padded_length, const char padding_character )
{
    std::stringstream string_stream;
    string_stream << input;
    std::string result = string_stream.str();
    while ( result.size() < padded_length )
        result = padding_character + result;
    return result;
}

// ********************************************************************************

std::string bool2string( const bool input )
{
    if ( input )
        return "1";
    return "0";
}

// ********************************************************************************

