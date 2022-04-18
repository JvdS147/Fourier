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

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string.h> // For strlen()


// ********************************************************************************

bool contains( const std::string & input, const std::string & word )
{
    return ( input.find( word ) != std::string::npos );
}

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

char to_upper( const char argument )
{
    return ( ( 'z' < argument ) || ( argument < 'a' ) ) ? argument : char( int(argument) + int('A') - int('a') );
}

// ********************************************************************************

char to_lower( const char argument )
{
    return ( ( 'Z' < argument ) || ( argument < 'A' ) ) ? argument : char( int(argument) + int('a') - int('A') );
}

// ********************************************************************************

std::string to_upper( const std::string & argument )
{
    std::string result;
    for ( size_t i( 0 ); i != argument.length(); ++i )
        result += to_upper( argument[i] );
    return result;
}

// ********************************************************************************

std::string to_lower( const std::string & argument )
{
    std::string result;
    for ( size_t i( 0 ); i != argument.length(); ++i )
        result += to_lower( argument[i] );
    return result;
}

// ********************************************************************************

bool is_upper_case_letter( const char argument )
{
    return ( ( 'A' <= argument ) && ( argument <= 'Z' ) );
}

// ********************************************************************************

bool is_lower_case_letter( const char argument )
{
    return ( ( 'a' <= argument ) && ( argument <= 'z' ) );
}

// ********************************************************************************

bool is_digit( const char argument )
{
    return ( ( '0' <= argument ) && ( argument <= '9' ) );
}

// ********************************************************************************

// Removes all occurences of char c from input string
std::string remove( const std::string & input, const char c )
{
    std::string result;
    size_t iPos( 0 );
    while ( iPos < input.length() )
    {
        if ( input[iPos] != c )
            result += input[iPos];
        ++iPos;
    }
    return result;
}

// ********************************************************************************

// Removes all occurences of char c from the start of the input string, i.e. until another character is found
std::string remove_from_start( const std::string & input, const char c )
{
    size_t iPos( 0 );
    while ( ( iPos < input.length() ) && ( input[iPos] == c ) )
    {
        ++iPos;
    }
    return input.substr( iPos );
}

// ********************************************************************************

std::string interlace( const std::string & input, const char c )
{
    std::string result;
    if ( ! input.empty() )
    {
        for ( size_t i( 0 ); i != input.length()-1; ++i )
        {
            result += input[i];
            result += c;
        }
        result += input[input.length()-1];
    }
    return result;
}

// ********************************************************************************

std::string strip( const std::string & input )
{
    // Absorb white space from start
    size_t iStart( 0 );
    while ( ( iStart < input.length() ) && ( ( input[iStart] == ' ' ) || ( input[iStart] == '\t' ) ) )
        ++iStart;
    if ( iStart == input.length() )
        return std::string( "" );
    // Absorb white space from end
    size_t iEnd( input.length() - 1 );
    while ( ( input[iEnd] == ' ' ) || ( input[iEnd] == '\t' ) )
        --iEnd;
    return input.substr( iStart, iEnd - iStart + 1 );
}

// ********************************************************************************

void strip( std::vector< std::string > & input )
{
    for ( size_t i( 0 ); i != input.size(); ++i )
        input[i] = strip( input[i] );    
}

// ********************************************************************************

std::string remove_delimiters( const std::string & input, const std::string & start_delimiter, const std::string & end_delimiter )
{
    if ( input.length() < ( start_delimiter.length() + end_delimiter.length() ) )
        return input;
    if ( ( input.substr( 0, start_delimiter.length() ) == start_delimiter ) && ( input.substr( input.length() - end_delimiter.length(), input.length() ) == end_delimiter ) )
        return input.substr( start_delimiter.length(), input.length() - ( start_delimiter.length() + end_delimiter.length() ) );
    return input;
}

// ********************************************************************************

size_t count_characters( const std::string & input, const char character )
{
    size_t result( 0 );
    size_t iPos( 0 );
    while ( iPos < input.length() )
    {
        if ( input[iPos] == character )
            ++result;
        ++iPos;
    }
    return result;
}

// ********************************************************************************

std::string centre( const std::string & input, const size_t padded_length, const char padding_character )
{
    std::string result( input );
    while ( result.size() < padded_length )
    {
        result = padding_character + result;
        if ( result.size() < padded_length )
            result = result + padding_character;
    }
    return result;
}

// ********************************************************************************

bool is_whole_word_character( const char c )
{
    return ( isalnum( c ) || ( c == '_') );
}

// ********************************************************************************

// Pads the string to e.g. "C3  "
// If the length of the input string is longer than the padded length, a string with
// the length of the input string is returned (so the string is never corrupted).
std::string pad( const std::string & input, const size_t padded_length, const char padding_character )
{
    std::string result( input );
    while ( result.size() < padded_length )
        result += padding_character;
    return result;
}

// ********************************************************************************

// Pads the string to e.g. " 1.123  "
// or                      "-1.123  "
// If the length of the input string is longer than the padded length, a string with
// the length of the input string is returned (so the string is never corrupted).
std::string pad_plus( const std::string & input, const size_t padded_length, const char padding_character )
{
    std::string result( input );
    if ( ( input.length() > 0 ) && ( input[0] != '-' ) && ( input[0] != '+' ) )
      result = padding_character + result;
    while ( result.size() < padded_length )
        result += padding_character;
    return result;
}

// ********************************************************************************

// Replaces a piece of text by another piece of text. The two do not have to have the same length.
// The second string can be empty
std::string replace( const std::string & input, const std::string & old_str, const std::string & new_str )
{
    std::string result( input );
    size_t iPos = result.find( old_str );
    while ( iPos != std::string::npos )
    {
        result.replace( iPos, old_str.length(), new_str );
        iPos = result.find( old_str );
    }
    return result;
}

// ********************************************************************************

std::vector< std::string > split( const std::string & input )
{
    std::vector< std::string > result;
    std::string current_word;
    size_t i( 0 );
    while ( i < input.length() )
    {
        // Absorb white space
        while ( ( i < input.length()) && ( (input[i] == ' ') || (input[i] == '\t') ) )
            ++i;
        // Parse "one word"
        if ( ( i < input.length()) && (input[i] == '"') )
        {
            ++i;
            while ( ( i < input.length() ) && ( input[i] != '"' ) )
            {
                current_word += input[i];
                ++i;
            }
            if ( i == input.length() )
                throw std::runtime_error( "split(): double quote is not terminated properly: |" + input + "|" );
            ++i; // Read past the quote
            // We must now hit the end of the line or whitespace
            if ( ( i < input.length() ) && ( ( input[i] != ' ' ) && ( input[i] != '\t' ) ) )
                throw std::runtime_error( "split(): quote inside string is not allowed: |" + input + "|" );
        }
        // Parse 'one word'
        else if ( ( i < input.length() ) && ( input[i] == '\'' ) )
        {
            ++i;
            while ( ( i < input.length() ) && ( input[i] != '\'' ) )
            {
                current_word += input[i];
                ++i;
            }
            if ( i == input.length() )
                throw std::runtime_error( "split(): single quote is not terminated properly: |" + input + "|" );
            ++i; // Read past the quote
            // We must now hit the end of the line or whitespace
            if ( ( i < input.length() ) && ( ( input[i] != ' ' ) && ( input[i] != '\t' ) ) )
                throw std::runtime_error( "split(): quote inside string is not allowed: |" + input + "|" );
        }
        // Collect non-quoted characters
        else
        {
            while ( ( i < input.length() ) && ( ( input[i] != ' ' ) && ( input[i] != '\t' ) ) )
            {
//                if ( ( input[i] == '"' ) || ( input[i] == '\'' ) )
//                    throw std::runtime_error( "split(): quote inside string is not allowed: |" + input + "|" );
                current_word += input[i];
                ++i;
            }
        }
        if ( current_word.length() != 0 )
        {
            result.push_back( current_word );
            current_word = "";
        }
    }
    return result;
}

// ********************************************************************************

std::vector< std::string > split( const std::string & input, const char delimiter )
{
    std::vector< std::string > result;
    std::string current_string;
    size_t iPos( 0 );
    while ( iPos != input.length() )
    {
        if ( input[iPos] == delimiter )
        {
            result.push_back( current_string );
            current_string = "";
        }
        else
            current_string += input[iPos];
        ++iPos;
    }
    result.push_back( current_string );
    return result;
}

// ********************************************************************************

std::vector< std::string > split_2( const std::string & input )
{
    std::vector< std::string > result;
    std::string current_word;
    size_t i( 0 );
    while ( i < input.length() )
    {
        // Absorb white space
        while ( ( i < input.length() ) && ( ( input[i] == ' ' ) || ( input[i] == '\t' ) ) )
            ++i;
        // Parse "one word"
        if ( ( i < input.length() ) && ( input[i] == '"' ) )
        {
            ++i;
            while ( ( i < input.length() ) && ( input[i] != '"' ) )
            {
                current_word += input[i];
                ++i;
            }
            if ( i == input.length() )
                throw std::runtime_error( "split(): double quote is not terminated properly: |" + input + "|" );
            ++i; // Read past the quote
            // We must now hit the end of the line or whitespace
            if ( ( i < input.length() ) && ( ( input[i] != ' ' ) && ( input[i] != '\t' ) ) )
                throw std::runtime_error( "split(): quote inside string is not allowed: |" + input + "|" );
        }
        // Collect non-quoted characters
        else
        {
            while ( ( i < input.length() ) && ( (input[i] != ' ' ) && ( input[i] != '\t' ) ) )
            {
//                if ( ( input[i] == '"' ) || ( input[i] == '\'' ) )
//                    throw std::runtime_error( "split(): quote inside string is not allowed: |" + input + "|" );
                current_word += input[i];
                ++i;
            }
        }
        if ( current_word.length() != 0 )
        {
            result.push_back( current_word );
            current_word = "";
        }
    }
    return result;
}

// ********************************************************************************

Splitter::Splitter(): merge_delimiters_(true),split_by_length_(false),split_length_(0)
{
    delimiters_.push_back( ' ' );
    delimiters_.push_back( '\t' );
}

// ********************************************************************************

void Splitter::split_by_length( const size_t split_length )
{
    if ( split_length == 0 )
        throw std::runtime_error( "Splitter::split_by_length(): length cannot be 0." );
    split_by_length_ = true;
    split_length_ = split_length;
}

// ********************************************************************************

Splitter::Splitter( const char* delimiters ): merge_delimiters_(true),split_by_length_(false),split_length_(0)
{
    for ( size_t i( 0 ); i != strlen( delimiters ); ++i )
        delimiters_.push_back( delimiters[i] );
}

// ********************************************************************************

std::vector< std::string > Splitter::split( const std::string & input ) const
{
    std::vector< std::string > result;
    if ( split_by_length_ )
    {
        if ( split_length_ == 0 )
            throw std::runtime_error( "Splitter::split(): programming error." );
        size_t iStart( 0 );
        while ( iStart < input.length() )
        {
            result.push_back( input.substr( iStart, split_length_) );
            iStart += split_length_;
        }
        return result;
    }
    std::string current_word;
    size_t i( 0 );
    if ( ! merge_delimiters_ )
    {

//N1      C2     ,   1.47566, 1.47769`, bond_width, bond_weight

        while ( i < input.length() )
        {
            // Parse "one word"
            if ( ( i < input.length() ) && (input[i] == '"') )
            {
                ++i;
                while ( ( i < input.length() ) && (input[i] != '"') )
                {
                    current_word += input[i];
                    ++i;
                }
                if ( i == input.length() )
                    throw std::runtime_error( "split(): double quote is not terminated properly: |" + input + "|" );
                ++i; // Read past the quote
                // We must now hit the end of the line or whitespace
                if ( ( i < input.length() ) && ( ! is_delimiter( input[i] ) ) )
                    throw std::runtime_error( "split(): quote inside string is not allowed: |" + input + "|" );
            }
            // Parse 'one word'
            else if ( ( i < input.length() ) && ( input[i] == '\'' ) )
            {
                ++i;
                while ( ( i < input.length() ) && ( input[i] != '\'' ) )
                {
                    current_word += input[i];
                    ++i;
                }
                if ( i == input.length() )
                    throw std::runtime_error( "split(): single quote is not terminated properly: |" + input + "|" );
                ++i; // Read past the quote
                // We must now hit the end of the line or whitespace
                if ( ( i < input.length() ) && ( ! is_delimiter( input[i] ) ) )
                    throw std::runtime_error( "split(): quote inside string is not allowed: |" + input + "|" );
            }
            // Collect non-quoted characters
            else
            {
                while ( ( i < input.length() ) && ( ! is_delimiter( input[i] ) ) )
                {
    //                if ( ( input[i] == '"' ) || ( input[i] == '\'' ) )
    //                    throw std::runtime_error( "split(): quote inside string is not allowed: |" + input + "|" );
                    current_word += input[i];
                    ++i;
                }
            }
            result.push_back( current_word );
            current_word = "";
            // When we are here, either i == input.length() or input[i] is a delimiter.
            ++i; // Read the delimiter
        }
    }
    else // Merge delimiters
    {
        while ( i < input.length() )
        {
            // Absorb white space
            while ( ( i < input.length() ) && is_delimiter( input[i] ) )
                ++i;
            // Parse "one word"
            if ( ( i < input.length() ) && ( input[i] == '"' ) )
            {
                ++i;
                while ( ( i < input.length() ) && ( input[i] != '"' ) )
                {
                    current_word += input[i];
                    ++i;
                }
                if ( i == input.length() )
                    throw std::runtime_error( "split(): double quote is not terminated properly: |" + input + "|" );
                ++i; // Read past the quote
                // We must now hit the end of the line or whitespace
                if ( ( i < input.length() ) && ( ! is_delimiter( input[i] ) ) )
                    throw std::runtime_error( "split(): quote inside string is not allowed: |" + input + "|" );
            }
            // Parse 'one word'
            else if ( ( i < input.length() ) && ( input[i] == '\'' ) )
            {
                ++i;
                while ( ( i < input.length() ) && ( input[i] != '\'' ) )
                {
                    current_word += input[i];
                    ++i;
                }
                if ( i == input.length() )
                    throw std::runtime_error( "split(): single quote is not terminated properly: |" + input + "|" );
                ++i; // Read past the quote
                // We must now hit the end of the line or whitespace
                if ( ( i < input.length() ) && ( ! is_delimiter( input[i] ) ) )
                    throw std::runtime_error( "split(): quote inside string is not allowed: |" + input + "|" );
            }
            // Collect non-quoted characters
            else
            {
                while ( ( i < input.length() ) && ( ! is_delimiter( input[i] ) ) )
                {
    //                if ( ( input[i] == '"' ) || ( input[i] == '\'' ) )
    //                    throw std::runtime_error( "split(): quote inside string is not allowed: |" + input + "|" );
                    current_word += input[i];
                    ++i;
                }
            }
            if ( ! current_word.empty() )
            {
                result.push_back( current_word );
                current_word = "";
            }
        }
    }
    return result;
}

// ********************************************************************************

bool Splitter::is_delimiter( const char c ) const
{
     for ( size_t i( 0 ); i != delimiters_.size(); ++i )
     {
         if ( delimiters_[i] == c )
             return true;
     }
     return false;
}

// ********************************************************************************

// If start and end delimiter are the same, the first occurrence from the start and the last occurrence from the end of the input are used.
// Returns empty string if delimiters not found or if the extracted string just happens to be empty.
std::string extract_delimited_text( const std::string & input, const std::string & start_delimiter, const std::string & end_delimiter )
{
    if ( start_delimiter.empty() || end_delimiter.empty() )
        throw std::runtime_error( "extract_delimited_text(): one of the delimiters is empty" );
    size_t iStart = input.find( start_delimiter );
    if ( iStart == std::string::npos )
        return "";
    // Finding the last match from the end is not trivial...
    size_t iPos = input.length();
    while ( ( iPos > iStart ) && ( input.find( end_delimiter, iPos ) == std::string::npos ) )
        --iPos;
    if ( iPos == iStart )
        return "";
    return input.substr( iStart+start_delimiter.length(), iPos-(iStart+start_delimiter.length()) );
}

// ********************************************************************************

double string2double_2( const std::string & input, const bool float_allowed )
{
    if ( input.empty() )
        throw std::runtime_error( "string2double_2(): input string is empty" );
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

// E.g. make_multiple( "# ", 3 ); returns "# # # "
std::string make_multiple( const std::string & input, const size_t number )
{
    std::string result;
    for ( size_t i( 0 ); i != number; ++i )
        result += input;
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

void check_if_quotes_correct( const std::string & input )
{
}

// ********************************************************************************

std::string remove_quotes( const std::string & input )
{
    return remove_delimiters( input, "\"", "\"" );
}

// ********************************************************************************

bool is_enclosed_in_quotes( const std::string & input )
{
    if ( input.empty() )
        return false;
    if ( input.length() == 1 )
    {
        if ( input == "\"" )
            throw std::runtime_error( "is_enclosed_in_quotes(): input consists of a single double quote." );
        else
            return false;
    }
    return ( ( input.substr( 0, 1 ) == "\"" ) && ( input.substr(input.length()-1,1) == "\"" ) );
    // Doesn't check if the string contains quotes within the input. Of course, those quotes could be escaped with a \.
}

// ********************************************************************************

std::string escape_slashes( const std::string & input )
{
    std::string result;
    for ( size_t i( 0 ); i != input.length(); ++i )
    {
       if ( input[i] == '\\' )
           result += '\\';
       result += input[i];
    }
    return result;
}

// ********************************************************************************

std::string append_backslash( const std::string & input )
{
    if ( input.empty() || ( input.substr( input.length() - 1 ) == "\\" ) )
        return input;
    return input + "\\";
}

// ********************************************************************************
