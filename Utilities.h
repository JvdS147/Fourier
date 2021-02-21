#ifndef UTILITIES_H
#define UTILITIES_H

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

#include <vector>
#include <string>
#include <cmath>

// We must have a separate StringFunctions.h

// E.g. initialise_with_sequential_numbers( 5 ) returns { 0, 1, 2, 3, 4 }
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

// Outputs between 0 and max_size (inclusive) number of c characters, proportional to value.
// Basically generates a quick and dirty ASCII histogram (rotated by 90 degrees).
std::string ASCII_histogram( const double min, const double max, const double value, const size_t max_size, const char c );

// It turns out to be surprisingly difficult to uppercase or lowercase a letter in C++

char to_upper( const char argument );

char to_lower( const char argument );

std::string to_upper( const std::string & argument );

std::string to_lower( const std::string & argument );

bool is_upper_case_letter( const char argument );

bool is_lower_case_letter( const char argument );

bool is_digit( const char argument );

// Removes all occurences of char c from input string
std::string remove( const std::string & input, const char c );

// Removes all occurences of char c from the start of the input string, i.e. until another character is found
std::string remove_from_start( const std::string & input, const char c );

// interlace( "1234", '#' ) returns "1#2#3#4"
std::string interlace( const std::string & input, const char c );

// Removes all spaces and tabs from the beginning and the end
std::string strip( const std::string & input );

// Removes all spaces and tabs from the beginning and the end
void strip( std::vector< std::string > & input );

// Removes delimiters if both present, removes exactly one from start and one from beginning
// If the delimiters overlap, returns the original string:
// remove_delimiters( "123", "12", "23" );
// returns "123", not "".
std::string remove_delimiters( const std::string & input, const std::string & start_delimiter, const std::string & end_delimiter );

// Counts number of occurrences of a specified character
size_t count_characters( const std::string & input, const char character );

// centre( "1", 5, '#') returns "##1##".
std::string centre( const std::string & input, const size_t padded_length = 0, const char padding_character = ' ' );

bool is_whole_word_character( const char c );

// Pads the string to e.g. " 1.123  "
// or                      "-1.123  "
// If the length of the input string is longer than the padded length, a string with
// the length of the input string is returned (so the string is never corrupted).
std::string pad_plus( const std::string & input, const size_t padded_length = 0, const char padding_character = ' ' );

// Pads the string to e.g. "C3  "
// If the length of the input string is longer than the padded length, a string with
// the length of the input string is returned (so the string is never corrupted).
std::string pad( const std::string & input, const size_t padded_length = 0, const char padding_character = ' ' );

// Replaces a piece of text by another piece of text. The two do not have to have the same length.
// The second string can be empty
std::string replace( const std::string & input, const std::string & old_str, const std::string & new_str );

// Splits a line into individual words, currently the separator is hard-coded to be a space or a tab
// "one word" and 'one word' are recognised as one word, the quotes are stripped.
// Empty words (either multiple spaces or e.g. "") are not retained.
// The following examples are not allowed:
//     one two th"ree
//     one two th'ree
//     one "two three
// @@ We need to allow escaped quotes, e.g. "He said \"yes\"."
// This kind of configurability suggests that this could be a class.
std::vector< std::string > split( const std::string & input );

// Very simple implementation: returns all strings that were separated by the delimiter.
// No special provisions for empty string (multiple delimiters in a row), whitespace or quotes.
// E.g. split( "|", '|' ) returns two empty strings.
std::vector< std::string > split( const std::string & input, const char delimiter );

// Same as split, but single quotes are not treated specially
// Necessary to read files where a single quote is a comment identifier (such as TOPAS .inp files)
std::vector< std::string > split_2( const std::string & input );

// We really have a problem here with empty words: allow or not? What should split( ",", "," ); return: two empty words?
// There should be an option to remove quotes from quoted fields
// There should be an option to configure what is used to delimit a quote
class Splitter
{
public:

    // Default constructor: delimiters are space and tab
    Splitter();

    explicit Splitter( const char* delimiters );

    // There are essentially two modes: merge delimiters and do not allow empty words or do not merge delimiters and do allow empty words
    bool merge_delimiters() { return merge_delimiters_; }
    void set_merge_delimiters( const bool merge_delimiters ) { merge_delimiters_ = merge_delimiters; }
    
    void split_by_length( const size_t split_length );

    std::vector< std::string > split( const std::string & input ) const;

private:
    std::vector< char > delimiters_;
    bool merge_delimiters_;
    bool split_by_length_;
    size_t split_length_;

    bool is_delimiter( const char c ) const;
};

// Allows easy extraction e.g. by typing:
// extract_delimited_text( "function(argument)", "(", ")" );
// extract_delimited_text( "He said \"quote\".", "\"", "\"" );
// If start and end delimiter are the same, the first occurrence from the start and the last occurrence from the end of the input are used.
// Returns empty string if delimiters not found or if the extracted string just happens to be empty.
std::string extract_delimited_text( const std::string & input, const std::string & start_delimiter, const std::string & end_delimiter );

// For internal use only.
double string2double_2( const std::string & input, const bool float_allowed );

// Recognises scientific notation with "E" or "e" such as -.234e-45
double string2double( std::string input );

int string2integer( const std::string & input );

// The configurability of double2string() suggests that a class is called for.
// The class could be fed all numbers to be printed so that it can determine optimum values for length and precision
// E.g. explicit "+"
// E.g. allow_scientific_notation
class Double2String
{
public:

    // Default constructor
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

// E.g. make_multiple( "# ", 3 ); returns "# # # "
std::string make_multiple( const std::string & input, const size_t number );

// Pads the string to e.g. "0001"
// If the length of the input value is longer than the padded length, a string with
// the length of the input value is returned (so the value is never corrupted).
std::string size_t2string( const size_t input, const size_t padded_length = 0, const char padding_character = '0' );

// Pads the string to e.g. "0001"
// If the length of the input value is longer than the padded length, a string with
// the length of the input value is returned (so the value is never corrupted).
std::string int2string( const int input, const size_t padded_length = 0, const char padding_character = '0' );

inline bool nearly_equal( const double lhs, const double rhs, const double tolerance = 0.0000001 )
{
    return ( std::abs( rhs - lhs ) < tolerance );
}

inline bool nearly_zero( const double lhs, const double tolerance = 0.0000001 )
{
    return ( std::abs( lhs ) < tolerance );
}

// Checks is the string is correctly enquoted, i.e.:
// There are no quotes at all.
// All quotes are escaped.
// There is exactly one non-escaped quote at the start and exactly one non-escaped quote at the end.
// Otherwise throws std::runtime_error
// returns true if there were quotes
void check_if_quotes_correct( const std::string & input );

  // @@ Could add "force_quotes()" because some applications need that. Perhaps better, add "add_quotes()" to Utilities.h

// Takes care of quotes if at least one string is in quotes.
// What to do if one is in quotes and the other is not in quotes but surrounded by whitespace? strip() or not?
// std::string add_strings( const std::string & lhs, const std::string & rhs );

// Assumes check_if_quotes_correct() has been called.
std::string remove_quotes( const std::string & input );

// Obsolete now, I think, should be replaced by check_if_quotes_correct()
bool is_enclosed_in_quotes( const std::string & input );

// Turns "C:\Data\file_name.txt" into "C:\\Data\\file_name.txt", necessary when writing input files for e.g. R
std::string escape_slashes( const std::string & input );

std::string append_backslash( const std::string & input );

#endif // UTILITIES_H

