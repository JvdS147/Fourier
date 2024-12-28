#ifndef STRINGFUNCTIONS_H
#define STRINGFUNCTIONS_H

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

#include <string>
#include <vector>

bool contains( const std::string & input, const std::string & word );

// It turns out to be surprisingly difficult to uppercase or lowercase a letter in C++.

char to_upper( const char argument );

char to_lower( const char argument );

std::string to_upper( const std::string & argument );

std::string to_lower( const std::string & argument );

bool is_upper_case_letter( const char argument );

bool is_lower_case_letter( const char argument );

bool is_letter( const char argument );

bool is_digit( const char argument );

// Removes all occurences of char c from input string.
std::string remove( const std::string & input, const char c );

// Removes all occurences of char c from the start of the input string, i.e. until another character is found.
std::string remove_from_start( const std::string & input, const char c );

// interlace( "1234", '#' ) returns "1#2#3#4"
std::string interlace( const std::string & input, const char c );

// Removes all spaces and tabs from the beginning and the end.
std::string strip( const std::string & input );

// Removes all spaces and tabs from the end.
std::string remove_trailing_whitespace( const std::string & input );

// Removes all spaces and tabs from the beginning and the end.
void strip( std::vector< std::string > & input );

// Removes delimiters if both present, removes exactly one from start and one from beginning.
// If the delimiters overlap, returns the original string:
// remove_delimiters( "123", "12", "23" );
// returns "123", not "".
std::string remove_delimiters( const std::string & input, const std::string & start_delimiter, const std::string & end_delimiter );

// Assumes that c is a character signifying the start of a comment, e.g. "'" in TOPAS or "#".
std::string remove_from( const std::string & input, const char c );

// Counts number of occurrences of a specified character.
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
// The second string can be empty.
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

// E.g. make_multiple( "# ", 3 ); returns "# # # "
std::string make_multiple( const std::string & input, const size_t number );

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

// Obsolete now, I think, should be replaced by check_if_quotes_correct().
bool is_enclosed_in_quotes( const std::string & input );

// Turns "C:\Data\file_name.txt" into "C:\\Data\\file_name.txt", necessary when writing input files for e.g. R.
std::string escape_slashes( const std::string & input );

// Expects lines like "zero-point error : 0.01", with Splitter(":"), returns "0.01" (without whitespace).
std::string extract_variable_value( const std::string & line, const Splitter & splitter );

#endif // STRINGFUNCTIONS_H

