#ifndef TEXTFILEREADER_H
#define TEXTFILEREADER_H

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

class FileName;

#include <fstream>
#include <string>
#include <vector>

// Reads one line at a time
// New lines:
// Windows       : \r\n
// Linux         : \n
// MacOS < 10    : \r
// MacOS >= 10   : \n
// C++ in principle uses (e.g. absorbs upon reading) \n
// This class deletes all \r characters from all input.

// Example code:
//        TextFileReader text_file_reader( FileName( "file_name.txt" ) );
//        text_file_reader.set_skip_empty_lines( false );
//        std::string line;
//        std::vector< std::string > words;
//        while ( text_file_reader.get_next_line( line ) )
//        while ( text_file_reader.get_next_line( words ) )
//        {
//            if ( words.size() != 4 )
//                std::cout << "Warning: words.size() != 4" << std::endl;
//            ...
//        }

class TextFileReader
{
public:

//    TextFileReader();

    explicit TextFileReader( const FileName & file_name );

    // For testing
 //   explicit TextFileReader( const std::vector< std::string > & lines );

    ~TextFileReader() { input_file_.close(); }

    bool get_next_line( std::vector< std::string > & words );

    bool get_next_line( std::string & line );

    // Returns the current line (e.g. for error reporting)
    std::string get_line() const { return line_; }

    // If any of the following strings is found at the start of a line, that line is skipped
    void set_comment_identifiers( const std::vector< std::string > & comment_identifiers ) { comment_identifiers_ = comment_identifiers; }

    bool skip_empty_lines() const { return skip_empty_lines_; }
    void set_skip_empty_lines( const bool skip_empty_lines ) { skip_empty_lines_ = skip_empty_lines; }

    bool skip_whitespace_only_lines() const { return skip_whitespace_only_lines_; }
    void set_skip_whitespace_only_lines( const bool skip_whitespace_only_lines ) { skip_whitespace_only_lines_ = skip_whitespace_only_lines; }

    bool allow_single_quotes() const { return allow_single_quotes_; }
    void set_allow_single_quotes( const bool allow_single_quotes ) { allow_single_quotes_ = allow_single_quotes; }

//    std::vector< std::string > words() const;

    // Causes the next read to return the current line again.
    // This can potentially lead to strange results if this function is called
    // after the last line has been read
    // Perhaps we should instead add a function "peek_next_line()",
    // which allows you to look at the following line without going to the next line on the next read.
    void push_back_last_line() const { push_back_last_line_ = true; }

private:
    std::ifstream input_file_;
    std::string line_;
    size_t line_number_;
    bool skip_empty_lines_;
    bool skip_whitespace_only_lines_;
    bool allow_single_quotes_; // Ugly name and quick hack to allow reading of .inp files without trying to interpret "'"
    std::vector< std::string > comment_identifiers_;
    mutable bool push_back_last_line_;
};

#endif // TEXTFILEREADER_H

