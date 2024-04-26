#ifndef TEXTFILE_H
#define TEXTFILE_H

/* *********************************************
Copyright (c) 2013-2023, Cornelis Jan (Jacco) van de Streek
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

// Reads whole file and keeps it in memory
// Only suitable for small files!
// Pretty much entirely identical in behaviour to a std::vector< std::string >
// In keeping with C++ convention: zero-based.
// New lines:
// Windows       : \r\n
// Linux         : \n
// MacOS < 10    : \r
// MacOS >= 10   : \n
// C++ in principle uses (e.g. absorbs upon reading) \n
// This class deletes all \r characters from all input.

// @@ It would really be nice if this class interface and the interface of TextFileReader were essentially the same.
// @@ Would be nice to include some kind of filter, like an on-the-fly purge?

class TextFile
{
public:

    TextFile() {}

    explicit TextFile( const FileName & file_name );

    // For testing
    explicit TextFile( const std::vector< std::string > & lines );

    void read( const FileName & file_name );

    size_t size() const { return lines_.size(); }

    // In keeping with C++ convention: zero-based.
    std::string line( const size_t i ) const;

    // Deletes entire lines.
    // comment_identifier must start at the start of the line.
    void purge_comment_lines( std::string comment_identifier, const bool case_sensitive = true );

    // Empty means empty, a line with only whitespace is not empty
    void purge_empty_lines();

    // Starts search from line i
    // return a line number or std::string::npos
    // An empty string matches nothing
    size_t find( const std::string & word, const size_t i_start = 0 ) const;

    // Starts search from line i
    // return a line number or std::string::npos
    // An empty string matches nothing
    size_t find_whole_word( const std::string & word, const size_t i_start = 0 ) const;

private:
    std::vector< std::string > lines_;
};

#endif // TEXTFILE

