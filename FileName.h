#ifndef FILENAME_H
#define FILENAME_H

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

/*
 Takes care of lots of stuff:
 If it finds a space in one of the components, encloses output in ".
 Takes care of \ after directory
 Takes care of dot before extension, erases dot if extension is empty.
 On input, directory may or may not end in a \.
 On input, file_name may or may not end in a dot.
 On input, extension may or may not begin with a dot.
 On input, extension may or may not be empty, or even just a dot.
 On input, file_name may or may not include the extension (but in that case extension must be empty).
 None of the components is allowed to start or end with whitespace.
 The file_name and the extension are not allowed to contain a \ or /.

Replaces all / by \
Replaces all \.\ by \
Replace all \directory\..\ by \

The following is currently not checked:
Is it "C:some_directory" or "C:\some_directory" ?

If we find "X:", where X is some letter, at the start of the string, can we always replace it by "X:\" ?

*/
class FileName
{
public:

    FileName();
    
    explicit FileName( std::string file_name );

    FileName( const std::string & directory, const std::string & file_name, const std::string & extension );

    // Always ends in backslash.
    std::string directory() const { return directory_; }

    // Extension is NOT included, cannot end in a dot.
    std::string file_name() const { return file_name_; }

// @@ We need a file_name_plus_extension() const;

    // Does not include dot.
    std::string extension() const { return extension_; }
    
    void set_directory( const std::string & directory );

    // Extension is NOT included.
    void set_file_name( const std::string & file_name );

    void set_extension( const std::string & extension );

    void set_full_name( const std::string & file_name );

  // @@ Could add "force_quotes()" because some applications need that. Perhaps better, add "add_quotes()" to Utilities.h

// @@ currently duplicate of assemble_file_name()
    std::string full_name() const;

    // Outputs full name with escaped slashes, i.e. turns "C:\Data\file_name.txt" into "C:\\Data\\file_name.txt",
    // necessary when writing input files for e.g. R.
    std::string escape_slashes() const;

    // Technically I guess it checks if the file exists and is readable.
    bool exists() const;

    char slash_character() const;
    
    void set_slash_character( const char slash_character );
    
private:

    std::string directory_; // Always ends in backslash
    std::string file_name_; // Extension is NOT included
    std::string extension_; // Does not include dot
    std::string slash_character_;
    
    std::string assemble_file_name() const;
    std::string correct_slashes( const std::string & input ) const;
};

FileName replace_extension( const FileName & file_name, const std::string & new_extension );

// append_to_file_name( "C:/directory/file.txt", "_2" ) => "C:/directory/file_2.txt"
FileName append_to_file_name( const FileName & file_name, const std::string & addition );

std::string append_backslash( const std::string & input );

//std::string change_all_slashes_to( std::string & input,  );

// Uses FileName::exists() and append_to_file_name( FileName ) to generate file names like
// directory/file_name_0001.txt
// The directory must be specified otherwise the check with FileName::exists() is meaningless.
// At the moment:
// 1. Even if directory/file_name.txt would be available, directory/file_name_0001.txt (or whatever is available) is returned
// 2. The length of the "0001" is not configurable.
// If the number is greater than "9999", then as many positions as needed are used, i.e. the number is never corrupted, i.e.:
// directory/file_name_1234567.txt (not directory/file_name_1234.txt)
// The smallest free number is returned. In a way this could be inefficient, and the method could be supplied with a hint as to where
// to start searching.
FileName generate_unique_file_name( const FileName & file_name );

// This is to sort file names that are given as arguments.
// Stable sort. Case insensitive.
std::vector< FileName > sort_file_names_by_extension( int argc, char** argv, std::vector< std::string > extensions );

#endif // FILENAME_H
