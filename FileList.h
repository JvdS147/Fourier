#ifndef FILELIST_H
#define FILELIST_H

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

#include "FileName.h"

#include <string>
#include <vector>

/*
  Under Windows, there is no easy, portable way to get a list of files from the operating system.
  This should be done with boost.
  As a hand-cranked alternative, read a file containing the file names. The file can be generated trivially by a command like:

  dir /B *.cif > cif_files.txt

  Files with spaces in them must be enclosed in double quotes.
  
  If no base directory is supplied but the file name of the file list itself contains a path,
  then that path is taken to be the base directory.
  
  If a base directory is supplied and the file name of the file list itself contains no path,
  then the base directory is taken to be the path of the file name.
  
  This class is becoming a bit of a Frankenclass. It is doing too many things and the implementation is not up to date with the documentation.
*/
class FileList
{
public:

    FileList();

    // The directory of the file_name will also be applied to the file names in the file.
    explicit FileList( const FileName & file_name );

    explicit FileList( const std::vector< FileName > & file_names );

    // The base_directory refers to the file names in the list *and* to the file name of the file containing the list,
    // unless that file already has a directory
    FileList( const std::string & base_directory, const FileName & file_name );

    FileList( const std::string & base_directory, const std::vector< FileName > & file_names );

    void initialise_from_file( const FileName & file_name );

    void push_back( const FileName & file_name ) { file_names_.push_back( file_name ); }

    void erase( const size_t i );

    void reserve( const size_t nvalues ) { file_names_.reserve( nvalues ); }
    size_t size() const { return file_names_.size(); }
    bool empty() const { return file_names_.empty(); }

    bool prepend_file_name_with_basedirectory() const { return prepend_file_name_with_basedirectory_; }
    void set_prepend_file_name_with_basedirectory( const bool prepend_file_name_with_basedirectory ) { prepend_file_name_with_basedirectory_ = prepend_file_name_with_basedirectory; }

    // Splits the file names into n lists of equal size that are stored as basename_i.
    // This provides a very quick and dirty parallelisation mechanism.
    void split( const size_t n, const std::string & basename ) const;

    // Guaranteed to end in a backslash
    std::string base_directory() const { return base_directory_; }

    // This is wasteful and an iterator would be much more elegant
    FileName value( const size_t i ) const;

    // @@ This really needs a flag to indicate if the directory name should explicitly be included / deleted / left as-is
    void save( const FileName & file_name ) const;

private:
    std::string base_directory_;
    std::vector< FileName > file_names_;
    bool prepend_file_name_with_basedirectory_;

    void initialise_from_file_2( const FileName & file_name );
};

FileList merge( const FileList & lhs, const FileList & rhs );

#endif // FILELIST_H

