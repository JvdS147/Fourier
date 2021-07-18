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

#include "FileList.h"
#include "FileName.h"
#include "TextFileReader.h"
#include "TextFileWriter.h"
#include "Utilities.h" // @@ Only for append_backslash(), that is silly

// ********************************************************************************

FileList::FileList() : prepend_file_name_with_basedirectory_(true)
{
}

// ********************************************************************************

FileList::FileList( const FileName & file_name ) : prepend_file_name_with_basedirectory_(true)
{
    initialise_from_file( file_name );
}

// ********************************************************************************

FileList::FileList( const std::vector< FileName > & file_names ) : file_names_(file_names), prepend_file_name_with_basedirectory_(true)
{
}

// ********************************************************************************

FileList::FileList( const std::string & base_directory, const FileName & file_name ) : base_directory_(base_directory), prepend_file_name_with_basedirectory_(true)
{
    base_directory_ = append_backslash( base_directory_ );
    if ( file_name.directory().empty() )
        initialise_from_file_2( FileName( base_directory_, file_name.file_name(), file_name.extension() ) );
    else
        initialise_from_file_2( file_name );
}

// ********************************************************************************

FileList::FileList( const std::string & base_directory, const std::vector< FileName > & file_names ) : base_directory_(base_directory), file_names_(file_names), prepend_file_name_with_basedirectory_(true)
{
    base_directory_ = append_backslash( base_directory_ );
}

// ********************************************************************************

FileName FileList::value( const size_t i ) const
{
    if ( prepend_file_name_with_basedirectory_ && file_names_[ i ].directory().empty() )
        return FileName( base_directory_, file_names_[ i ].file_name(), file_names_[ i ].extension() );
    return file_names_[ i ];
}

// ********************************************************************************

void FileList::save( const FileName & file_name ) const
{
    TextFileWriter file_list_writer( file_name );
    for ( size_t i( 0 ); i != size(); ++i )
        file_list_writer.write_line( file_names_[ i ].full_name() );
}

// ********************************************************************************

void FileList::initialise_from_file( const FileName & file_name )
{
    
    initialise_from_file_2( file_name );
    base_directory_ = file_name.directory();
}

// ********************************************************************************

void FileList::initialise_from_file_2( const FileName & file_name )
{
    TextFileReader text_file_reader( file_name );
    text_file_reader.set_skip_empty_lines( true );
    std::vector< std::string > comment_identifiers;
    comment_identifiers.push_back( "#" );
    text_file_reader.set_comment_identifiers( comment_identifiers );
    std::vector< std::string > words;
    while ( text_file_reader.get_next_line( words ) )
    {
        if ( words.size() != 1 )
        {
            //std::cout << << std::endl;
        }
        else
            file_names_.push_back( FileName( words[0] ) );
    }
}

// ********************************************************************************

