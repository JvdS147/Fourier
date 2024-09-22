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
#include "StringFunctions.h"
#include "Utilities.h"

#include <fstream>
#include <stdexcept>

#include <iostream> // For debugging

const char default_slash_character = '/';

// ********************************************************************************

FileName::FileName()
{
    slash_character_ = default_slash_character;
    set_directory( "" );
    set_file_name( "" );
    set_extension( "" );
}

// ********************************************************************************

FileName::FileName( std::string input )
{
    slash_character_ = default_slash_character;
    check_if_quotes_correct( input ); // @@ Has not been programmed yet
    input = strip( input );
    input = remove_delimiters( input, "\"", "\"" );
    input = replace( input, "/", "\\" );
    input = replace( input, "\\.\\", "\\" );
    // Replace all occurrences of "\directory\..\" by "\"
    bool changed( true );
    while ( changed )
    {
        changed = false;
        size_t iPos2 = input.find( "\\..\\" );
        if ( iPos2 != std::string::npos )
        {
            if ( iPos2 == 0 )
                throw std::runtime_error( "FileName::FileName( std::string ): directory starts with \\..\\." );
            // Find "\" starting at iPos2 and searching backwards
            size_t iPos1 = input.find_last_of( '\\', iPos2-1 );
            if ( iPos1 != std::string::npos )
            {
                input = input.substr( 0, iPos1 ) + input.substr( iPos2 + 3 );
                changed = true;
            }
        }
    }
    // Find the last occurrence of "\\"
    size_t iPos1 = input.find_last_of( "\\" );
    if ( iPos1 == std::string::npos )
    {
        iPos1 = 0;
        set_directory( "" );
    }
    else
    {
        ++iPos1;
        set_directory( "\"" + input.substr( 0, iPos1 ) + "\"" );
    }
    // Find the last occurrence of "."
    // Three problem cases here:
    // The last "." is in the directory name, not in the file_name
    // The file_name ends with "."
    // There is no "."
    size_t iPos2 = input.find_last_of( "." );
    if ( ( iPos2 == std::string::npos ) ||
         ( iPos2 < iPos1 ) )
    {
        set_file_name( "\"" + input.substr( iPos1 ) + "\"" );
        set_extension( "" );
    }
    else
    {
        set_file_name( "\"" + input.substr( iPos1, iPos2-iPos1 ) + "\"" );
        set_extension( "\"" + input.substr( iPos2+1 ) + "\"" );
    }
}

// ********************************************************************************

FileName::FileName( const std::string & directory, const std::string & file_name, const std::string & extension )
{
    slash_character_ = default_slash_character;
    set_directory( directory );
    set_file_name( file_name );
    set_extension( extension );
}

// ********************************************************************************

void FileName::set_full_name( const std::string & file_name )
{
    *this = FileName( file_name );
}
    
// ********************************************************************************

std::string FileName::full_name() const
{
    return assemble_file_name();
}

// ********************************************************************************

void FileName::set_directory( const std::string & directory )
{
    directory_ = strip( directory );
    directory_ = replace( directory_, "/", "\\" );
    directory_ = replace( directory_, "\\.\\", "\\" );
    if ( is_enclosed_in_quotes( directory_ ) )
        directory_ = extract_delimited_text( directory_, "\"", "\"" );
    if ( ! ( directory_.empty() || ( directory_.substr( directory_.length() - 1 ) == "\\" ) ) )
        directory_ += "\\";

    // @@ Currently, "\" \"" would give major problems

}

// ********************************************************************************

// Extension is NOT included
void FileName::set_file_name( const std::string & file_name )
{
    file_name_ = strip( file_name );
    file_name_ = replace( file_name_, "/", "\\" );
    if ( is_enclosed_in_quotes( file_name_ ) )
        file_name_ = extract_delimited_text( file_name_, "\"", "\"" );
    if ( file_name_.substr(0,1) == "\\" )
        file_name_ = file_name_.substr( 1 );
    size_t iPos = file_name_.find( "\\" );
    if ( iPos != std::string::npos )
        throw std::runtime_error( "FileName::set_file_name(): file name contains \\: " + file_name );

    //    @@ We should check here for other things that are not allowed in file names
        
    if ( file_name_.empty() )
        return;
    if ( file_name_.substr( file_name_.length()-1, 1 ) == "." )
        file_name_ = file_name.substr( 0, file_name_.length()-1 );
    if ( file_name_.empty() )
        return;
    if ( file_name_.substr( file_name_.length()-1, 1 ) == "." )
        throw std::runtime_error( "FileName::set_file_name(): file name ends in multiple dots." );
}

// ********************************************************************************

void FileName::set_extension( const std::string & extension )
{
    extension_ = strip( extension );
    if ( is_enclosed_in_quotes( extension_ ) )
        extension_ = extract_delimited_text( extension_, "\"", "\"" );
        // @@ Error here if length == 0???
    if ( extension_.substr( 0, 1 ) == "." )
        extension_ = extension_.substr( 1 );
    if ( extension_.find( "." ) != std::string::npos )
        throw std::runtime_error( "FileName::set_extension(): extension contains multiple dots." );
}

// ********************************************************************************

// Turns "C:\Data\file_name.txt" into "C:\\Data\\file_name.txt", necessary when writing input files for e.g. R
std::string FileName::escape_slashes() const
{
    return ::escape_slashes( full_name() );
}

// ********************************************************************************

// Technically I guess it checks if the file exists and is readable.
bool FileName::exists() const
{
    std::ifstream input_file( full_name().c_str() );
    if ( ! input_file )
       return false;
    input_file.close();
    return true;
}

// ********************************************************************************

char FileName::slash_character() const
{
    return slash_character_[0];
}

// ********************************************************************************

void FileName::set_slash_character( const char slash_character )
{
    if ( slash_character == '/' )
        slash_character_ = "/";
    else if ( slash_character == '\\' )
        slash_character_ = "\\";
    else
        throw std::runtime_error( "FileName::set_slash_character(): slash character must be / or \\." );
}

// ********************************************************************************

std::string FileName::assemble_file_name() const
{
    size_t iPos;
    bool contains_space( false );
    iPos = directory_.find( " " );
    if ( iPos != std::string::npos )
        contains_space = true;
    else
    {
        iPos = file_name_.find( " " );
        if ( iPos != std::string::npos )
            contains_space = true;
        else
        {
            iPos = extension_.find( " " );
            if ( iPos != std::string::npos )
                contains_space = true;
        }
    }
    std::string result( directory_ );
    result += file_name_;
    if ( ! extension_.empty() )
        result += "." + extension_;
    if ( contains_space )
        result = "\"" + result + "\"";
    return correct_slashes( result );
}

// ********************************************************************************

FileName replace_extension( const FileName & file_name, const std::string & new_extension )
{
    return FileName( file_name.directory(), file_name.file_name(), new_extension );
}

// ********************************************************************************

// append_to_file_name( "C:\directory\file.txt", "_2" ) => "C:\directory\file_2.txt"
FileName append_to_file_name( const FileName & file_name, const std::string & addition )
{
    return FileName( file_name.directory(), file_name.file_name() + addition, file_name.extension() );
}

// ********************************************************************************

std::string append_backslash( const std::string & input )
{
    if ( input.empty() || ( input.substr( input.length() - 1 ) == "\\" ) )
        return input;
    return input + "\\";
}

// ********************************************************************************

std::string FileName::correct_slashes( const std::string & input ) const
{
    return replace( input, "\\", slash_character_ );
}

// ********************************************************************************

FileName generate_unique_file_name( const FileName & file_name )
{
    size_t i( 1 );
    while ( append_to_file_name( file_name, "_" + size_t2string( i, 4, '0' ) ).exists() )
        ++i;
    return append_to_file_name( file_name, "_" + size_t2string( i, 4, '0' ) );
}

// ********************************************************************************

std::vector< FileName > sort_file_names_by_extension( int argc, char** argv, std::vector< std::string > extensions )
{
    std::string expected_extensions;
    for ( size_t i( 0 ); i != extensions.size(); ++i )
        expected_extensions += " " + extensions[i];
    if ( argc != extensions.size() + 1 )
        throw std::runtime_error( "sort_file_names_by_extension(): Error: expected" + expected_extensions + "." );
    std::vector< FileName > result;
    for ( size_t i( 0 ); i != extensions.size(); ++i )
        result.push_back( FileName( argv[ i+1 ] ) );
    for ( size_t i( 0 ); i != extensions.size(); ++i )
    {
        if ( to_upper( result[i].extension() ) == to_upper( extensions[i] ) )
            continue;
        bool match_found( false );
        for ( size_t j( i+1 ); j != extensions.size(); ++j )
        {
            if ( to_upper( result[j].extension() ) == to_upper( extensions[i] ) )
            {
                std::swap( result[i], result[j] );
                std::swap( extensions[i], extensions[j] );
                match_found = true;
                break;
            }
        }
        if ( ! match_found )
            throw std::runtime_error( "sort_file_names_by_extension(): Error: expected" + expected_extensions + "." );
    }
    return result;
}

// ********************************************************************************

