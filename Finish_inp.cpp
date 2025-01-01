/* *********************************************
Copyright (c) 2013-2025, Cornelis Jan (Jacco) van de Streek
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

#include "Finish_inp.h"
#include "CopyTextFile.h"
#include "Element.h"
#include "FileName.h"
#include "StringFunctions.h"
#include "TextFileReader.h"
#include "TextFileWriter.h"
#include "TOPAS.h"

#include <iostream> // For debugging
#include <stdexcept>

// ********************************************************************************

void finish_inp( const FileName & input_file_name )
{
    { // Scoping brackets for the TextFileReader and the TextFileWriter, to ensure files are closed before being copied.
    TextFileReader text_file_reader( input_file_name );
    text_file_reader.set_skip_empty_lines( true );
    text_file_reader.set_allow_single_quotes( true );
    TextFileWriter text_file_writer( append_to_file_name( input_file_name, "_2" ) );
    size_t max_label_length = 0;
    Splitter splitter( "," );
    splitter.set_merge_delimiters( false );
    std::vector< std::string > words;
    while ( text_file_reader.get_next_line( words ) )
    {
        
        if ( words[0] == "#include" )
        {
        }
        else if ( ( words[0] == "do_errors" ) || ( words[0] == "'do_errors" ) )
        {
            text_file_writer.write_line( "'do_errors_include_penalties" );
            text_file_writer.write_line( "prm JvdS_shift = Get(refine_ls_shift_on_su_max); : 0");
            text_file_writer.write_line( "prm JvdS_numpar = Get(number_independent_parameters); : 0");
        }
        else if ( words[0] == "weighted_Durbin_Watson" )
        {
        }
        else if ( words[0] == "gof" )
        {
            text_file_writer.write_line( text_file_reader.get_line() );
            text_file_writer.write_line( "'continue_after_convergence" );
        }
        else if ( words[0] == "Zero_Error(" )
        {
            text_file_writer.write_line( insert_at_sign( text_file_reader.get_line() ) );
        }
        else if ( words[0] == "LP_Factor(" )
        {
            text_file_writer.write_line( "'Synchrotron use: LP_Factor( 90 )" );
            text_file_writer.write_line( "'Neutrons use: LP_Factor( 90 )" );
            text_file_writer.write_line( "'No monochromator use: LP_Factor( 0 )" );
            text_file_writer.write_line( "'Ge Monochromator, Cu radiation, use LP_Factor( 27.3 )" );
            text_file_writer.write_line( "'Graphite Monochromator, Cu radiation, use LP_Factor( 26.4 )" );
            text_file_writer.write_line( "'Quartz Monochromator, Cu radiation, use LP_Factor( 26.6 )" );
            text_file_writer.write_line( text_file_reader.get_line() );
        }
        else if ( words[0] == "filament_length" )
        {
            text_file_writer.write_line( "    " + words[0] + " @ " + words[1] );
        }
        else if ( words[0] == "sample_length" )
        {
            text_file_writer.write_line( "    " + words[0] + " @ " + words[1] );
        }
        else if ( words[0] == "receiving_slit_length" )
        {
            text_file_writer.write_line( "    " + words[0] + " @ " + words[1] );
        }
        else if ( ( words.size() > 1 ) && ( words[1] == "217.5" ) )
        {
        }
        else if ( words[0] == "axial_n_beta" )
        {
            text_file_writer.write_line( "    axial_n_beta 50" );
        }
        else if ( words[0] == "axial_del" )
        {
        }
        else if ( words[0] == "str" )
        {
            text_file_writer.write_line( text_file_reader.get_line() );
            text_file_writer.write_line( "    r_bragg 0.0" );
        }
        else if ( words[0] == "CS_L(" )
        {
            text_file_writer.write_line( insert_at_sign( text_file_reader.get_line() ) );
        }
        else if ( words[0] == "CS_G(" )
        {
            text_file_writer.write_line( insert_at_sign( text_file_reader.get_line() ) );
        }
        else if ( words[0] == "Strain_G(" )
        {
            text_file_writer.write_line( insert_at_sign( text_file_reader.get_line() ) );
        }
        else if ( words[0] == "Strain_L(" )
        {
            text_file_writer.write_line( insert_at_sign( text_file_reader.get_line() ) );
        }
        else if ( ( words.size() > 1 ) && ( words[1] == "!sh_scale" ) )
        {
            text_file_writer.write_line( text_file_reader.get_line() );
        }
        else if ( ( words.size() > 1 ) && ( words[1] == "ref_flag" ) )
        {
            text_file_writer.write_line( "    macro ref_flag {   }" );
        }
        else if ( words[0] == "site" )
        {
            if ( words[1].length() > max_label_length )
                max_label_length = words[1].length();
            text_file_writer.write_line( text_file_reader.get_line() );
        }
        else if ( ( words.size() > 1 ) && ( words[1] == "!flatten_weight" ) )
        {
            text_file_writer.write_line( "    prm !flatten_weight    100000" );
        }
        else if ( words[0].substr(0,17) == "Distance_Restrain" )
        {
            // Distance_Restrain(N1      C2     ,   1.47566, 1.47769`, bond_width, bond_weight)
            std::string line = extract_delimited_text( text_file_reader.get_line(), "(", ")" );
            std::vector< std::string > words_2 = splitter.split( line );
            std::vector< std::string > words_3 = split( words_2[0] );
            if ( words_3.size() != 2 )
                throw std::runtime_error( "Programming error." );
            std::string new_line = "    Distance_Restrain(";
            for ( size_t i( 0 ); i != words_3.size(); ++i )
                new_line += " " + pad( words_3[i], max_label_length, ' ' );
            if ( ( element_from_atom_label( words_3[0] ).atomic_number() == 1 ) ||
                 ( element_from_atom_label( words_3[1] ).atomic_number() == 1 ) )
                words_2[1] = "0.95000";
            for ( size_t i( 1 ); i != words_2.size(); ++i )
                new_line += ", " + strip( words_2[i] );
            new_line += " )";
            text_file_writer.write_line( new_line );
        }
        else if ( words[0].substr(0,14) == "Angle_Restrain" )
        {
            // Angle_Restrain(C2      N1      C3     , 112.42979, 112.44385, angle_width, angle_weight)
            std::string line = extract_delimited_text( text_file_reader.get_line(), "(", ")" );
            std::vector< std::string > words_2 = splitter.split( line );
            std::vector< std::string > words_3 = split( words_2[0] );
            if ( words_3.size() != 3 )
                throw std::runtime_error( "Programming error." );
            std::string new_line = "    Angle_Restrain(";
            for ( size_t i( 0 ); i != words_3.size(); ++i )
                new_line += " " + pad( words_3[i], max_label_length, ' ' );
            for ( size_t i( 1 ); i != words_2.size(); ++i )
                new_line += ", " + strip( words_2[i] );
            new_line += " )";
            text_file_writer.write_line( new_line );
        }
        else if ( words[0].substr(0,7) == "Flatten" )
        {
            // Flatten( C7      C2      C15     O24     C23     H34     C32    ,, 87.2962215, flatten_width, flatten_weight)
            std::string line = extract_delimited_text( text_file_reader.get_line(), "(", ")" );
            std::vector< std::string > words_2 = splitter.split( line );
            std::vector< std::string > words_3 = split( words_2[0] );
            std::string new_line = "    Flatten(";
            for ( size_t i( 0 ); i != words_3.size(); ++i )
                new_line += " " + pad( words_3[i], max_label_length, ' ' );
            for ( size_t i( 1 ); i != words_2.size(); ++i )
                new_line += ", " + strip( words_2[i] );
            new_line += " )";
            text_file_writer.write_line( new_line );
        }
        else if ( words[0].substr(0,16) == "Out_CIF_STR_Uiso" ) // Change Out_CIF_STR_Uiso("") to Out_CIF_STR( "" )
        {
            std::string line = extract_delimited_text( text_file_reader.get_line(), "(", ")" );
            FileName file_name_FN( line );
            std::string new_file_name = FileName( file_name_FN.directory(), file_name_FN.name(), "cif" ).full_name();
            if ( ! is_enclosed_in_quotes( new_file_name ) )
                new_file_name = "\"" + new_file_name + "\"";
            text_file_writer.write_line( "    Out_CIF_STR( " + new_file_name + " )" );
        }
        else if ( words[0].substr(0,11) == "Out_Profile" )
        {
            // Out_Profile("C:\GD\name_plot.pro")
            // Keep everything between ( and )
            std::string line = extract_delimited_text( text_file_reader.get_line(), "(", ")" );
            FileName file_name_FN( line );
            std::string file_name = file_name_FN.name();
            if ( file_name.substr( file_name.length() - 5, 5 ) == "_plot" )
                file_name = file_name.substr( 0, file_name.length() - 4 ) + "profile";
            std::string new_file_name = FileName( file_name_FN.directory(), file_name, "txt" ).full_name();
            if ( ! is_enclosed_in_quotes( new_file_name ) )
                new_file_name = "\"" + new_file_name + "\"";
            text_file_writer.write_line( "    xdd_out " + new_file_name + " load out_record out_fmt out_eqn" );
            text_file_writer.write_line( "    {" );
            text_file_writer.write_line( "        \" %11.5f \" = X;" );
            text_file_writer.write_line( "        \" %11.5f \" = Yobs;" );
            text_file_writer.write_line( "        \" %11.5f \" = Ycalc;" );
            text_file_writer.write_line( "        \" %11.5f\\n\" = SigmaYobs;" );
            text_file_writer.write_line( "    }" );
        }
        else if ( words[0].substr(0,8) == "Out_Tick" )
        {
            // Out_Tick("C:\GD\name_plot.tic")
            std::string line = extract_delimited_text( text_file_reader.get_line(), "(", ")" );
            FileName file_name_FN( line );
            std::string file_name = file_name_FN.name();
            std::string fcf_file_name;
            if ( file_name.substr( file_name.length() - 5, 5 ) == "_plot" )
            {
                fcf_file_name = file_name.substr( 0, file_name.length() - 5 );
                file_name = file_name.substr( 0, file_name.length() - 4 ) + "tickmarks";
            }
            else
                fcf_file_name = file_name;
            std::string directory = file_name_FN.directory();
            FileName new_file_name_FN( directory, file_name, "txt" );
            std::string new_file_name = new_file_name_FN.full_name();
            if ( ! is_enclosed_in_quotes( new_file_name ) )
                new_file_name = "\"" + new_file_name + "\"";
            text_file_writer.write_line( "    phase_out " + new_file_name + " load out_record out_fmt out_eqn" );
            text_file_writer.write_line( "    {" );
            text_file_writer.write_line( "        \" %11.5f -200\\n\" = 2.0 * Rad * Th;" );
            text_file_writer.write_line( "    }" );
            fcf_file_name = FileName( directory, fcf_file_name, "fcf" ).full_name();
            if ( ! is_enclosed_in_quotes( fcf_file_name ) )
                fcf_file_name = "\"" + fcf_file_name + "\"";
            text_file_writer.write_line( "    Out_FCF( " + fcf_file_name + " )" );
        }
        else if ( words[0].substr(0,19) == "Out_PowderDataBlock" )
        {
        }
        else
            text_file_writer.write_line( text_file_reader.get_line() );
    }
    } // Scoping brackets for the TextFileReader and the TextFileWriter, to ensure files are closed before being copied.
    copy_text_file( input_file_name, append_to_file_name( input_file_name, "_old" ) );
    copy_text_file( append_to_file_name( input_file_name, "_2" ), input_file_name );
}

// ********************************************************************************

