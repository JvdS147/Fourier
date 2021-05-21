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

#include "ReadCif.h"
#include "CheckFoundItem.h"
#include "CrystalStructure.h"
#include "TextFileReader.h"
#include "TextFileWriter.h"
#include "Utilities.h"

#include <stdexcept>

#include <iostream> // For debugging

namespace {

class AtomLineInterpreter
{
public:
    AtomLineInterpreter( const std::vector< std::string > & loop_items )
    {
        loop_items_size_ = loop_items.size();
        site_label_index_       = loop_items_size_;
        site_type_symbol_index_ = loop_items_size_;
        x_coordinate_index_     = loop_items_size_;
        y_coordinate_index_     = loop_items_size_;
        z_coordinate_index_     = loop_items_size_;
        Uiso_index_             = loop_items_size_;
        charge_index_           = loop_items_size_;
        occupancy_index_        = loop_items_size_;
        for ( size_t i( 0 ); i != loop_items_size_; ++i )
        {
            if ( loop_items[i] == "_atom_site_label" )
                site_label_index_ = i;
            else if ( loop_items[i] == "_atom_site_type_symbol" )
                site_type_symbol_index_ = i;
            else if ( loop_items[i] == "_atom_site_fract_x" )
                x_coordinate_index_ = i;
            else if ( loop_items[i] == "_atom_site_fract_y" )
                y_coordinate_index_ = i;
            else if ( loop_items[i] == "_atom_site_fract_z" )
                z_coordinate_index_ = i;
            else if ( loop_items[i] == "_atom_site_U_iso_or_equiv" )
                Uiso_index_ = i;
            else if ( loop_items[i] == "_atom_site_charge" )
                charge_index_ = i;
            else if ( loop_items[i] == "_atom_site_occupancy" )
                occupancy_index_ = i;
        }
        all_items_found_ = true;
        if ( x_coordinate_index_ == loop_items_size_ )
            all_items_found_ = false; // throw std::runtime_error( "read_cif(): _atom_site_fract_x missing from atom loop_." );
        if ( y_coordinate_index_ == loop_items_size_ )
            all_items_found_ = false; // throw std::runtime_error( "read_cif(): _atom_site_fract_y missing from atom loop_." );
        if ( z_coordinate_index_ == loop_items_size_ )
            all_items_found_ = false; // throw std::runtime_error( "read_cif(): _atom_site_fract_z missing from atom loop_." );
        if ( ( site_label_index_ == loop_items_size_ ) && ( site_type_symbol_index_ == loop_items_size_ ) )
            all_items_found_ = false; // throw std::runtime_error( "read_cif(): need one of _atom_site_label and _atom_site_type_symbol." );
    }
    
    void interpret( const std::vector< std::string > & words, CrystalStructure & crystal_structure ) const
    {
        if ( ! all_items_found_ )
            return;
        if ( words.size() != loop_items_size_ )
        {
//            std::cout << "loop_items_size_ = " << loop_items_size_ << std::endl;
//            for ( size_t i( 0 ); i != words.size(); ++i )
//            {
//                std::cout << words[i] << std::endl;
//            }
            throw std::runtime_error( "read_cif(): atom line must have same number of items as specified in loop." );
        }
        double x = string2double( words[x_coordinate_index_] );
        double y = string2double( words[y_coordinate_index_] );
        double z = string2double( words[z_coordinate_index_] );
        std::string label;
        if ( site_label_index_ != loop_items_size_ )
            label = words[site_label_index_];
        std::string element_string;
        if ( site_type_symbol_index_ != loop_items_size_ )
            element_string = words[site_type_symbol_index_];
        if ( site_label_index_ == loop_items_size_ )
            label = element_string;
        if ( site_type_symbol_index_ == loop_items_size_ )
        {
            if ( label.size() == 1 )
                element_string = label;
            else
            {
                if ( isalpha( label[1] ) )
                    element_string = label.substr( 0, 2 );
                else
                    element_string = label.substr( 0, 1 );
            }
        }
        Atom new_atom( Element( element_string ), Vector3D( x, y, z ), label );
        if ( charge_index_ != loop_items_size_ )
            new_atom.set_charge( string2double( words[charge_index_] ) );
        if ( Uiso_index_ != loop_items_size_ )
            new_atom.set_Uiso( string2double( words[Uiso_index_] ) );
        if ( occupancy_index_ != loop_items_size_ )
            new_atom.set_occupancy( string2double( words[occupancy_index_] ) );
        crystal_structure.add_atom( new_atom );
    }

private:
    size_t loop_items_size_;
    size_t site_label_index_;
    size_t site_type_symbol_index_;
    size_t x_coordinate_index_;
    size_t y_coordinate_index_;
    size_t z_coordinate_index_;
    size_t Uiso_index_;
    size_t charge_index_;
    size_t occupancy_index_;
    
    bool all_items_found_; // @@ major hack
};

class AnisoLineInterpreter
{
public:
    AnisoLineInterpreter( const std::vector< std::string > & loop_items )
    {
        loop_items_size_ = loop_items.size();
        label_index_ = loop_items_size_;
        U11_index_   = loop_items_size_;
        U22_index_   = loop_items_size_;
        U33_index_   = loop_items_size_;
        U12_index_   = loop_items_size_;
        U13_index_   = loop_items_size_;
        U23_index_   = loop_items_size_;
        for ( size_t i( 0 ); i != loop_items_size_; ++i )
        {
            if ( loop_items[i] == "_atom_site_aniso_label" )
                label_index_ = i;
            else if ( loop_items[i] == "_atom_site_aniso_U_11" )
                U11_index_ = i;
            else if ( loop_items[i] == "_atom_site_aniso_U_22" )
                U22_index_ = i;
            else if ( loop_items[i] == "_atom_site_aniso_U_33" )
                U33_index_ = i;
            else if ( loop_items[i] == "_atom_site_aniso_U_12" )
                U12_index_ = i;
            else if ( loop_items[i] == "_atom_site_aniso_U_13" )
                U13_index_ = i;
            else if ( loop_items[i] == "_atom_site_aniso_U_23" )
                U23_index_ = i;
        }
        if ( label_index_ == loop_items_size_ )
            throw std::runtime_error( "read_cif(): _atom_site_aniso_label missing from _atom_site_aniso loop_." );
        if ( U11_index_ == loop_items_size_ )
            throw std::runtime_error( "read_cif(): _atom_site_aniso_U_11 missing from _atom_site_aniso loop_." );
        if ( U22_index_ == loop_items_size_ )
            throw std::runtime_error( "read_cif(): _atom_site_aniso_U_22 missing from _atom_site_aniso loop_." );
        if ( U33_index_ == loop_items_size_ )
            throw std::runtime_error( "read_cif(): _atom_site_aniso_U_33 missing from _atom_site_aniso loop_." );
        if ( U12_index_ == loop_items_size_ )
            throw std::runtime_error( "read_cif(): _atom_site_aniso_U_12 missing from _atom_site_aniso loop_." );
        if ( U13_index_ == loop_items_size_ )
            throw std::runtime_error( "read_cif(): _atom_site_aniso_U_13 missing from _atom_site_aniso loop_." );
        if ( U23_index_ == loop_items_size_ )
            throw std::runtime_error( "read_cif(): _atom_site_aniso_U_23 missing from _atom_site_aniso loop_." );
    }

    void interpret( const std::vector< std::string > & words, CrystalStructure & crystal_structure ) const
    {
        if ( words.size() != loop_items_size_ )
        {
            std::cout << "loop_items_size_ = " << loop_items_size_ << std::endl;
            for ( size_t i( 0 ); i != words.size(); ++i )
            {
                std::cout << words[i] << std::endl;
            }
            throw std::runtime_error( "read_cif(): atom aniso line must have same number of items as specified in loop." );
        }
        double U11 = string2double( words[U11_index_] );
        double U22 = string2double( words[U22_index_] );
        double U33 = string2double( words[U33_index_] );
        double U12 = string2double( words[U12_index_] );
        double U13 = string2double( words[U13_index_] );
        double U23 = string2double( words[U23_index_] );
        SymmetricMatrix3D U_cif( U11, U22, U33, U12, U13, U23 );
        SymmetricMatrix3D U_cart = U_cif_2_U_cart( U_cif, crystal_structure.crystal_lattice() );
        AnisotropicDisplacementParameters adps = AnisotropicDisplacementParameters( U_cart );
        size_t i = crystal_structure.find_label( words[label_index_] );
        if ( i == crystal_structure.natoms() )
            throw std::runtime_error( "read_cif(): atom not found." );
        Atom new_atom = crystal_structure.atom( i );
        // Check if the atom already has ADPs?
        new_atom.set_anisotropic_displacement_parameters( adps );
        crystal_structure.set_atom( i, new_atom );
    }

private:
    size_t loop_items_size_;
    size_t label_index_;
    size_t U11_index_;
    size_t U22_index_;
    size_t U33_index_;
    size_t U12_index_;
    size_t U13_index_;
    size_t U23_index_;

};

// ********************************************************************************

void deal_with_atom_loop( TextFileReader & text_file_reader, const std::vector< std::string > & loop_items, CrystalStructure & crystal_structure )
{
// We need:
//_atom_site_label
//_atom_site_type_symbol
//_atom_site_fract_x
//_atom_site_fract_y
//_atom_site_fract_z

    AtomLineInterpreter atom_line_interpreter( loop_items );
    std::vector< std::string > words;
    do // Read the atoms
    {
        if ( ! text_file_reader.get_next_line( words ) )
            return;
        if ( ( words[0][0] == '_' ) || ( words[0] == "loop_" ) )
        {
            text_file_reader.push_back_last_line();
            return;
        }
        // When we are here, we have an atom line
        atom_line_interpreter.interpret( words, crystal_structure );
    } while ( true );
}

// ********************************************************************************

void deal_with_symmetry_loop( TextFileReader & text_file_reader, const std::vector< std::string > & loop_items, CrystalStructure & crystal_structure )
{
//loop_
//_symmetry_equiv_pos_site_id
//_symmetry_equiv_pos_as_xyz     // This definition has been superseded and is retained here only for archival purposes. Use instead '_space_group_symop_operation_xyz'
// Example: -y+x,-y,1/3+z
//1 x,y,z

    bool skip_whitespace_only_lines = text_file_reader.skip_whitespace_only_lines();
    text_file_reader.set_skip_whitespace_only_lines( true );
    size_t symmetry_equiv_pos_as_xyz_index = loop_items.size();
    for ( size_t i(0); i != loop_items.size(); ++i )
    {
        if ( loop_items[i] == "_symmetry_equiv_pos_as_xyz" )
            symmetry_equiv_pos_as_xyz_index = i;
    }
    if ( symmetry_equiv_pos_as_xyz_index == loop_items.size() )
        throw std::runtime_error( "read_cif(): _symmetry_equiv_pos_as_xyz must be present." );
    std::vector< SymmetryOperator > symmetry_operators;
    std::vector< std::string > words;
    bool finished( false );
    do // Read the symmetry operators
    {
        if ( ! text_file_reader.get_next_line( words ) )
            finished = true;
        else if ( ( words[0][0] == '_' ) || ( words[0] == "loop_" ) )
        {
            text_file_reader.push_back_last_line();
            finished = true;
        }
        else // When we are here, we have a symmetry line
        {
            if ( words.size() != loop_items.size() )
                throw std::runtime_error( "read_cif(): symmetry line must have same number of items as specified in loop." );
            SymmetryOperator symmetry_operator( words[symmetry_equiv_pos_as_xyz_index] );
            symmetry_operators.push_back( symmetry_operator );
        }
    } while ( ! finished );
    text_file_reader.set_skip_whitespace_only_lines( skip_whitespace_only_lines );
    SpaceGroup space_group( symmetry_operators, "P21/c" );
    crystal_structure.set_space_group( space_group );
}

// ********************************************************************************

void deal_with_aniso_loop( TextFileReader & text_file_reader, const std::vector< std::string > & loop_items, CrystalStructure & crystal_structure )
{
//    loop_
//    _atom_site_aniso_label
//    _atom_site_aniso_U_11
//    _atom_site_aniso_U_22
//    _atom_site_aniso_U_33
//    _atom_site_aniso_U_12
//    _atom_site_aniso_U_13
//    _atom_site_aniso_U_23
//    Cl1 0.1071(7) 0.1580(10) 0.0482(5) -0.0068(7) 0.0000(5) -0.0108(5)

    AnisoLineInterpreter aniso_line_interpreter( loop_items );
    std::vector< std::string > words;
    do // Read the ADPs
    {
        if ( ! text_file_reader.get_next_line( words ) )
            return;
        if ( ( words[0][0] == '_' ) || ( words[0] == "loop_" ) || ( words[0] == "#END" ) )
        {
            text_file_reader.push_back_last_line();
            return;
        }
        // When we are here, we have an aniso line
        aniso_line_interpreter.interpret( words, crystal_structure );
    } while ( true );
}

// ********************************************************************************

void get_loop_items( TextFileReader & text_file_reader, std::vector< std::string > & loop_items )
{
    std::vector< std::string > words;
    bool OK( true );
    do // Read the loop item keywords
    {
        if ( ! text_file_reader.get_next_line( words ) )
            throw std::runtime_error( "read_cif(): loop empty." );
        if ( words[0][0] == '_' )
        {
            if ( words.size() != 1 )
                throw std::runtime_error( "read_cif(): loop_ item cannot have value." );
            loop_items.push_back( words[0] );
        }
        else
        {
            text_file_reader.push_back_last_line();
            return;
        }
    } while ( OK );
}

// ********************************************************************************

void deal_with_loop( TextFileReader & text_file_reader, CrystalStructure & crystal_structure )
{
    std::vector< std::string > words;
    if ( ! text_file_reader.get_next_line( words ) )
        throw std::runtime_error( "read_cif(): loop empty." );
    if ( words[0].length() < 5 )
        throw std::runtime_error( "read_cif(): unrecognised loop keyword." );
    if ( words[0][0] != '_' )
        throw std::runtime_error( "read_cif(): no keyword after loop keyword." );
    text_file_reader.push_back_last_line();
    std::vector< std::string > loop_items;
    get_loop_items( text_file_reader, loop_items );
    if ( words[0].substr( 0, 16 ) == "_atom_site_aniso" )
    {
        deal_with_aniso_loop( text_file_reader, loop_items, crystal_structure );
        return;
    }
    if ( words[0].substr( 0, 5 ) == "_atom" )
    {
        deal_with_atom_loop( text_file_reader, loop_items, crystal_structure );
        return;
    }
    if ( words[0].substr( 0, 9 ) == "_symmetry" )
    {
        deal_with_symmetry_loop( text_file_reader, loop_items, crystal_structure );
        return;
    }
    // For the moment we ignore everything else
}

// ********************************************************************************

} // namespace

//data_MD
//_symmetry_cell_setting           triclinic
//_symmetry_space_group_name_H-M   'P 1'
//_symmetry_Int_Tables_number      1
//loop_
//_symmetry_equiv_pos_site_id
//_symmetry_equiv_pos_as_xyz
//1 x,y,z
//_cell_length_a                   50.0
//_cell_length_b                   50.0
//_cell_length_c                   50.0
//_cell_angle_alpha                90
//_cell_angle_beta                 90
//_cell_angle_gamma                90
//loop_
//_atom_site_label
//_atom_site_type_symbol
//_atom_site_fract_x
//_atom_site_fract_y
//_atom_site_fract_z
//C C 0.024    0.15    0.04
//C C 0.099    0.05    0.04
//C C 0.027    0.10    0.13
//C C 0.093    0.00    0.13
//C C 0.215    0.05    0.04
//C C 0.206    0.00    0.13
//C C 0.148    0.10    0.13
//C C 0.158    0.15    0.04
//C C 0.090    0.20    0.13
//C C 0.299    0.15    0.04

// ********************************************************************************

// A very simple cif reader, can essentially only read cifs from MD trajectories
// from Materials Studio.
void read_cif( const FileName & file_name, CrystalStructure & crystal_structure )
{
    crystal_structure = CrystalStructure();
    TextFileReader text_file_reader( file_name );
    text_file_reader.set_skip_empty_lines( true ); // This is crucial
    std::vector< std::string > comment_identifiers;
    comment_identifiers.push_back( "#" );
    text_file_reader.set_comment_identifiers( comment_identifiers );
    std::string name;
    std::string space_group_str;
    bool found_a( false );
    bool found_b( false );
    bool found_c( false );
    bool found_alpha( false );
    bool found_beta( false );
    bool found_gamma( false );
    double a( 0.0 ); // Stupid initialisation to silence compiler warnings
    double b( 0.0 ); // Stupid initialisation to silence compiler warnings
    double c( 0.0 ); // Stupid initialisation to silence compiler warnings
    Angle alpha;
    Angle beta;
    Angle gamma;
    std::vector< std::string > words;
    while ( text_file_reader.get_next_line( words ) )
    {

        // According to the cif standard, data_ and loop_ are case-insensitive.

        if ( ( words[0].length() > 4 ) && ( to_lower( words[0].substr( 0, 5 ) ) == "data_" ) )
        {
            if ( words.size() != 1 )
                throw std::runtime_error( "read_cif(): data_ keyword cannot have a value." );
            name = words[0].substr( 5, std::string::npos );
            continue;
        }
        if ( to_lower( words[0] == "loop_" ) )
        {
            if ( words.size() != 1 )
                throw std::runtime_error( "read_cif(): loop_ keyword cannot have a value." );
            deal_with_loop( text_file_reader, crystal_structure );
            continue;
        }
        if ( words.size() == 2 )
        {
            if ( words[0] == "_symmetry_space_group_name_H-M" )
            {
                space_group_str = words[1];
                continue;
            }
            if ( words[0] == "_cell_length_a" )
            {
                a = string2double( words[1] );
                found_a = true;
                continue;
            }
            if ( words[0] == "_cell_length_b" )
            {
                b = string2double( words[1] );
                found_b = true;
                continue;
            }
            if ( words[0] == "_cell_length_c" )
            {
                c = string2double( words[1] );
                found_c = true;
                continue;
            }
            if ( words[0] == "_cell_angle_alpha" )
            {
                alpha = Angle::from_degrees( string2double( words[1] ) );
                found_alpha = true;
                continue;
            }
            if ( words[0] == "_cell_angle_beta" )
            {
                beta = Angle::from_degrees( string2double( words[1] ) );
                found_beta = true;
                continue;
            }
            if ( words[0] == "_cell_angle_gamma" )
            {
                gamma = Angle::from_degrees( string2double( words[1] ) );
                found_gamma = true;
                continue;
            }
        }
    }
    if ( ! ( found_a && found_b && found_c && found_alpha && found_beta && found_gamma ) )
        throw std::runtime_error( "read_cif(): not all cell parameters found." );
    CrystalLattice crystal_lattice( a, b, c, alpha, beta, gamma );
    crystal_structure.set_crystal_lattice( crystal_lattice );
    crystal_structure.set_name( name );
    if ( space_group_str != "" )
    {
        SpaceGroup space_group = crystal_structure.space_group();
        space_group.set_name( space_group_str );
        crystal_structure.set_space_group( space_group );
    }
}

// ********************************************************************************

// Entirely text based: removes all lines with five fields or more of which the first field starts with H, the second field is "H" and the third fourth and fifth field are floating point numbers
void remove_hydrogen_atoms( const FileName & input_file_name, const FileName & output_file_name )
{
    TextFileReader text_file_reader( input_file_name );
    TextFileWriter output_file( output_file_name );
    
    std::vector< std::string > words;
    while ( text_file_reader.get_next_line( words ) )
    {
        bool keep( true );
        if ( words.size() > 4 )
        {
            if ( ( to_upper( words[1] ) == "H" ) ||
                 ( to_upper( words[1] ) == "D" ) )
            {
                std::string element_string;
                if ( words[0].size() == 1 )
                    element_string = words[0];
                else
                {
                    if ( isalpha(words[0][1]) )
                        element_string = words[0].substr( 0, 2 );
                    else
                        element_string = words[0].substr( 0, 1 );
                }
                if ( ( to_upper( element_string ) == "H" ) ||
                     ( to_upper( element_string ) == "D" ) )
                    keep = false;
            }
        }
        if ( keep )
            output_file.write_line( text_file_reader.get_line() );
    }
}

// ********************************************************************************

// As above. The output file_name is the input file name with ".cif" replaced by "_noH.cif"
void remove_hydrogen_atoms( const FileName & input_file_name )
{
    remove_hydrogen_atoms( input_file_name, append_to_file_name( input_file_name, "_noH" ) );
}

// ********************************************************************************
