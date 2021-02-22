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

#include "ReadXSD.h"
#include "3DCalculations.h"
#include "CrystalStructure.h"
#include "TextFileReader.h"
#include "Utilities.h"

#include <iostream>
#include <stdexcept>

// ********************************************************************************

void read_xsd( const FileName & file_name, CrystalStructure & crystal_structure )
{

// This is what we are looking for:
//                        <Atom3d ID="42" HasProperties="1" Mapping="104" Parent="2" Children="442,443,444,445" Name="Cl39" UserID="39" XYZ="-0.544406997224045,-0.119468460872243,-0.582538283348663" Connections="92" TemperatureType="Isotropic" AnisotropicTemperature="0,0,0,0,0,0,0,0,0" Components="Cl" Visible="0">
//                            <Properties NMRShielding="563.67" EFGQuadrupolarCoupling="-72.78" EFGAsymmetry="0.22" Force="-0.214441360892478,1.32464677237851,1.59841007334931"/>
//                        </Atom3d>
//
// If the user only assigned the NMR shieldings and not the EFG etc., the lines look like this:
//                        <Atom3d ID="4" HasProperties="1" Mapping="135" Parent="2" Children="137,393,394,395" Name="C1" UserID="1" XYZ="0.004362191714,0.464087446874054,0.655927752421914" Connections="69,70,100" TemperatureType="Isotropic" AnisotropicTemperature="0,0,0,0,0,0,0,0,0" Components="C">
//                            <Properties NMRShielding="8.93"/>
//                        </Atom3d>
//
// Beware of lines of this type:
//                        <Atom3d ID="566" Mapping="952" ImageOf="4"/>
//

    TextFileReader text_file_reader( file_name );
    text_file_reader.set_skip_empty_lines( true );
    Splitter splitter_comma( "," );
    std::vector< std::string > words;
    std::vector< Atom > atoms;
    bool unit_cell_found( false );
    while ( text_file_reader.get_next_line( words ) )
    {
        if ( words[0] == "<Atom3d" )
        {
            if ( ( words.size() == 4 ) && ( words[3].substr( words[3].length() - 2 ) == "/>" ) )
                continue;
            std::string total_line = text_file_reader.get_line();
            std::string line;
            size_t line_counter( 1 );
            while ( text_file_reader.get_next_line( line ) )
            {
                ++line_counter;
                line = strip( line );
                if ( line == "</Atom3d>" )
                    break;
                else
                    total_line += " " + line;
            }
            if ( line_counter != 3 )
                std::cout << "WARNING: number of lines for one atom not equal three." << std::endl;
            words = split( total_line );

//		<Atom3d ID="42" HasProperties="1" Mapping="355" Parent="2" Name="H23" XYZ="-0.129725828307252,0.460163887977488,0.106616484219332" Connections="84" ForcefieldType="h1" Charge="0.053" TemperatureType="Isotropic" AnisotropicTemperature="0,0,0,0,0,0,0,0,0" Components="H" Visible="0">
//			<Properties Force="0.0051914800946029,-0.00371354061100781,-0.0309962205299153"/>
//		</Atom3d>

            // Find "Name="
            bool name_found( false );
            std::string label;
            for ( size_t i( 0 ); i != words.size(); ++i )
            {
                if ( words[i].substr( 0, 5 ) == "Name=" )
                {
                    label = words[i].substr( 6, words[i].length() - 7 );
                    name_found = true;
                    break;
                }
            }
            if ( ! name_found )
                throw std::runtime_error( "read_xsd(): name not found." );
            // Find "XYZ="
            bool xyz_found( false );
            Vector3D position;
            for ( size_t i( 0 ); i != words.size(); ++i )
            {
                if ( words[i].substr( 0, 4 ) == "XYZ=" )
                {
                    std::string xyz_string = extract_delimited_text( words[i], "\"", "\"" );
                    std::vector< std::string > words_2 = splitter_comma.split( xyz_string );
                    if ( words_2.size() != 3 )
                        throw std::runtime_error( "read_xsd(): xyz coordinates not found 2." );
                    position = Vector3D( string2double( words_2[0] ), string2double( words_2[1] ), string2double( words_2[2] ) );
                    xyz_found = true;
                    break;
                }
            }
            if ( ! xyz_found )
                throw std::runtime_error( "read_xsd(): xyz coordinates not found 1." );
            atoms.push_back( Atom( element_from_atom_label( label ), position, label ) );
        }
        else if ( words[0] == "<SpaceGroup" )
        {
            
//  <SpaceGroup ID="352" Parent="2" Children="356" Name="P1" DisplayStyle="Solid" Color="255,255,255,255"
// AVector="14.4767680982048,0,-8.07646527115057" BVector="0,5.24594832057182,0" CVector="0,0,16.3107460377921" OrientationBase="C along Z, B in YZ plane"
// Centering="3D Primitive-Centered" Lattice="3D Triclinic" GroupName="P1" Operators="1,0,0,0,0,1,0,0,0,0,1,0" DisplayRange="0,1,0,1,0,1" LineThickness="2"
// CylinderRadius="0.2" LabelAxes="1" ActiveSystem="2" ITNumber="1" LongName="P 1" Qualifier="Origin-1" SchoenfliesName="C1-1" System="Triclinic" Class="1"/>

            bool avector_found( false );
            Vector3D a_vector;
            for ( size_t i( 0 ); i != words.size(); ++i )
            {
                if ( words[i].substr( 0, 8 ) == "AVector=" )
                {
                    std::string xyz_string = extract_delimited_text( words[i], "\"", "\"" );
                    std::vector< std::string > words_2 = splitter_comma.split( xyz_string );
                    if ( words_2.size() != 3 )
                        throw std::runtime_error( "read_xsd(): AVector not found 2." );
                    a_vector = Vector3D( string2double( words_2[0] ), string2double( words_2[1] ), string2double( words_2[2] ) );
                    avector_found = true;
                    break;
                }
            }
            if ( ! avector_found )
                throw std::runtime_error( "read_xsd(): AVector not found 1." );
            bool bvector_found( false );
            Vector3D b_vector;
            for ( size_t i( 0 ); i != words.size(); ++i )
            {
                if ( words[i].substr( 0, 8 ) == "BVector=" )
                {
                    std::string xyz_string = extract_delimited_text( words[i], "\"", "\"" );
                    std::vector< std::string > words_2 = splitter_comma.split( xyz_string );
                    if ( words_2.size() != 3 )
                        throw std::runtime_error( "read_xsd(): BVector not found 2." );
                    b_vector = Vector3D( string2double( words_2[0] ), string2double( words_2[1] ), string2double( words_2[2] ) );
                    bvector_found = true;
                    break;
                }
            }
            if ( ! bvector_found )
                throw std::runtime_error( "read_xsd(): BVector not found 1." );
            bool cvector_found( false );
            Vector3D c_vector;
            for ( size_t i( 0 ); i != words.size(); ++i )
            {
                if ( words[i].substr( 0, 8 ) == "CVector=" )
                {
                    std::string xyz_string = extract_delimited_text( words[i], "\"", "\"" );
                    std::vector< std::string > words_2 = splitter_comma.split( xyz_string );
                    if ( words_2.size() != 3 )
                        throw std::runtime_error( "read_xsd(): CVector not found 2." );
                    c_vector = Vector3D( string2double( words_2[0] ), string2double( words_2[1] ), string2double( words_2[2] ) );
                    cvector_found = true;
                    break;
                }
            }
            if ( ! cvector_found )
                throw std::runtime_error( "read_xsd(): CVector not found 1." );
            crystal_structure.set_crystal_lattice( CrystalLattice( a_vector.length(), b_vector.length(), c_vector.length(), angle( b_vector, c_vector ), angle( a_vector, c_vector ), angle( a_vector, b_vector ) ) );
            unit_cell_found = true;
        }
    }
    if ( ! unit_cell_found )
        throw std::runtime_error( "read_xsd(): unit cell not found." );
    crystal_structure.add_atoms( atoms );
}

// ********************************************************************************

