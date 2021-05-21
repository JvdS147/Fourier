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

#include "TLSWriter.h"
#include "CopyTextFile.h"
#include "CrystalStructure.h"
#include "ReadCif.h"
#include "TextFileReader_2.h"
#include "TextFileWriter.h"
#include "TOPAS.h"
#include "Utilities.h"

#include <iostream>
#include <stdexcept>

// ********************************************************************************

// From .cif file
void TLSWriter( const FileName & input_file_name )
{
    bool is_on_inversion_at_origin( true );
    TextFileWriter text_file_writer( replace_extension( append_to_file_name( input_file_name, "_TLS" ), "inp" ) );
    bool write_bond_restraints = FileName( input_file_name.directory(), input_file_name.file_name() + "_bond_restraints", "txt" ).exists();
    bool write_angle_restraints = FileName( input_file_name.directory(), input_file_name.file_name() + "_angle_restraints", "txt" ).exists();
    std::vector< std::string > bond_labels_1;
    std::vector< std::string > bond_labels_2;
    std::vector< double > bond_target_values;
    std::vector< std::string > angle_labels_1;
    std::vector< std::string > angle_labels_2;
    std::vector< std::string > angle_labels_3;
    std::vector< double > angle_target_values;
    if ( write_bond_restraints )
    {
        std::cout << "Bond restraints file found, bond restraints will be written out." << std::endl;
        TextFileReader_2 file_bond_restraints( FileName( input_file_name.directory(), input_file_name.file_name() + "_bond_restraints", "txt" ) );
        std::vector< std::string > words;
        for ( size_t i( 1 ); i != file_bond_restraints.size(); ++i )
        {
            words = split( file_bond_restraints.line( i ) );
            if ( words.size() != 4 )
                throw std::runtime_error( "TLS .inp writer: unexpected format in bond restraints file: >" + file_bond_restraints.line( i ) + "<" );
            bond_labels_1.push_back( words[1] );
            bond_labels_2.push_back( words[2] );
            bond_target_values.push_back( string2double( words[3] ) );
        }
    }
    else
        std::cout << "No bond restraints file found, no bond restraints will be written out." << std::endl;
    if ( write_angle_restraints )
    {
        std::cout << "Angle restraints file found, angle restraints will be written out." << std::endl;
        TextFileReader_2 file_angle_restraints( FileName( input_file_name.directory(), input_file_name.file_name() + "_angle_restraints", "txt" ) );
        std::vector< std::string > words;
        for ( size_t i( 1 ); i != file_angle_restraints.size(); ++i )
        {
            words = split( file_angle_restraints.line( i ) );
            if ( words.size() != 5 )
                throw std::runtime_error( "TLS .inp writer: unexpected format in angle restraints file: >" + file_angle_restraints.line( i ) + "<" );
            angle_labels_1.push_back( words[1] );
            angle_labels_2.push_back( words[2] );
            angle_labels_3.push_back( words[3] );
            angle_target_values.push_back( string2double( words[4] ) );
        }
    }
    else
        std::cout << "No angle restraints file found, no angle restraints will be written out." << std::endl;
    // We must find the unique combinations, e.g.:
    // angle r1-r2-r3 and bond r3-r2 only require two combinations: r1-r2 and r3-r2
    // Create list of unique bonds
    std::vector< std::string > unique_labels_1;
    std::vector< std::string > unique_labels_2;
    for ( size_t i( 0 ); i != angle_labels_1.size(); ++i )
    {
        bool found( false );
        // Find bond 1
        for ( size_t j( 0 ); j != unique_labels_1.size(); ++j )
        {
            if ( ( unique_labels_1[j] == angle_labels_1[i] ) && ( unique_labels_2[j] == angle_labels_2[i] ) )
            {
                 found = true;
                 break;
            }
        }
        if ( ! found )
        {
            unique_labels_1.push_back( angle_labels_1[i] );
            unique_labels_2.push_back( angle_labels_2[i] );
        }
        found = false;
        // Find bond 2
        for ( size_t j( 0 ); j != unique_labels_1.size(); ++j )
        {
            if ( ( unique_labels_1[j] == angle_labels_3[i] ) && ( unique_labels_2[j] == angle_labels_2[i] ) )
            {
                 found = true;
                 break;
            }
        }
        if ( ! found )
        {
            unique_labels_1.push_back( angle_labels_3[i] );
            unique_labels_2.push_back( angle_labels_2[i] );
        }
    }
//    for ( size_t i( 0 ); i != bond_labels_1.size(); ++i )
//    {
//        bool found( false );
//        // Find bond
//        for ( size_t j( 0 ); j != unique_labels_1.size(); ++j )
//        {
//            if ( ( unique_labels_1[j] == bond_labels_1[i] ) && ( unique_labels_2[j] == bond_labels_2[i] ) )
//            {
//                 found = true;
//                 break;
//            }
//        }
//        if ( ! found )
//        {
//            unique_labels_1.push_back( bond_labels_1[i] );
//            unique_labels_2.push_back( bond_labels_2[i] );
//        }
//        found = false;
//
//    }
    CrystalStructure crystal_structure;
    std::cout << "Now reading cif... " + input_file_name.full_name() << std::endl;
    read_cif( input_file_name, crystal_structure );
    if ( ! is_on_inversion_at_origin )
    {
        text_file_writer.write_line( "    ' Origin of the molecule (fractional coordinates)." );
        Vector3D c_o_m = crystal_structure.centre_of_mass();
        text_file_writer.write_line( "    prm centx " + double2string( c_o_m.x() ) );
        text_file_writer.write_line( "    prm centy " + double2string( c_o_m.y() ) );
        text_file_writer.write_line( "    prm centz " + double2string( c_o_m.z() ) );
        text_file_writer.write_line();
    }
    text_file_writer.write_line( "    ' T11 etc. are elements of the T, L and S tensors." );
    text_file_writer.write_line( "    prm T11  0.05" );
    text_file_writer.write_line( "    prm T22  0.05" );
    text_file_writer.write_line( "    prm T33  0.05" );
    text_file_writer.write_line( "    prm T12  0.0" );
    text_file_writer.write_line( "    prm T13  0.0" );
    text_file_writer.write_line( "    prm T23  0.0" );
    text_file_writer.write_line( "    prm L11  0.0" );
    text_file_writer.write_line( "    prm L22  0.0" );
    text_file_writer.write_line( "    prm L33  0.0" );
    text_file_writer.write_line( "    prm L12  0.0" );
    text_file_writer.write_line( "    prm L13  0.0" );
    text_file_writer.write_line( "    prm L23  0.0" );
    if ( ! is_on_inversion_at_origin )
    {
        text_file_writer.write_line( "    prm S11  0.0" );
        text_file_writer.write_line( "    prm S12  0.0" );
        text_file_writer.write_line( "    prm S13  0.0" );
        text_file_writer.write_line( "    prm S22  0.0" );
        text_file_writer.write_line( "    prm S23  0.0" );
        text_file_writer.write_line();
        text_file_writer.write_line( "    ' Choice of origin to make S symmetric (Downs, p80-81)." );
        text_file_writer.write_line( "    prm S21 = S12;" );
        text_file_writer.write_line( "    prm S31 = S13;" );
        text_file_writer.write_line( "    prm S32 = S23;" );
        text_file_writer.write_line();
        text_file_writer.write_line( "    ' Standard constraint to make equations determinate." );
        text_file_writer.write_line( "    prm S33 = -S11 - S22; : 0.0" );
    }
    text_file_writer.write_line();
    text_file_writer.write_line( "    ' The libration is sqrt(trace(L))." );
    text_file_writer.write_line( "    prm libration = Sqrt(L11 + L22 + L33) * Rad; : 0.0" );
    text_file_writer.write_line();
    text_file_writer.write_line( "    ' Conversion of fractional coordinates to Cartesian." );
    text_file_writer.write_line( "    ' xN1   = fractional x coordinate of atom N1" );
    text_file_writer.write_line( "    ' rxN1  = Cartesian coordinate x of atom N1. 'r' always indicates a Cartesian coordinate." );
    text_file_writer.write_line( "    ' drxN1 = Correction for libration to Cartesian x coordinate of atom N1." );
    text_file_writer.write_line( "    ' crxN1 = Cartesian x coordinate of atom N1 corrected for libration:" );
    text_file_writer.write_line( "    ' crxN1 = rxN1 + drxN1" );
    text_file_writer.write_line();
    text_file_writer.write_line( "    ' Unit-cell parameters are assumed to be constant." );
    text_file_writer.write_line();
    Matrix3D f2c = crystal_structure.crystal_lattice().fractional_to_orthogonal_matrix();

   // // Many matrix elements evaluate to e.g. -1.0E-18, these are set to 0.0 here.
  //  for ( size_t i( 0 ); i < 3; ++i )
  //  {
  //      for ( size_t j( 0 ); j < 3; ++j )
  //      {
   //         if ( nearly_zero( f2c.value( i, j ) ) )
   //             f2c.set_value( i, j, 0.0 );
   //     }
   // }

    if ( is_on_inversion_at_origin )
    {
        for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
         {
            std::string label = crystal_structure.atom( i ).label();
            std::string line;
            line = "    prm !rx" + label + " = x" + label + "*" + double2string( f2c.value(0,0) );
            if ( ! nearly_zero( f2c.value(0,1) ) )
                 line += " + y" + label + "*" + double2string( f2c.value(0,1) );
            if ( ! nearly_zero( f2c.value(0,2) ) )
                 line += " + z" + label + "*" + double2string( f2c.value(0,2) );
            line += ";";
            text_file_writer.write_line( line );
            line = "    prm !ry" + label + " = y" + label + "*" + double2string( f2c.value(1,1) );
            if ( ! nearly_zero( f2c.value(1,2) ) )
                 line += " + z" + label + "*" + double2string( f2c.value(1,2) );
            line += ";";
            text_file_writer.write_line( line );
            line = "    prm !rz" + label + " = z" + label + "*" + double2string( f2c.value(2,2) );
            line += ";";
            text_file_writer.write_line( line );
        }
    }
    else
    {
        for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
        {
            std::string label = crystal_structure.atom( i ).label();
            text_file_writer.write_line( "    prm !rx" + label + " = (x" + label + " - centx)*" + double2string( f2c.value(0,0) ) + " + (y" + label + " - centy)*" + double2string( f2c.value(0,1) ) + " + (z" + label + " - centz)*" + double2string( f2c.value(0,2) ) + ";" );
            text_file_writer.write_line( "    prm !ry" + label + " = (y" + label + " - centy)*" + double2string( f2c.value(1,1) ) + " + (z" + label + " - centz)*" + double2string( f2c.value(1,2) ) + ";" );
            text_file_writer.write_line( "    prm !rz" + label + " = (z" + label + " - centz)*" + double2string( f2c.value(2,2) ) + ";" );
        }
    }
    text_file_writer.write_line();
    text_file_writer.write_line( "    ' Calculation of Uij values in Cartesian framework. See Dunitz page 251, Eqn 5.42 and Table 5.1 page 252." );
    text_file_writer.write_line();
    if ( is_on_inversion_at_origin )
    {
        for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
        {
            std::string label = crystal_structure.atom( i ).label();
            text_file_writer.write_line( "    prm ru11" + label + " = L22*rz" + label + "^2 + L33*ry" + label + "^2 - 2*ry" + label + "*rz" + label + "*L23 + T11; : 0.0" );
            text_file_writer.write_line( "    prm ru22" + label + " = L11*rz" + label + "^2 + L33*rx" + label + "^2 - 2*rx" + label + "*rz" + label + "*L13 + T22; : 0.0" );
            text_file_writer.write_line( "    prm ru33" + label + " = L11*ry" + label + "^2 + L22*rx" + label + "^2 - 2*rx" + label + "*ry" + label + "*L12 + T33; : 0.0" );
            text_file_writer.write_line( "    prm ru12" + label + " = -rx" + label + "*ry" + label + "*L33 + rx" + label + "*rz" + label + "*L23 + ry" + label + "*rz" + label + "*L13 - rz" + label + "^2*L12 + T12; : 0.0" );
            text_file_writer.write_line( "    prm ru13" + label + " = -rx" + label + "*rz" + label + "*L22 + rx" + label + "*ry" + label + "*L23 - ry" + label + "^2*L13 + ry" + label + "*rz" + label + "*L12 + T13; : 0.0" );
            text_file_writer.write_line( "    prm ru23" + label + " = -ry" + label + "*rz" + label + "*L11 - rx" + label + "^2*L23 + rx" + label + "*ry" + label + "*L13 + rx" + label + "*rz" + label + "*L12 + T23; : 0.0" );
            text_file_writer.write_line();
        }
    }
    else
    {
        for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
        {
            std::string label = crystal_structure.atom( i ).label();
            text_file_writer.write_line( "    prm ru11" + label + " = L22*rz" + label + "^2 + L33*ry" + label + "^2 - 2*ry" + label + "*rz" + label + "*L23 - 2*ry" + label + "*S31 + 2*rz" + label + "*S21 + T11; : 0.0" );
            text_file_writer.write_line( "    prm ru22" + label + " = L11*rz" + label + "^2 + L33*rx" + label + "^2 - 2*rx" + label + "*rz" + label + "*L13 - 2*rz" + label + "*S12 + 2*rx" + label + "*S32 + T22; : 0.0" );
            text_file_writer.write_line( "    prm ru33" + label + " = L11*ry" + label + "^2 + L22*rx" + label + "^2 - 2*rx" + label + "*ry" + label + "*L12 - 2*rx" + label + "*S23 + 2*ry" + label + "*S13 + T33; : 0.0" );
            text_file_writer.write_line( "    prm ru12" + label + " = -rx" + label + "*ry" + label + "*L33 + rx" + label + "*rz" + label + "*L23 + ry" + label + "*rz" + label + "*L13 - rz" + label + "^2*L12 - rz" + label +
                                     "*S11 + rz" + label + "*S22 + rx" + label + "*S31 - ry" + label + "*S32 + T12; : 0.0" );
            text_file_writer.write_line( "    prm ru13" + label + " = -rx" + label + "*rz" + label + "*L22 + rx" + label + "*ry" + label + "*L23 - ry" + label + "^2*L13 + ry" + label + "*rz" + label + "*L12 + ry" + label +
                                     "*S11 - ry" + label + "*S33 + rz" + label + "*S23 - rx" + label + "*S21 + T13; : 0.0" );
            text_file_writer.write_line( "    prm ru23" + label + " = -ry" + label + "*rz" + label + "*L11 - rx" + label + "^2*L23 + rx" + label + "*ry" + label + "*L13 + rx" + label + "*rz" + label + "*L12 - rx" + label +
                                     "*S22 + rx" + label + "*S33 + ry" + label + "*S12 - rz" + label + "*S13 + T23; : 0.0" );
            text_file_writer.write_line();
        }
    }
    text_file_writer.write_line( "    ' ru11N1 = u11 of atom N1 in Cartesian coordinates." );
    text_file_writer.write_line( "    ' u11N1 = u11 of atom N1 in cif coordinates." );
    text_file_writer.write_line( "    ' See R. W. Grosse-Kunstleve & P. D. Adams (2002). J. Appl. Cryst. 35, 477-480." );
    text_file_writer.write_line();
    for ( size_t i( 0 ); i < crystal_structure.natoms(); ++i )
    {
        std::string label = crystal_structure.atom( i ).label();
        std::vector<double> k;
        k.push_back( 1.0 / crystal_structure.crystal_lattice().a_star() );
        k.push_back( 1.0 / crystal_structure.crystal_lattice().b_star() );
        k.push_back( 1.0 / crystal_structure.crystal_lattice().c_star() );
        Matrix3D c2f = crystal_structure.crystal_lattice().orthogonal_to_fractional_matrix();
        size_t r;
        size_t s;
        // Many matrix elements evaluate to e.g. -1.0E-18, these are set to 0.0 here.
        for ( r = 0; r < 3; ++r )
        {
            for ( s = 0; s < 3; ++s )
            {
                if ( nearly_zero( c2f.value( r, s ) ) )
                    c2f.set_value( r, s, 0.0 );
            }
        }
        // TOPAS cannot cope with "prm n = 4 + -2;", a workaround is to write "prm n = 4 + (-2);".
        std::vector< size_t > r_s; // "r's", multiple of "r"
        std::vector< size_t > s_s; // "s's", multiple of "s"
        r_s.push_back( 1 );
        s_s.push_back( 1 );
        r_s.push_back( 2 );
        s_s.push_back( 2 );
        r_s.push_back( 3 );
        s_s.push_back( 3 );
        r_s.push_back( 1 );
        s_s.push_back( 2 );
        r_s.push_back( 1 );
        s_s.push_back( 3 );
        r_s.push_back( 2 );
        s_s.push_back( 3 );


//    prm u11H13 = 37.9275 * (
//        0.160488 * ( ru11H13 * 0.160488 + ru12H13 * 0 + ru13H13 * 0.0246915 ) +
//        (0) * ( ru12H13 * 0.160488 + ru22H13 * 0 + ru23H13 * 0.0246915 ) +
//        (0.0246915) * ( ru13H13 * 0.160488 + ru23H13 * 0 + ru33H13 * 0.0246915 ) ); : 0.0


        for ( size_t j(0); j != r_s.size(); ++j )
        {
            r = r_s[j];
            s = s_s[j];
            text_file_writer.write_line( "    prm u" + size_t2string(r) + size_t2string(s) + label + " = " + double2string( k[r-1] * k[s-1] ) + " * ( " );
            text_file_writer.write_line( "        " + double2string( c2f.value(r-1,0) ) + " * ( ru11" + label + " * " + double2string( c2f.value(s-1,0) ) + " + ru12" + label + " * " + double2string( c2f.value(s-1,1) ) + " + ru13" + label + " * " + double2string( c2f.value(s-1,2) ) + " ) + " );
            text_file_writer.write_line( "        (" + double2string( c2f.value(r-1,1) ) + ") * ( ru12" + label + " * " + double2string( c2f.value(s-1,0) ) + " + ru22" + label + " * " + double2string( c2f.value(s-1,1) ) + " + ru23" + label + " * " + double2string( c2f.value(s-1,2) ) + " ) + " );
            text_file_writer.write_line( "        (" + double2string( c2f.value(r-1,2) ) + ") * ( ru13" + label + " * " + double2string( c2f.value(s-1,0) ) + " + ru23" + label + " * " + double2string( c2f.value(s-1,1) ) + " + ru33" + label + " * " + double2string( c2f.value(s-1,2) ) + " )" + " ); : 0.0" );
        }
        text_file_writer.write_line();
    }

    for ( size_t i( 0 ); i < crystal_structure.natoms(); ++i )
    {
        std::string label = crystal_structure.atom( i ).label();
        text_file_writer.write_line( "    site " + label + " x = x" + label + ";" + " y = y" + label + ";" + " z = z" + label + ";" +
                                     " occ " + crystal_structure.atom( i ).element().symbol() + " 1 ADPs_Keep_PD u11=u11" + label + "; u22=u22" + label + "; u33=u33" + label + "; u12=u12" + label + "; u13=u13" + label + "; u23=u23" + label + ";" );
    }
    text_file_writer.write_line();
    text_file_writer.write_line( "    ' Calculate corrected Cartesian coordinates. See Downs, page 83." );
    text_file_writer.write_line( "    ' Note that only intramolecular corrected distances can be trusted, the intermolecular corrected distances have no physical meaning." );
    for ( size_t i( 0 ); i < crystal_structure.natoms(); ++i )
    {
        std::string label = crystal_structure.atom( i ).label();
        text_file_writer.write_line( "    prm !drx" + label + " = 0.5*( (L22+L33)*rx" + label + "        -L12*ry" + label + "        -L13*rz" + label + " ) ; : 0.0" );
        text_file_writer.write_line( "    prm !dry" + label + " = 0.5*(      -L12*rx" + label + " + (L11+L33)*ry" + label + "        -L23*rz" + label + " ) ; : 0.0" );
        text_file_writer.write_line( "    prm !drz" + label + " = 0.5*(      -L13*rx" + label + "        -L23*ry" + label + " + (L11+L22)*rz" + label + " ) ; : 0.0" );
        text_file_writer.write_line( "    prm !crx" + label + " =                 rx" + label + " +          drx" + label + " ; : 0.0" );
        text_file_writer.write_line( "    prm !cry" + label + " =                 ry" + label + " +          dry" + label + " ; : 0.0" );
        text_file_writer.write_line( "    prm !crz" + label + " =                 rz" + label + " +          drz" + label + " ; : 0.0" );
    }
    text_file_writer.write_line();
    if ( write_bond_restraints )
    {
        text_file_writer.write_line( "    ' Calculate current bond lengths." );
        for ( size_t i( 0 ); i < bond_labels_1.size(); ++i )
        {
            text_file_writer.write_line( "    prm !curr_dist_" + bond_labels_1[i] + "_" + bond_labels_2[i] + " = Sqrt( (rx" + bond_labels_1[i] + "-rx" + bond_labels_2[i] + ")^2 + " +
                                                                                                                      "(ry" + bond_labels_1[i] + "-ry" + bond_labels_2[i] + ")^2 + " +
                                                                                                                      "(rz" + bond_labels_1[i] + "-rz" + bond_labels_2[i] + ")^2 ); : 0.0" );
        }
        text_file_writer.write_line();
        text_file_writer.write_line( "    ' Calculate corrected relative Cartesian coordinates. See Downs, page 83." );
        for ( size_t i( 0 ); i < bond_labels_1.size(); ++i )
        {
            text_file_writer.write_line( "    prm !corr_dist_" + bond_labels_1[i] + "_" + bond_labels_2[i] + " = Sqrt( (crx" + bond_labels_1[i] + "-crx" + bond_labels_2[i] + ")^2 + " +
                                                                                                                      "(cry" + bond_labels_1[i] + "-cry" + bond_labels_2[i] + ")^2 + " +
                                                                                                                      "(crz" + bond_labels_1[i] + "-crz" + bond_labels_2[i] + ")^2 ); : 0.0" );
        }
        text_file_writer.write_line();
        text_file_writer.write_line( "    prm !bond_width   0" );
        text_file_writer.write_line( "    prm !bond_weight  10000" );
        text_file_writer.write_line();
        text_file_writer.write_line( "    ' Restrain the corrected bond lengths." );
        for ( size_t i( 0 ); i < bond_labels_1.size(); ++i )
        {
            text_file_writer.write_line( "    Angle_Distance_Restrain( corr_dist_" + bond_labels_1[i] + "_" + bond_labels_2[i] + ", " + double2string( bond_target_values[i] ) + ", 0.0, bond_width, bond_weight )" );
        }
        text_file_writer.write_line();
    }
    if ( write_angle_restraints )
    {

        text_file_writer.write_line( "    ' Calculate corrected valence angles." );
        for ( size_t i( 0 ); i < unique_labels_1.size(); ++i )
        {
            text_file_writer.write_line( "    prm !diff_x_" + unique_labels_1[i] + "_" + unique_labels_2[i] + " = crx" + unique_labels_1[i] + " - crx" + unique_labels_2[i] + "; : 0.0" );
            text_file_writer.write_line( "    prm !diff_y_" + unique_labels_1[i] + "_" + unique_labels_2[i] + " = cry" + unique_labels_1[i] + " - cry" + unique_labels_2[i] + "; : 0.0" );
            text_file_writer.write_line( "    prm !diff_z_" + unique_labels_1[i] + "_" + unique_labels_2[i] + " = crz" + unique_labels_1[i] + " - crz" + unique_labels_2[i] + "; : 0.0" );
        }
        for ( size_t i( 0 ); i < angle_labels_1.size(); ++i )
        {
            text_file_writer.write_line( "    prm !corr_ang_" + angle_labels_1[i] + "_" + angle_labels_2[i] + "_" + angle_labels_3[i] + " = Rad * ArcCos( ( diff_x_" + angle_labels_1[i] + "_" + angle_labels_2[i] + " * diff_x_" + angle_labels_3[i] + "_" + angle_labels_2[i] + " + " +
                                                                                                                                                           "diff_y_" + angle_labels_1[i] + "_" + angle_labels_2[i] + " * diff_y_" + angle_labels_3[i] + "_" + angle_labels_2[i] + " + " +
                                                                                                                                                           "diff_z_" + angle_labels_1[i] + "_" + angle_labels_2[i] + " * diff_z_" + angle_labels_3[i] + "_" + angle_labels_2[i] + " ) /" +
                                                                                                                                                    " (Sqrt(diff_x_" + angle_labels_1[i] + "_" + angle_labels_2[i] + "^2+" +
                                                                                                                                                           "diff_y_" + angle_labels_1[i] + "_" + angle_labels_2[i] + "^2+" +
                                                                                                                                                           "diff_z_" + angle_labels_1[i] + "_" + angle_labels_2[i] + "^2) * " +
                                                                                                                                                      "Sqrt(diff_x_" + angle_labels_3[i] + "_" + angle_labels_2[i] + "^2+" +
                                                                                                                                                           "diff_y_" + angle_labels_3[i] + "_" + angle_labels_2[i] + "^2+" +
                                                                                                                                                           "diff_z_" + angle_labels_3[i] + "_" + angle_labels_2[i] + "^2)) ); : 0.0" );
        }
        text_file_writer.write_line();
        text_file_writer.write_line( "    prm !angle_width  0" );
        text_file_writer.write_line( "    prm !angle_weight 1" );
        text_file_writer.write_line();
        text_file_writer.write_line( "    ' Restrain the corrected valence angles." );
        for ( size_t i( 0 ); i < angle_labels_1.size(); ++i )
        {
           text_file_writer.write_line( "    Angle_Distance_Restrain( corr_ang_" + angle_labels_1[i] + "_" + angle_labels_2[i] + "_" + angle_labels_3[i] + ", " + double2string( angle_target_values[i] ) + ", 0.0, angle_width, angle_weight )" );
        }
    }

    text_file_writer.write_line();
    text_file_writer.write_line( "    Out_CIF_ADPs_TLS( \"filename.cif\" )" );
    text_file_writer.write_line();
    text_file_writer.write_line( "    xdd_out \"filename_profile.txt\" load out_record out_fmt out_eqn" );
    text_file_writer.write_line( "    {" );
    text_file_writer.write_line( "        \" %11.5f \" = X;" );
    text_file_writer.write_line( "        \" %11.5f \" = Yobs;" );
    text_file_writer.write_line( "        \" %11.5f \" = Ycalc;" );
    text_file_writer.write_line( "        \" %11.5f\\n\" = SigmaYobs;" );
    text_file_writer.write_line( "    }" );
    text_file_writer.write_line();
    text_file_writer.write_line( "    phase_out \"filename_tickmarks.txt\" load out_record out_fmt out_eqn" );
    text_file_writer.write_line( "    {" );
    text_file_writer.write_line( "        \" %11.5f -200\\n\" = 2.0 * Rad * Th;" );
    text_file_writer.write_line( "    }" );
    text_file_writer.write_line();

    text_file_writer.write_line( "    macro Out_CIF_ADPs_TLS( file )" );
    text_file_writer.write_line( "    {" );
    text_file_writer.write_line( "        out file" );
    text_file_writer.write_line( "        Out_String(\"data_\\n\")" );
    text_file_writer.write_line( "        Out(Get(a), \"_cell_length_a  %V\")" );
    text_file_writer.write_line( "        Out(Get(b), \"\\n_cell_length_b  %V\")" );
    text_file_writer.write_line( "        Out(Get(c), \"\\n_cell_length_c  %V\")" );
    text_file_writer.write_line( "        Out(Get(al), \"\\n_cell_angle_alpha %V\")" );
    text_file_writer.write_line( "        Out(Get(be), \"\\n_cell_angle_beta  %V\")" );
    text_file_writer.write_line( "        Out(Get(ga), \"\\n_cell_angle_gamma %V\")" );
    text_file_writer.write_line( "        Out(Get(cell_volume), \"\\n_cell_volume %V\")" );
    text_file_writer.write_line( "        Out(Get(sp_grp_char), \"\\n_symmetry_space_group_name_H-M %s\")" );
    text_file_writer.write_line( "        Out_String(\"\\nloop_\")" );
    text_file_writer.write_line( "        Out_String(\"\\n_symmetry_equiv_pos_as_xyz\")" );
    text_file_writer.write_line( "            Out(Get(sp_xyzs_txt),  \"%s\")" );
    text_file_writer.write_line( "        Out_String(\"\\nloop_\")" );
    text_file_writer.write_line( "        Out_String(\"\\n_atom_site_label\")" );
    text_file_writer.write_line( "        Out_String(\"\\n_atom_site_type_symbol\")" );
    text_file_writer.write_line( "        Out_String(\"\\n_atom_site_symmetry_multiplicity\")" );
    text_file_writer.write_line( "        Out_String(\"\\n_atom_site_fract_x\")" );
    text_file_writer.write_line( "        Out_String(\"\\n_atom_site_fract_y\")" );
    text_file_writer.write_line( "        Out_String(\"\\n_atom_site_fract_z\")" );
    text_file_writer.write_line( "        Out_String(\"\\n_atom_site_occupancy\")" );
    text_file_writer.write_line( "        atom_out file append" );
    text_file_writer.write_line( "            load out_record out_fmt out_eqn" );
    text_file_writer.write_line( "            {" );
    text_file_writer.write_line( "                \"\\n%s\" = Get_From_String(Get(current_atom), site);" );
    text_file_writer.write_line( "                \" %s\" = Get_From_String(Get(current_atom), atom);" );
    text_file_writer.write_line( "                \" %3.0f\" = Get_From_String(Get(current_atom), num_posns);" );
    text_file_writer.write_line( "                \" %V\" = Get_From_String(Get(current_atom), x);" );
    text_file_writer.write_line( "                \" %V\" = Get_From_String(Get(current_atom), y);" );
    text_file_writer.write_line( "                \" %V\" = Get_From_String(Get(current_atom), z);" );
    text_file_writer.write_line( "                \" %V\" = Get_From_String(Get(current_atom), occ);" );
    text_file_writer.write_line( "            }" );
    text_file_writer.write_line( "        out file append" );
    text_file_writer.write_line( "        Out_String(\"\\nloop_\")" );
    text_file_writer.write_line( "        Out_String(\"\\n_atom_site_aniso_label\")" );
    text_file_writer.write_line( "        Out_String(\"\\n_atom_site_aniso_U_11\")" );
    text_file_writer.write_line( "        Out_String(\"\\n_atom_site_aniso_U_22\")" );
    text_file_writer.write_line( "        Out_String(\"\\n_atom_site_aniso_U_33\")" );
    text_file_writer.write_line( "        Out_String(\"\\n_atom_site_aniso_U_12\")" );
    text_file_writer.write_line( "        Out_String(\"\\n_atom_site_aniso_U_13\")" );
    text_file_writer.write_line( "        Out_String(\"\\n_atom_site_aniso_U_23\")" );
    text_file_writer.write_line( "        atom_out file append adps" );
    text_file_writer.write_line( "            load out_record out_fmt out_eqn" );
    text_file_writer.write_line( "            {" );
    text_file_writer.write_line( "                \"\\n%s\" = Get_From_String(Get(current_atom), site);" );
    text_file_writer.write_line( "                \" %V\" = Get_From_String(Get(current_atom), u11);" );
    text_file_writer.write_line( "                \" %V\" = Get_From_String(Get(current_atom), u22);" );
    text_file_writer.write_line( "                \" %V\" = Get_From_String(Get(current_atom), u33);" );
    text_file_writer.write_line( "                \" %V\" = Get_From_String(Get(current_atom), u12);" );
    text_file_writer.write_line( "                \" %V\" = Get_From_String(Get(current_atom), u13);" );
    text_file_writer.write_line( "                \" %V\" = Get_From_String(Get(current_atom), u23);" );
    text_file_writer.write_line( "            }" );
    text_file_writer.write_line( "        out file append" );
    text_file_writer.write_line( "        Out_String(\"\\n\")" );
    text_file_writer.write_line( "    }" );
    copy_text_file( replace_extension( append_to_file_name( input_file_name, "_TLS" ), "inp" ), replace_extension( append_to_file_name( input_file_name, "_TLS" ), "org" ) );
}

// ********************************************************************************

// From .inp file
void TLSWriter_2( const FileName & input_file_name )
{
    TextFileReader_2 input_file( input_file_name );
    TextFileWriter text_file_writer( append_to_file_name( input_file_name, "_TLS" ) );
    text_file_writer.write_line( "/* TLS refinement." );
    text_file_writer.write_line( "Based on code by Simon Parsons, University of Edinburgh." );
    text_file_writer.write_line( "Generalised by Jacco van de Streek, University of Frankfurt." );
    text_file_writer.write_line();
    text_file_writer.write_line( "References are" );
    text_file_writer.write_line();
    text_file_writer.write_line( "R.T. Downs. High Temperature and High Pressure Crystal Chemistry. Rev. In Mineralogy & Geochemistry. Vol 41. Pages 61-87." );
    text_file_writer.write_line( "This excellent paper can be downloaded from http://www.geo.arizona.edu/xtal/group/pdf/chapter7.pdf" );
    text_file_writer.write_line();
    text_file_writer.write_line( "J. Dunitz: X-Ray Analysis and the Structures of Organic Molecules. Ch. 5" );
    text_file_writer.write_line( "This book has been republished by Wiley, and contains a nice introduction" );
    text_file_writer.write_line( "to TLS analysis applied to molecular systems. Highly recommended." );
    text_file_writer.write_line();
    text_file_writer.write_line( "*/" );
    bool found( false );
    size_t iLine( 0 );
    std::vector< std::string > words;
    while ( ( ! found ) && ( iLine != input_file.size() ) )
    {
        words = split_2( input_file.line( iLine ) );
        if ( words.size() == 2 )
        {
            if ( words[0] == "a" )
            {
                try
                {
                    /* double dummy = */ TOPASstring2double( words[1] );
                    found = true;
                }
                catch ( std::exception & e ) {}
            }
        }
        else if ( words.size() == 3 )
        {
            if ( ( words[0] == "a" ) && ( words[1] == "@" ) )
            {
                try
                {
                    /* double dummy = */ TOPASstring2double( words[2] );
                    found = true;
                }
                catch ( std::exception & e ) {}
            }
        }
        if ( ! found )
        {
            text_file_writer.write_line( input_file.line( iLine ) );
            ++iLine;
        }
    }
    // When we are here, either we are at the end of the file or we are at the line "a @ 8.743".
    if ( ! found )
        throw std::runtime_error( "TLSWriter_2(): could not find lattice parameter a." );
    // Assume that the following lines are a, b, c, al, be and ga
    double a( 0.0 ); // Stupid initialisation to silence compiler warnings
    double b( 0.0 );
    double c( 0.0 );
    double al( 0.0 );
    double be( 0.0 );
    double ga( 0.0 );
    words = split_2( input_file.line( iLine ) );
    if ( words.size() == 2 )
    {
        if ( words[0] == "a" )
            a = TOPASstring2double( words[1] );
    }
    else if ( words.size() == 3 )
    {
        if ( ( words[0] == "a" ) && ( words[1] == "@" ) )
            a = TOPASstring2double( words[2] );
    }
    else
        throw std::runtime_error( "TLSWriter_2(): could not find lattice parameters." );
    text_file_writer.write_line( input_file.line( iLine ) );
    ++iLine;
    words = split_2( input_file.line( iLine ) );
    if ( words.size() == 2 )
    {
        if ( words[0] == "b" )
            b = TOPASstring2double( words[1] );
    }
    else if ( words.size() == 3 )
    {
        if ( ( words[0] == "b" ) && ( words[1] == "@" ) )
            b = TOPASstring2double( words[2] );
    }
    else
        throw std::runtime_error( "TLSWriter_2(): could not find lattice parameters." );
    text_file_writer.write_line( input_file.line( iLine ) );

    ++iLine;
    words = split_2( input_file.line( iLine ) );
    if ( words.size() == 2 )
    {
        if ( words[0] == "c" )
            c = TOPASstring2double( words[1] );
    }
    else if ( words.size() == 3 )
    {
        if ( ( words[0] == "c" ) && ( words[1] == "@" ) )
            c = TOPASstring2double( words[2] );
    }
    else
        throw std::runtime_error( "TLSWriter_2(): could not find lattice parameters." );
    text_file_writer.write_line( input_file.line( iLine ) );

    ++iLine;
    words = split_2( input_file.line( iLine ) );
    if ( words.size() == 2 )
    {
        if ( words[0] == "al" )
            al = TOPASstring2double( words[1] );
    }
    else if ( words.size() == 3 )
    {
        if ( ( words[0] == "al" ) && ( words[1] == "@" ) )
            al = TOPASstring2double( words[2] );
    }
    else
        throw std::runtime_error( "TLSWriter_2(): could not find lattice parameters." );
    text_file_writer.write_line( input_file.line( iLine ) );

    ++iLine;
    words = split_2( input_file.line( iLine ) );
    if ( words.size() == 2 )
    {
        if ( words[0] == "be" )
            be = TOPASstring2double( words[1] );
    }
    else if ( words.size() == 3 )
    {
        if ( ( words[0] == "be" ) && ( words[1] == "@" ) )
            be = TOPASstring2double( words[2] );
    }
    else
        throw std::runtime_error( "TLSWriter_2(): could not find lattice parameters." );
    text_file_writer.write_line( input_file.line( iLine ) );

    ++iLine;
    words = split_2( input_file.line( iLine ) );
    if ( words.size() == 2 )
    {
        if ( words[0] == "ga" )
            ga = TOPASstring2double( words[1] );
    }
    else if ( words.size() == 3 )
    {
        if ( ( words[0] == "ga" ) && ( words[1] == "@" ) )
            ga = TOPASstring2double( words[2] );
    }
    else
        throw std::runtime_error( "TLSWriter_2(): could not find lattice parameters." );
    text_file_writer.write_line( input_file.line( iLine ) );
    // Read the lines between the lattice parameters and the "site"s
    ++iLine;
    found = false;
    while ( ( ! found ) && ( iLine != input_file.size() ) )
    {
        words = split_2( input_file.line( iLine ) );
        if ( ( words.size() > 6 ) && ( words[0] == "site" ) )
            found = true;
        else
        {
            text_file_writer.write_line( input_file.line( iLine ) );
            ++iLine;
        }
    }
    // When we are here, either we are at the end of the file or we are at the line "site ..."
    if ( ! found )
        throw std::runtime_error( "TLSWriter_2(): could not find \"site\"." );
    CrystalStructure crystal_structure;
    crystal_structure.set_crystal_lattice( CrystalLattice( a, b, c, Angle::from_degrees( al ), Angle::from_degrees( be ), Angle::from_degrees( ga ) ) );
    // Read the "site" lines
    found = true;
    while ( found && ( iLine != input_file.size() ) )
    {
        words = split_2( input_file.line( iLine ) );
        if ( ( words.size() == 17 ) && ( words[0] == "site" ) )
        {
            crystal_structure.add_atom( Atom( Element( words[12] ), Vector3D( TOPASstring2double( words[4] ), TOPASstring2double( words[7] ), TOPASstring2double( words[10] ) ), words[1] ) );
            ++iLine;
        }
        else
            found = false;
    }
    text_file_writer.write_line( "    ' Origin of the molecule (fractional coordinates)." );
    Vector3D c_o_m = crystal_structure.centre_of_mass();
    text_file_writer.write_line( "    prm centx " + double2string( c_o_m.x() ) );
    text_file_writer.write_line( "    prm centy " + double2string( c_o_m.y() ) );
    text_file_writer.write_line( "    prm centz " + double2string( c_o_m.z() ) );
    text_file_writer.write_line();
    text_file_writer.write_line( "    ' T11 etc. are elements of the T, L and S tensors." );
    text_file_writer.write_line( "    prm T11  0.05" );
    text_file_writer.write_line( "    prm T22  0.05" );
    text_file_writer.write_line( "    prm T33  0.05" );
    text_file_writer.write_line( "    prm T12  0.0" );
    text_file_writer.write_line( "    prm T13  0.0" );
    text_file_writer.write_line( "    prm T23  0.0" );
    text_file_writer.write_line( "    prm L11  0.0" );
    text_file_writer.write_line( "    prm L22  0.0" );
    text_file_writer.write_line( "    prm L33  0.0" );
    text_file_writer.write_line( "    prm L12  0.0" );
    text_file_writer.write_line( "    prm L13  0.0" );
    text_file_writer.write_line( "    prm L23  0.0" );
    text_file_writer.write_line( "    prm S11  0.0" );
    text_file_writer.write_line( "    prm S12  0.0" );
    text_file_writer.write_line( "    prm S13  0.0" );
    text_file_writer.write_line( "    prm S22  0.0" );
    text_file_writer.write_line( "    prm S23  0.0" );
    text_file_writer.write_line();
    text_file_writer.write_line( "    ' Choice of origin to make S symmetric (Downs, p80-81)." );
    text_file_writer.write_line( "    prm S21 = S12;" );
    text_file_writer.write_line( "    prm S31 = S13;" );
    text_file_writer.write_line( "    prm S32 = S23;" );
    text_file_writer.write_line();
    text_file_writer.write_line( "    ' Standard constraint to make equations determinate." );
    text_file_writer.write_line( "    prm S33 = -S11 - S22; : 0.0" );
    text_file_writer.write_line();
    text_file_writer.write_line( "    ' The libration is sqrt(trace(L))." );
    text_file_writer.write_line( "    prm libration = Sqrt(L11 + L22 + L33) * Rad; : 0.0" );
    text_file_writer.write_line();
    text_file_writer.write_line( "    ' Conversion of fractional coordinates to Cartesian." );
    text_file_writer.write_line( "    ' xN1   = fractional x coordinate of atom N1" );
    text_file_writer.write_line( "    ' rxN1  = Cartesian coordinate x of atom N1. 'r' always indicates a Cartesian coordinate." );
    text_file_writer.write_line( "    ' drxN1 = Correction for libration to Cartesian x coordinate of atom N1." );
    text_file_writer.write_line( "    ' crxN1 = Cartesian x coordinate of atom N1 corrected for libration:" );
    text_file_writer.write_line( "    ' crxN1 = rxN1 + drxN1" );
    text_file_writer.write_line();
    text_file_writer.write_line( "    ' Unit-cell parameters are assumed to be constant." );
    text_file_writer.write_line();
    Matrix3D f2c = crystal_structure.crystal_lattice().fractional_to_orthogonal_matrix();

    // Many matrix elements evaluate to e.g. -1.0E-18, these are set to 0.0 here.
    for ( size_t i( 0 ); i < 3; ++i )
    {
        for ( size_t j( 0 ); j < 3; ++j )
        {
            if ( nearly_zero( f2c.value( i, j ) ) )
                f2c.set_value( i, j, 0.0 );
        }
    }

    for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
    {
        std::string label = crystal_structure.atom( i ).label();
        text_file_writer.write_line( "    prm !rx" + label + " = (x" + label + " - centx)*" + double2string( f2c.value(0,0) ) + " + (y" + label + " - centy)*" + double2string( f2c.value(0,1) ) + " + (z" + label + " - centz)*" + double2string( f2c.value(0,2) ) + ";" );
        text_file_writer.write_line( "    prm !ry" + label + " = (y" + label + " - centy)*" + double2string( f2c.value(1,1) ) + " + (z" + label + " - centz)*" + double2string( f2c.value(1,2) ) + ";" );
        text_file_writer.write_line( "    prm !rz" + label + " = (z" + label + " - centz)*" + double2string( f2c.value(2,2) ) + ";" );
    }
    text_file_writer.write_line();
    text_file_writer.write_line( "    ' Calculation of Uij values in Cartesian framework. See Dunitz page 251, Eqn 5.42 and Table 5.1 page 252." );
    text_file_writer.write_line();
    for ( size_t i( 0 ); i < crystal_structure.natoms(); ++i )
    {
        std::string label = crystal_structure.atom( i ).label();
        text_file_writer.write_line( "    prm ru11" + label + " = L22*rz" + label + "^2 + L33*ry" + label + "^2 - 2*ry" + label + "*rz" + label + "*L23 - 2*ry" + label + "*S31 + 2*rz" + label + "*S21 + T11; : 0.0" );
        text_file_writer.write_line( "    prm ru22" + label + " = L11*rz" + label + "^2 + L33*rx" + label + "^2 - 2*rx" + label + "*rz" + label + "*L13 - 2*rz" + label + "*S12 + 2*rx" + label + "*S32 + T22; : 0.0" );
        text_file_writer.write_line( "    prm ru33" + label + " = L11*ry" + label + "^2 + L22*rx" + label + "^2 - 2*rx" + label + "*ry" + label + "*L12 - 2*rx" + label + "*S23 + 2*ry" + label + "*S13 + T33; : 0.0" );
        text_file_writer.write_line( "    prm ru12" + label + " = -rx" + label + "*ry" + label + "*L33 + rx" + label + "*rz" + label + "*L23 + ry" + label + "*rz" + label + "*L13 - rz" + label + "^2*L12 - rz" + label +
                                     "*S11 + rz" + label + "*S22 + rx" + label + "*S31 - ry" + label + "*S32 + T12; : 0.0" );
        text_file_writer.write_line( "    prm ru13" + label + " = -rx" + label + "*rz" + label + "*L22 + rx" + label + "*ry" + label + "*L23 - ry" + label + "^2*L13 + ry" + label + "*rz" + label + "*L12 + ry" + label +
                                     "*S11 - ry" + label + "*S33 + rz" + label + "*S23 - rx" + label + "*S21 + T13; : 0.0" );
        text_file_writer.write_line( "    prm ru23" + label + " = -ry" + label + "*rz" + label + "*L11 - rx" + label + "^2*L23 + rx" + label + "*ry" + label + "*L13 + rx" + label + "*rz" + label + "*L12 - rx" + label +
                                     "*S22 + rx" + label + "*S33 + ry" + label + "*S12 - rz" + label + "*S13 + T23; : 0.0" );
        text_file_writer.write_line();
    }
    text_file_writer.write_line( "    ' ru11N1 = u11 of atom N1 in Cartesian coordinates." );
    text_file_writer.write_line( "    ' u11N1 = u11 of atom N1 in cif coordinates." );
    text_file_writer.write_line( "    ' See R. W. Grosse-Kunstleve & P. D. Adams (2002). J. Appl. Cryst. 35, 477-480." );
    text_file_writer.write_line();
    for ( size_t i( 0 ); i < crystal_structure.natoms(); ++i )
    {
        std::string label = crystal_structure.atom( i ).label();
        std::vector<double> k;
        k.push_back( 1.0 / crystal_structure.crystal_lattice().a_star() );
        k.push_back( 1.0 / crystal_structure.crystal_lattice().b_star() );
        k.push_back( 1.0 / crystal_structure.crystal_lattice().c_star() );
        Matrix3D c2f = crystal_structure.crystal_lattice().orthogonal_to_fractional_matrix();
        size_t r;
        size_t s;
        // Many matrix elements evaluate to e.g. -1.0E-18, these are set to 0.0 here.
        for ( r = 0; r < 3; ++r )
        {
            for ( s = 0; s < 3; ++s )
            {
                if ( nearly_zero( c2f.value( r, s ) ) )
                    c2f.set_value( r, s, 0.0 );
            }
        }
        // TOPAS cannot cope with "prm n = 4 + -2;", a workaround is to write "prm n = 4 + (-2);".
        std::vector< size_t > r_s; // "r's", multiple of "r"
        std::vector< size_t > s_s; // "s's", multiple of "s"
        r_s.push_back( 1 );
        s_s.push_back( 1 );
        r_s.push_back( 2 );
        s_s.push_back( 2 );
        r_s.push_back( 3 );
        s_s.push_back( 3 );
        r_s.push_back( 1 );
        s_s.push_back( 2 );
        r_s.push_back( 1 );
        s_s.push_back( 3 );
        r_s.push_back( 2 );
        s_s.push_back( 3 );
        for ( size_t j(0); j != r_s.size(); ++j )
        {
            r = r_s[j];
            s = s_s[j];
            text_file_writer.write_line( "    prm u" + size_t2string(r) + size_t2string(s) + label + " = " + double2string( k[r-1] * k[s-1] ) + " * ( " );
            text_file_writer.write_line( "        "  + double2string( c2f.value(r-1,0) ) + "  * ( ru11" + label + " * " + double2string( c2f.value(s-1,0) ) + " + ru12" + label + " * " + double2string( c2f.value(s-1,1) ) + " + ru13" + label + " * " + double2string( c2f.value(s-1,2) ) + " ) + " );
            text_file_writer.write_line( "        (" + double2string( c2f.value(r-1,1) ) + ") * ( ru12" + label + " * " + double2string( c2f.value(s-1,0) ) + " + ru22" + label + " * " + double2string( c2f.value(s-1,1) ) + " + ru23" + label + " * " + double2string( c2f.value(s-1,2) ) + " ) + " );
            text_file_writer.write_line( "        (" + double2string( c2f.value(r-1,2) ) + ") * ( ru13" + label + " * " + double2string( c2f.value(s-1,0) ) + " + ru23" + label + " * " + double2string( c2f.value(s-1,1) ) + " + ru33" + label + " * " + double2string( c2f.value(s-1,2) ) + " ) ); : 0.0" );
        }
        text_file_writer.write_line();
    }
    for ( size_t i( 0 ); i < crystal_structure.natoms(); ++i )
    {
        std::string label = crystal_structure.atom( i ).label();
        bool keep_ADPs_positive_definite = false;
        if ( keep_ADPs_positive_definite )
            text_file_writer.write_line( "    site " + label + " x = x" + label + "; " +
                                                               " y = y" + label + "; " +
                                                               " z = z" + label + "; " +
                                         " occ " + crystal_structure.atom( i ).element().symbol() + " 1 ADPs_Keep_PD u11=u11" + label + "; u22=u22" + label + "; u33=u33" + label + "; u12=u12" + label + "; u13=u13" + label + "; u23=u23" + label + ";" );
        else
            text_file_writer.write_line( "    site " + label + " x = x" + label + "; " +
                                                               " y = y" + label + "; " +
                                                               " z = z" + label + "; " +
                                         " occ " + crystal_structure.atom( i ).element().symbol() + " 1 ADPs { =u11" + label + "; =u22" + label + "; =u33" + label + "; =u12" + label + "; =u13" + label + "; =u23" + label + "; }" );
    }
    text_file_writer.write_line();
    text_file_writer.write_line( "    ' Calculate corrected relative Cartesian coordinates. See Downs, page 83." );
    text_file_writer.write_line( "    ' Note that only intramolecular corrected distances can be trusted, the intermolecular corrected distances have no physical meaning." );
    for ( size_t i( 0 ); i < crystal_structure.natoms(); ++i )
    {
        std::string label = crystal_structure.atom( i ).label();
        text_file_writer.write_line( "    prm !drx" + label + " = 0.5*( (L22+L33)*rx" + label + "        -L12*ry" + label + "        -L13*rz" + label + " ) ; : 0.0" );
        text_file_writer.write_line( "    prm !dry" + label + " = 0.5*(      -L12*rx" + label + " + (L11+L33)*ry" + label + "        -L23*rz" + label + " ) ; : 0.0" );
        text_file_writer.write_line( "    prm !drz" + label + " = 0.5*(      -L13*rx" + label + "        -L23*ry" + label + " + (L11+L22)*rz" + label + " ) ; : 0.0" );
        text_file_writer.write_line( "    prm !crx" + label + " =                 rx" + label + " +          drx" + label + " ; : 0.0" );
        text_file_writer.write_line( "    prm !cry" + label + " =                 ry" + label + " +          dry" + label + " ; : 0.0" );
        text_file_writer.write_line( "    prm !crz" + label + " =                 rz" + label + " +          drz" + label + " ; : 0.0" );
    }
    text_file_writer.write_line();
    // Read the lines between the "site" and "restrain". This should be:

//    site H75 x ref_flag -0.04362`_0.00021 y ref_flag -0.24684`_0.00184 z ref_flag  1.35597`_0.00243 occ H 1 beq = bh;
//    prm !bond_width 0
//    prm !bond_weight 10000
//    prm !angle_width 1
//    prm !angle_weight 1
//    prm !flatten_width 0
//    prm !flatten_weight 100000
//    Distance_Restrain(   C1   C2, 1.459, 1.41782`_0.00994, bond_width, bond_weight)

    ++iLine;
    while ( ( ! found ) && ( iLine != input_file.size() ) )
    {
        size_t iPos = input_file.line( iLine ).find( "Restrain" );
        if ( iPos != std::string::npos )
            found = true;
        else
            ++iLine;
    }
    text_file_writer.write_line( "    prm !bond_width 0" );
    text_file_writer.write_line( "    prm !bond_weight 1000" );
    text_file_writer.write_line( "    prm !angle_width 0" );
    text_file_writer.write_line( "    prm !angle_weight 1" );
    text_file_writer.write_line( "    prm !flatten_width 0" );
    text_file_writer.write_line( "    prm !flatten_weight 100000" );
    // When we are here, either we are at the end of the file or we are at the line "Restrain".
    // Read the "restrain" lines
    std::vector< std::string > bond_labels_1;
    std::vector< std::string > bond_labels_2;
    std::vector< double > bond_target_values;
    std::vector< std::string > angle_labels_1;
    std::vector< std::string > angle_labels_2;
    std::vector< std::string > angle_labels_3;
    std::vector< double > angle_target_values;
    Splitter splitter( "," );
    found = true;
    while ( found && ( iLine != input_file.size() ) )
    {
        found = false;
        size_t iPos = input_file.line( iLine ).find( "Restrain" );
        if ( iPos == std::string::npos )
            break;
        else
        {
            found = true;
            iPos = input_file.line( iLine ).find( "Distance_Restrain" );
            if ( iPos == std::string::npos )
            {
                iPos = input_file.line( iLine ).find( "Angle_Restrain" );
                if ( iPos == std::string::npos )
                    throw std::runtime_error( "TLSWriter_2(): could not interpret \"Restrain\" >" + input_file.line( iLine ) + "<" );
                // Angle_Restrain( C2 C1 N3, 118.359, 118.34201`_0.68292, angle_width, angle_weight)
                std::string line = extract_delimited_text( input_file.line( iLine ), "(", ")" );
                words = splitter.split( line );
                if ( words.size() != 5 )
                    throw std::runtime_error( "TLSWriter_2(): could not interpret Angle_Restrain line. >" + input_file.line( iLine ) + "<" );
                angle_target_values.push_back( string2double( strip( words[1] ) ) );
                words = split_2( words[0] );
                if ( words.size() != 3 )
                    throw std::runtime_error( "TLSWriter_2(): could not interpret Angle_Restrain line. >" + input_file.line( iLine ) + "<" );
                angle_labels_1.push_back( words[0] );
                angle_labels_2.push_back( words[1] );
                angle_labels_3.push_back( words[2] );
            }
            else
            {
                // Distance_Restrain( C63 H75, 0.95 , 0.95015`_0.01357, bond_width, bond_weight)
                std::string line = extract_delimited_text( input_file.line( iLine ), "(", ")" );
                words = splitter.split( line );
                if ( words.size() != 5 )
                    throw std::runtime_error( "TLSWriter_2(): could not interpret Distance_Restrain line. >" + input_file.line( iLine ) + "<" );
                bond_target_values.push_back( string2double( strip( words[1] ) ) );
                words = split_2( words[0] );
                if ( words.size() != 2 )
                    throw std::runtime_error( "TLSWriter_2(): could not interpret Distance_Restrain line. >" + input_file.line( iLine ) + "<" );
                bond_labels_1.push_back( words[0] );
                bond_labels_2.push_back( words[1] );
            }
            ++iLine;
        }
    }

    // We must find the unique combinations, e.g.:
    // angle r1-r2-r3 and bond r3-r2 only require two combinations: r1-r2 and r3-r2
    // Create list of unique bonds
    std::vector< std::string > unique_labels_1;
    std::vector< std::string > unique_labels_2;
    for ( size_t i( 0 ); i != angle_labels_1.size(); ++i )
    {
        bool found( false );
        // Find bond 1
        for ( size_t j( 0 ); j != unique_labels_1.size(); ++j )
        {
            if ( ( unique_labels_1[j] == angle_labels_1[i] ) && ( unique_labels_2[j] == angle_labels_2[i] ) )
            {
                 found = true;
                 break;
            }
        }
        if ( ! found )
        {
            unique_labels_1.push_back( angle_labels_1[i] );
            unique_labels_2.push_back( angle_labels_2[i] );
        }
        found = false;
        // Find bond 2
        for ( size_t j( 0 ); j != unique_labels_1.size(); ++j )
        {
            if ( ( unique_labels_1[j] == angle_labels_3[i] ) && ( unique_labels_2[j] == angle_labels_2[i] ) )
            {
                 found = true;
                 break;
            }
        }
        if ( ! found )
        {
            unique_labels_1.push_back( angle_labels_3[i] );
            unique_labels_2.push_back( angle_labels_2[i] );
        }
    }

//    if ( write_bond_restraints )
    {
        text_file_writer.write_line( "    ' Calculate current bond lengths." );
        for ( size_t i( 0 ); i < bond_labels_1.size(); ++i )
        {
            text_file_writer.write_line( "    prm !curr_dist_" + bond_labels_1[i] + "_" + bond_labels_2[i] + " = Sqrt( (rx" + bond_labels_1[i] + "-rx" + bond_labels_2[i] + ")^2 + " +
                                                                                                                      "(ry" + bond_labels_1[i] + "-ry" + bond_labels_2[i] + ")^2 + " +
                                                                                                                      "(rz" + bond_labels_1[i] + "-rz" + bond_labels_2[i] + ")^2 ); : 0.0" );
        }
        text_file_writer.write_line();
        text_file_writer.write_line( "    ' Calculate corrected relative Cartesian coordinates. See Downs, page 83." );
        text_file_writer.write_line( "    ' Note that only intramolecular corrected distances can be trusted, the intermolecular corrected distances have no physical meaning." );
        for ( size_t i( 0 ); i < bond_labels_1.size(); ++i )
        {
            text_file_writer.write_line( "    prm !corr_dist_" + bond_labels_1[i] + "_" + bond_labels_2[i] + " = Sqrt( (crx" + bond_labels_1[i] + "-crx" + bond_labels_2[i] + ")^2 + " +
                                                                                                                      "(cry" + bond_labels_1[i] + "-cry" + bond_labels_2[i] + ")^2 + " +
                                                                                                                      "(crz" + bond_labels_1[i] + "-crz" + bond_labels_2[i] + ")^2 ); : 0.0" );
        }
        text_file_writer.write_line();
        text_file_writer.write_line( "    ' Restrain the corrected bond lengths." );
        for ( size_t i( 0 ); i < bond_labels_1.size(); ++i )
        {
            text_file_writer.write_line( "    Angle_Distance_Restrain( corr_dist_" + bond_labels_1[i] + "_" + bond_labels_2[i] + ", " + double2string( bond_target_values[i] ) + ", 0.0, bond_width, bond_weight )" );
        }
        text_file_writer.write_line();
    }
//    if ( write_angle_restraints )
    {

        text_file_writer.write_line( "    ' Calculate corrected valence angles." );
        for ( size_t i( 0 ); i < unique_labels_1.size(); ++i )
        {
            text_file_writer.write_line( "    prm !diff_x_" + unique_labels_1[i] + "_" + unique_labels_2[i] + " = crx" + unique_labels_1[i] + " - crx" + unique_labels_2[i] + "; : 0.0" );
            text_file_writer.write_line( "    prm !diff_y_" + unique_labels_1[i] + "_" + unique_labels_2[i] + " = cry" + unique_labels_1[i] + " - cry" + unique_labels_2[i] + "; : 0.0" );
            text_file_writer.write_line( "    prm !diff_z_" + unique_labels_1[i] + "_" + unique_labels_2[i] + " = crz" + unique_labels_1[i] + " - crz" + unique_labels_2[i] + "; : 0.0" );
        }
        for ( size_t i( 0 ); i < angle_labels_1.size(); ++i )
        {
            text_file_writer.write_line( "    prm !corr_ang_" + angle_labels_1[i] + "_" + angle_labels_2[i] + "_" + angle_labels_3[i] + " = Rad * ArcCos( ( diff_x_" + angle_labels_1[i] + "_" + angle_labels_2[i] + " * diff_x_" + angle_labels_3[i] + "_" + angle_labels_2[i] + " + " +
                                                                                                                                                           "diff_y_" + angle_labels_1[i] + "_" + angle_labels_2[i] + " * diff_y_" + angle_labels_3[i] + "_" + angle_labels_2[i] + " + " +
                                                                                                                                                           "diff_z_" + angle_labels_1[i] + "_" + angle_labels_2[i] + " * diff_z_" + angle_labels_3[i] + "_" + angle_labels_2[i] + " ) /" +
                                                                                                                                                    " (Sqrt(diff_x_" + angle_labels_1[i] + "_" + angle_labels_2[i] + "^2+" +
                                                                                                                                                           "diff_y_" + angle_labels_1[i] + "_" + angle_labels_2[i] + "^2+" +
                                                                                                                                                           "diff_z_" + angle_labels_1[i] + "_" + angle_labels_2[i] + "^2) * " +
                                                                                                                                                      "Sqrt(diff_x_" + angle_labels_3[i] + "_" + angle_labels_2[i] + "^2+" +
                                                                                                                                                           "diff_y_" + angle_labels_3[i] + "_" + angle_labels_2[i] + "^2+" +
                                                                                                                                                           "diff_z_" + angle_labels_3[i] + "_" + angle_labels_2[i] + "^2)) ); : 0.0" );
        }
        text_file_writer.write_line();
        text_file_writer.write_line( "    ' Restrain the corrected valence angles." );
        for ( size_t i( 0 ); i < angle_labels_1.size(); ++i )
        {
           text_file_writer.write_line( "    Angle_Distance_Restrain( corr_ang_" + angle_labels_1[i] + "_" + angle_labels_2[i] + "_" + angle_labels_3[i] + ", " + double2string( angle_target_values[i] ) + ", 0.0, angle_width, angle_weight )" );
        }
    }

    while ( iLine != input_file.size() )
    {
        // Find Out_CIF_STR_Uiso and replace by Out_CIF_ADPs_TLS
        size_t iPos = input_file.line( iLine ).find( "Out_CIF_STR_Uiso" );
        if ( iPos == std::string::npos )
            text_file_writer.write_line( input_file.line( iLine ) );
        else
            text_file_writer.write_line( input_file.line( iLine ).substr( 0, iPos ) + "Out_CIF_ADPs_TLS" + input_file.line( iLine ).substr( iPos + 16 ) );
        ++iLine;
    }
    text_file_writer.write_line( "    macro Out_CIF_ADPs_TLS( file )" );
    text_file_writer.write_line( "    {" );
    text_file_writer.write_line( "        out file" );
    text_file_writer.write_line( "        Out_String(\"data_\\n\")" );
    text_file_writer.write_line( "        Out(Get(a), \"_cell_length_a  %V\")" );
    text_file_writer.write_line( "        Out(Get(b), \"\\n_cell_length_b  %V\")" );
    text_file_writer.write_line( "        Out(Get(c), \"\\n_cell_length_c  %V\")" );
    text_file_writer.write_line( "        Out(Get(al), \"\\n_cell_angle_alpha %V\")" );
    text_file_writer.write_line( "        Out(Get(be), \"\\n_cell_angle_beta  %V\")" );
    text_file_writer.write_line( "        Out(Get(ga), \"\\n_cell_angle_gamma %V\")" );
    text_file_writer.write_line( "        Out(Get(cell_volume), \"\\n_cell_volume %V\")" );
    text_file_writer.write_line( "        Out(Get(sp_grp_char), \"\\n_symmetry_space_group_name_H-M %s\")" );
    text_file_writer.write_line( "        Out_String(\"\\nloop_\")" );
    text_file_writer.write_line( "        Out_String(\"\\n_symmetry_equiv_pos_as_xyz\")" );
    text_file_writer.write_line( "            Out(Get(sp_xyzs_txt),  \"%s\")" );
    text_file_writer.write_line( "        Out_String(\"\\nloop_\")" );
    text_file_writer.write_line( "        Out_String(\"\\n_atom_site_label\")" );
    text_file_writer.write_line( "        Out_String(\"\\n_atom_site_type_symbol\")" );
    text_file_writer.write_line( "        Out_String(\"\\n_atom_site_symmetry_multiplicity\")" );
    text_file_writer.write_line( "        Out_String(\"\\n_atom_site_fract_x\")" );
    text_file_writer.write_line( "        Out_String(\"\\n_atom_site_fract_y\")" );
    text_file_writer.write_line( "        Out_String(\"\\n_atom_site_fract_z\")" );
    text_file_writer.write_line( "        Out_String(\"\\n_atom_site_occupancy\")" );
    text_file_writer.write_line( "        atom_out file append" );
    text_file_writer.write_line( "            load out_record out_fmt out_eqn" );
    text_file_writer.write_line( "            {" );
    text_file_writer.write_line( "                \"\\n%s\" = Get_From_String(Get(current_atom), site);" );
    text_file_writer.write_line( "                \" %s\" = Get_From_String(Get(current_atom), atom);" );
    text_file_writer.write_line( "                \" %3.0f\" = Get_From_String(Get(current_atom), num_posns);" );
    text_file_writer.write_line( "                \" %V\" = Get_From_String(Get(current_atom), x);" );
    text_file_writer.write_line( "                \" %V\" = Get_From_String(Get(current_atom), y);" );
    text_file_writer.write_line( "                \" %V\" = Get_From_String(Get(current_atom), z);" );
    text_file_writer.write_line( "                \" %V\" = Get_From_String(Get(current_atom), occ);" );
    text_file_writer.write_line( "            }" );
    text_file_writer.write_line( "        out file append" );
    text_file_writer.write_line( "        Out_String(\"\\nloop_\")" );
    text_file_writer.write_line( "        Out_String(\"\\n_atom_site_aniso_label\")" );
    text_file_writer.write_line( "        Out_String(\"\\n_atom_site_aniso_U_11\")" );
    text_file_writer.write_line( "        Out_String(\"\\n_atom_site_aniso_U_22\")" );
    text_file_writer.write_line( "        Out_String(\"\\n_atom_site_aniso_U_33\")" );
    text_file_writer.write_line( "        Out_String(\"\\n_atom_site_aniso_U_12\")" );
    text_file_writer.write_line( "        Out_String(\"\\n_atom_site_aniso_U_13\")" );
    text_file_writer.write_line( "        Out_String(\"\\n_atom_site_aniso_U_23\")" );
    text_file_writer.write_line( "        atom_out file append adps" );
    text_file_writer.write_line( "            load out_record out_fmt out_eqn" );
    text_file_writer.write_line( "            {" );
    text_file_writer.write_line( "                \"\\n%s\" = Get_From_String(Get(current_atom), site);" );
    text_file_writer.write_line( "                \" %V\" = Get_From_String(Get(current_atom), u11);" );
    text_file_writer.write_line( "                \" %V\" = Get_From_String(Get(current_atom), u22);" );
    text_file_writer.write_line( "                \" %V\" = Get_From_String(Get(current_atom), u33);" );
    text_file_writer.write_line( "                \" %V\" = Get_From_String(Get(current_atom), u12);" );
    text_file_writer.write_line( "                \" %V\" = Get_From_String(Get(current_atom), u13);" );
    text_file_writer.write_line( "                \" %V\" = Get_From_String(Get(current_atom), u23);" );
    text_file_writer.write_line( "            }" );
    text_file_writer.write_line( "        out file append" );
    text_file_writer.write_line( "        Out_String(\"\\n\")" );
    text_file_writer.write_line( "    }" );
    copy_text_file( append_to_file_name( input_file_name, "_TLS" ), replace_extension( append_to_file_name( input_file_name, "_TLS" ), "org" ) );
}

// ********************************************************************************
