#ifndef TLS_H
#define TLS_H

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

class Angle;
class AnisotropicDisplacementParameters;

#include "Matrix3D.h"
#include "SymmetricMatrix3D.h"
#include "Vector3D.h"

#include <string>
#include <vector>
/*
    try // Calculate libration eigenvectors
    {
        if ( argc != 2 )
            throw std::runtime_error( "Please give the name of a .inp file that needs to be converted to _EV.txt." );
        FileName input_file_name( argv[ 1 ] );
        TextFileReader_2 input_file( input_file_name );
        TextFileWriter text_file_writer( append_to_file_name( input_file_name, "_EV" ) );
 //   ' Origin of the molecule.
//    prm !centx 0.0
//    prm !centy 0.0
//    prm !centz 0.0
        double centx = read_keyword( "centx", input_file );
        double centy = read_keyword( "centy", input_file );
        double centz = read_keyword( "centz", input_file );
        Vector3D origin( centx, centy, centz );

//	macro VV { val_on_continue = Val * Rand(0.5, 2); }
//    ' T11 etc are elements of T, L and S tensors.
//    prm !T11 0.05    VV
//    prm !T12 0.0    VV
//    prm !T13 0.0    VV
//    prm !T22 0.05    VV
//    prm !T23 0.0    VV
//    prm !T33 0.05    VV
//    prm !L11 0.0    VV
//    prm !L12 0.0    VV
//    prm !L13 0.0    VV
//    prm !L22 0.02    VV
//    prm !L23 0.0    VV
//    prm !L33 0.02    VV
//    prm !S11  0.0
//    prm !S12  0.0
//    prm !S13  0.0
//    prm !S22  0.0
//    prm !S23  0.0

        double T11 = read_keyword( "T11", input_file );
        double T22 = read_keyword( "T22", input_file );
        double T33 = read_keyword( "T33", input_file );
        double T12 = read_keyword( "T12", input_file );
        double T13 = read_keyword( "T13", input_file );
        double T23 = read_keyword( "T23", input_file );
        double L11 = read_keyword( "L11", input_file );
        double L22 = read_keyword( "L22", input_file );
        double L33 = read_keyword( "L33", input_file );
        double L12 = read_keyword( "L12", input_file );
        double L13 = read_keyword( "L13", input_file );
        double L23 = read_keyword( "L23", input_file );

        SymmetricMatrix3D L = SymmetricMatrix3D( L11, L22, L33, L12, L13, L23 );
        std::vector< double > eigenvalues;
        std::vector< NormalisedVector3D > eigenvectors;
        calculate_eigenvalues( L, eigenvalues, eigenvectors );
        std::cout << "Eigenvalues : " + double2string( eigenvalues[ 0 ] ) + ", " + double2string( eigenvalues[ 1 ] ) + ", " + double2string( eigenvalues[ 2 ] ) << std::endl;
        text_file_writer.write_line( "Eigenvalues : " + double2string( eigenvalues[ 0 ] ) + ", " + double2string( eigenvalues[ 1 ] ) + ", " + double2string( eigenvalues[ 2 ] ) );
        std::cout << "Eigenvectors:" << std::endl;
        for ( size_t i( 0 ); i != eigenvectors.size(); ++i )
            eigenvectors[i].show();
        text_file_writer.write_line( "Eigenvectors:" );
        for ( size_t i( 0 ); i != eigenvectors.size(); ++i )
            text_file_writer.write_line( eigenvectors[i].to_string() );
        std::cout << "Origin = " << std::endl;
        origin.show();
        text_file_writer.write_line( "Origin = " + origin.to_string() );
        // Convert from Cartesian to fractional coordinates
        CrystalLattice crystal_lattice = read_lattice_parameters( input_file );
        // Scale eigenvalues such that the largest = 1
        eigenvalues[0] /= eigenvalues[2];
        eigenvalues[1] /= eigenvalues[2];
        eigenvalues[2] /= eigenvalues[2];
        std::vector< Vector3D > ev_fractional;
        for ( size_t i( 0 ); i != eigenvectors.size(); ++i )
        {
            ev_fractional.push_back( origin + crystal_lattice.orthogonal_to_fractional_matrix() * ( eigenvalues[i] * eigenvectors[i] ) );
        }
        text_file_writer.write_line( "O Un 4 " + double2string( origin.value( 0 ) ) + " " + double2string( origin.value( 1 ) ) + " " + double2string( origin.value( 2 ) ) + " 0 0.3" );
        for ( size_t i( 0 ); i != eigenvectors.size(); ++i )
        {
            text_file_writer.write_line( "EV" + size_t2string( i ) + " Un 4 " + double2string( ev_fractional[i].value( 0 ) ) + " " +
                                                                                double2string( ev_fractional[i].value( 1 ) ) + " " +
                                                                                double2string( ev_fractional[i].value( 2 ) ) +  " 0 0.3" );
        }
    MACRO_END_GAME
*/


/*

  Units: Angstrom and radians.

  It is up to the user to impose S33 = -S11 - S22.
  
  It is up to the user to calculate the centre of reaction, then to move the origin and then to impose S21 = S12, S31 = S13, S32 = S23
  to make S symmetric.
  
  If on inversion:
  
  origin is [ 0.0, 0.0, 0.0 ]
  
  All elements of S are 0.0.

*/
class TLS
{
public:

    // Default constructor
    TLS();

    TLS( const SymmetricMatrix3D & T, const SymmetricMatrix3D & L, const Matrix3D & S );

    SymmetricMatrix3D T() const { return T_; }
    SymmetricMatrix3D L() const { return L_; }
    Matrix3D S() const { return S_; }

    bool is_on_inversion() const { return is_on_inversion_; }
    
    // Throws if is_on_inversion is set to true and S is not the zero matrix.
    void set_is_on_inversion( const bool value );

    AnisotropicDisplacementParameters U( const Vector3D & r ) const;

    // Origin choice that makes S symmetric
    Vector3D centre_of_reaction() const;

    // sqrt( trace( L ) )
    Angle libration() const;

    std::vector< std::string > TOPAS_lines( const Vector3D & r ) const;

    // Note that only intramolecular corrected distances can be trusted, the intermolecular corrected distances have no physical meaning.
    Vector3D corrected_relative_Cartesian_coordinate( const Vector3D & r ) const;
    // To enforce that this function can only be used to correct distances within a rigid body,
    // the current interface only allows you to calculate a corrected distance between two points.
    // Of course, this is very inefficient, but because we generally only have 60 or so atoms, this should not matter a lot.
    // Other possible implementation: convert a Vector3D to a corrected Vector3D or return the Matrix3D that corrects a Vector3D
//    double corrected_distance( const Vector3D & lhs, const Vector3D & rhs ) const;
    
private:

    Vector3D origin_;
    SymmetricMatrix3D T_;
    SymmetricMatrix3D L_;
    Matrix3D S_;
    SymmetricMatrix3D correction_matrix_;
    bool is_on_inversion_;
};

#endif // TLS_H

