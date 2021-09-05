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

#include "MC_alkanes.h"
#include "3DCalculations.h"
#include "Angle.h"
#include "FileName.h"
#include "NormalisedVector3D.h"
#include "TextFileWriter.h"
#include "Utilities.h"
#include "Vector3D.h"

// ********************************************************************************

std::vector< Vector3D > build_alkane( const size_t n, const std::vector< Angle > & torsion_angles )
{
//    size_t n; // Number of carbon atoms
    Angle C_C_C = Angle::from_degrees( 113.5 );
    std::vector< Vector3D > coordinates;
    coordinates.reserve( n );
    coordinates.push_back( Vector3D() );
    coordinates.push_back( Vector3D( 1.54, 0.0, 0.0 ) );
    coordinates.push_back( Vector3D( 1.54 + ( C_C_C - Angle::from_degrees( 90.0 ) ).sine() * 1.54, ( C_C_C - Angle::from_degrees( 90.0 ) ).cosine() * 1.54, 0.0 ) );
    for ( size_t i( 3 ); i != n; ++i )
    {
        Vector3D new_coordinate = coordinates[i-1] + 1.54 * normalised_vector( coordinates[i-2] - coordinates[i-3] );
        NormalisedVector3D axis = normalised_vector( coordinates[i-2] - coordinates[i-1] );
        new_coordinate = rotate_point_about_axis( new_coordinate, coordinates[i-1], axis, torsion_angles[i] );
        coordinates.push_back( new_coordinate );
    }
    return coordinates;
}

// ********************************************************************************

void save_as_xyz( const std::vector< Vector3D > & coordinates, const FileName & file_name )
{
    TextFileWriter text_file_writer( file_name );
    text_file_writer.write_line( size_t2string( coordinates.size() ) );
    text_file_writer.write_line( "Comment" );
    for ( size_t i( 0 ); i != coordinates.size(); ++i )
    {
        text_file_writer.write_line( "C " + double2string( coordinates[i].x() ) + " " + double2string( coordinates[i].y() ) + " " + double2string( coordinates[i].z() ) );
    }
}

// ********************************************************************************

bool there_is_overlap( const std::vector< Vector3D > & coordinates, const double exclusion_distance )
{
    size_t n = coordinates.size();
    double exclusion_distance2 = square( exclusion_distance );
    // Find all interatomc distances, except those between nearest neightbours (because they are bonded)
    for ( size_t i( 0 ); i != n - 3; ++i )
    {
        for ( size_t j( i + 3 ); j != n; ++j )
        {
            double distance2 = ( coordinates[i] - coordinates[j] ).norm2();
            if ( distance2 < exclusion_distance2 )
                return true;
        }
    }
    return false;
}

// ********************************************************************************

