/* *********************************************
Copyright (c) 2013-2020, Cornelis Jan (Jacco) van de Streek
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

#include "AnalyseRings.h"
//#include "CollectionOfPoints.h"
#include "Angle.h"
#include "CyclicInteger.h"
#include "MathFunctions.h"
#include "Plane.h"
#include "Vector3D.h"
#include "Sort.h"
#include "Vector3DCalculations.h"
#include "3DCalculations.h"

#include <stdexcept>
#include <iostream> // For testing only
#include <cmath>

// ********************************************************************************

FiveMemberedRingAnalyser::FiveMemberedRingAnalyser():
is_planar_(false),
is_envelope_(false),
is_double_envelope_(false)
{}

// ********************************************************************************

FiveMemberedRingAnalyser::FiveMemberedRingAnalyser( const std::vector< Vector3D > & points )
{
    analyse( points );
}

// ********************************************************************************

void FiveMemberedRingAnalyser::analyse( const std::vector< Vector3D > & points )
{
    if ( points.size() != 5 )
        throw std::runtime_error( "FiveMemberedRingAnalyser(): number of atoms must be five." );
    is_planar_ = false;
    is_envelope_ = false;
    is_double_envelope_ = false;
    points_ =  points;
    Plane plane( points );
    // Completely flat = 0.003
    // Normal chair: 0.1
    // Borderline cases: 0.045
    if ( root_mean_square_devation_from_mean_plane( points, plane ) < 0.025 )
    {
        is_planar_ = true;
        return;
    }

    distances_from_plane_.clear();
    rmsds_from_mean_plane_.clear();
    
    // Taking each point in turn as the speculative odd one out, calculate their distance from the plane formed by the other points.
    for ( size_t i( 0 ); i != 5; ++i )
    {
        std::vector< Vector3D > points2;
        CyclicInteger ci( 0, 4, i+1 );
        points2.push_back( points[ci.next_value()] );
        points2.push_back( points[ci.next_value()] );
        points2.push_back( points[ci.next_value()] );
        points2.push_back( points[ci.next_value()] );
        Plane plane_1234( points2 );
        rmsds_from_mean_plane_.push_back( root_mean_square_devation_from_mean_plane( points2, plane_1234 ) );
        distances_from_plane_.push_back( plane_1234.distance( points[i] ) );
    }
    sorted_map_ = sort( distances_from_plane_ );

    // For the two greatest distances, determine the two signed distances to the plane through the remaining three points
    // But only of the RMSD from the plane is very low (< 0.075).
    
    std::vector< Vector3D > points2;
    for ( size_t i( 0 ); i != 5; ++i )
    {
        if ( ( i != sorted_map_[4] ) && ( i != sorted_map_[3] ) )
            points2.push_back( points[i] );
    }
    Plane plane_123( points2 );
    double dist1 = plane_123.signed_distance( points[ sorted_map_[4] ] );
    double dist2 = plane_123.signed_distance( points[ sorted_map_[3] ] );
    if ( ( sign( dist1 ) * sign( dist2 ) ) < 0 )
    {
        // The points lie on different sides of the plane.
        // The following should probably be done based on dist1 and dist2, but I do not have the proper limits
        // from tests I get the impression that the old criterion is much more specific
        if ( std::abs( distances_from_plane_[ sorted_map_[4] ] - distances_from_plane_[ sorted_map_[3] ] ) < 0.0275 )
        {
            if ( rmsds_from_mean_plane_[ sorted_map_[3] ] < 0.075 )
            {
                is_double_envelope_ = true;
                root_mean_square_devation_from_mean_plane_ = rmsds_from_mean_plane_[ sorted_map_[4] ];
                distance_ = distances_from_plane_[ sorted_map_[4] ];
                unique_envelope_point_1_ = sorted_map_[4];
                unique_envelope_point_2_ = sorted_map_[3];
                return;
            }
        }
    }
    // When we are here, the ring is not flat and it is not a double envelope.
    is_envelope_ = true;
    unique_envelope_point_1_ = sorted_map_[4];
    root_mean_square_devation_from_mean_plane_ = rmsds_from_mean_plane_[ sorted_map_[4] ];
    distance_ = distances_from_plane_[ sorted_map_[4] ];
    return;
    
    root_mean_square_devation_from_mean_plane_ = 1.0;
    distance_ = 0.0;
    for ( size_t i( 0 ); i != 5; ++i )
    {
        std::vector< Vector3D > points2;
        CyclicInteger ci( 0, 4, i+1 );
        points2.push_back( points[ci.next_value()] );
        points2.push_back( points[ci.next_value()] );
        points2.push_back( points[ci.next_value()] );
        points2.push_back( points[ci.next_value()] );
        Plane plane_1234( points2 );
        double root_mean_square_devation_from_mean_plane = ::root_mean_square_devation_from_mean_plane( points2, plane_1234 );
        double distance = plane_1234.distance( points[i] );
        if ( root_mean_square_devation_from_mean_plane > 0.06 )
            continue;
        if ( std::abs( root_mean_square_devation_from_mean_plane - root_mean_square_devation_from_mean_plane_ ) < 0.03 )
        {
            if ( is_double_envelope_ )
                throw std::runtime_error( "FiveMemberedRingAnalyser(): triple envelope found." );
            is_double_envelope_ = true;
            if ( ! is_envelope_ )
                throw std::runtime_error( "FiveMemberedRingAnalyser(): second envelope found before first." );
            is_envelope_ = false;
            if ( root_mean_square_devation_from_mean_plane < root_mean_square_devation_from_mean_plane_ )
            {
                root_mean_square_devation_from_mean_plane_ = root_mean_square_devation_from_mean_plane;
                distance_ = distance;
                unique_envelope_point_2_ = unique_envelope_point_1_;
                unique_envelope_point_1_ = i;
            }
            else
            {
                unique_envelope_point_2_ = i;
            }
        }
        else if ( root_mean_square_devation_from_mean_plane < root_mean_square_devation_from_mean_plane_ )
        {
            root_mean_square_devation_from_mean_plane_ = root_mean_square_devation_from_mean_plane;
            distance_ = distance;
            is_envelope_ = true;
            unique_envelope_point_1_ = i;
        }
    }
}

// ********************************************************************************

FiveMemberedRingAnalyser::GeometryType FiveMemberedRingAnalyser::axial_or_equatorial( const Vector3D & point ) const
{
    if ( ! is_envelope() )
        return NONE;
    // Find the point among the five points that is closest to the new point
    size_t smallest_distance_index = 0;
    double smallest_distance = ( points_[0] - point ).length();
    for ( size_t i( 1 ); i != 5; ++i )
    {
        double distance = ( points_[i] - point ).length();
        if ( distance < smallest_distance )
        {
            smallest_distance = distance;
            smallest_distance_index = i;
        }
    }
    if ( smallest_distance < 0.1 )
        throw std::runtime_error( "FiveMemberedRingAnalyser::axial_or_equatorial(): points too close together." );
    if ( smallest_distance > 3.0 )
        throw std::runtime_error( "FiveMemberedRingAnalyser::axial_or_equatorial(): atoms not bound." );
    if ( smallest_distance_index != unique_envelope_point_1_ )
    {
        std::cout << "Warning: FiveMemberedRingAnalyser::axial_or_equatorial(): nearest ring atom is not the unique envelope point." << std::endl;
        return NONE;
    }
    NormalisedVector3D r = normalised_vector( point - points_[smallest_distance_index] );
    std::vector< Vector3D > points2;
    CyclicInteger ci( 0, 4, smallest_distance_index+1 );
    points2.push_back( points_[ci.next_value()] );
    points2.push_back( points_[ci.next_value()] );
    ++ci;
    points2.push_back( points_[ci.next_value()] );
    points2.push_back( points_[ci.next_value()] );
    Plane plane_1245( points2 );
    if ( ( plane_1245.normal()*r ) < 0.0 )
        r = -r;
    if ( angle( plane_1245.normal(), r ) < Angle::from_degrees( 30.0 ) )
        return AXIAL;
    return EQUATORIAL;

    if ( angle( plane_1245.normal(), r ) > Angle::from_degrees( 60.0 ) )
        return EQUATORIAL;
    throw std::runtime_error( "FiveMemberedRingAnalyser::axial_or_equatorial(): geometry is envelope, but neighbour is neither axial nor equatorial." );
}

// ********************************************************************************

SixMemberedRingAnalyser::SixMemberedRingAnalyser():
is_planar_(false),
is_chair_(false),
is_boat_(false),
is_twisted_boat_(false)
{}

// ********************************************************************************

SixMemberedRingAnalyser::SixMemberedRingAnalyser( const std::vector< Vector3D > & points )
{
    analyse( points );
}

// ********************************************************************************

// The actual tolerances probably depend on the elements, so a piperazine ring will require different tolerances than a hexane ring
void SixMemberedRingAnalyser::analyse( const std::vector< Vector3D > & points )
{
    if ( points.size() != 6 )
        throw std::runtime_error( "SixMemberedRingAnalyser(): number of atoms must be six." );
    is_planar_ = false;
    is_chair_ = false;
    is_boat_ = false;
    is_twisted_boat_ = false;
    points_ = points;
    Plane plane( points );
    if ( root_mean_square_devation_from_mean_plane( points, plane ) < 0.1 )
    {
        is_planar_ = true;
    }
    
    double in_plane_tolerance( 1.5 );
    double out_of_plane_tolerance( 0.25 );
    CyclicInteger ci( 0, 5, 0 );
    if ( ( std::abs( plane.signed_distance( points[0] ) ) < in_plane_tolerance ) && ( std::abs( plane.signed_distance( points[ci.plus_n(3)] ) ) < in_plane_tolerance ) &&
         ( std::abs( plane.signed_distance( points[ci.plus_n(1)] ) ) > out_of_plane_tolerance ) && ( std::abs( plane.signed_distance( points[ci.plus_n(2)] ) ) > out_of_plane_tolerance ) &&
         ( std::abs( plane.signed_distance( points[ci.plus_n(4)] ) ) > out_of_plane_tolerance ) && ( std::abs( plane.signed_distance( points[ci.plus_n(5)] ) ) > out_of_plane_tolerance ) &&
         ( sign( plane.signed_distance( points[ci.plus_n(1)] ) ) == -sign( plane.signed_distance( points[ci.plus_n(2)] ) ) ) &&
         ( sign( plane.signed_distance( points[ci.plus_n(4)] ) ) == -sign( plane.signed_distance( points[ci.plus_n(5)] ) ) ) &&
         ( sign( plane.signed_distance( points[ci.plus_n(1)] ) ) == -sign( plane.signed_distance( points[ci.plus_n(5)] ) ) ) )
    {
        is_twisted_boat_ = true;
    }
    ++ci;
    if ( ( std::abs( plane.signed_distance( points[1] ) ) < in_plane_tolerance ) && ( std::abs( plane.signed_distance( points[ ci.plus_n(3)] ) ) < in_plane_tolerance ) &&
         ( std::abs( plane.signed_distance( points[ci.plus_n(1)] ) ) > out_of_plane_tolerance ) && ( std::abs( plane.signed_distance( points[ci.plus_n(2)] ) ) > out_of_plane_tolerance ) &&
         ( std::abs( plane.signed_distance( points[ci.plus_n(4)] ) ) > out_of_plane_tolerance ) && ( std::abs( plane.signed_distance( points[ci.plus_n(5)] ) ) > out_of_plane_tolerance ) &&
         ( sign( plane.signed_distance( points[ci.plus_n(1)] ) ) == -sign( plane.signed_distance( points[ci.plus_n(2)] ) ) ) &&
         ( sign( plane.signed_distance( points[ci.plus_n(4)] ) ) == -sign( plane.signed_distance( points[ci.plus_n(5)] ) ) ) &&
         ( sign( plane.signed_distance( points[ci.plus_n(1)] ) ) == -sign( plane.signed_distance( points[ci.plus_n(5)] ) ) ) )
    {
        is_twisted_boat_ = true;
    }
    ++ci;
    if ( ( std::abs( plane.signed_distance( points[2] ) ) < in_plane_tolerance ) && ( std::abs( plane.signed_distance( points[ci.plus_n(3)] ) ) < in_plane_tolerance ) &&
         ( std::abs( plane.signed_distance( points[ci.plus_n(1)] ) ) > out_of_plane_tolerance ) && ( std::abs( plane.signed_distance( points[ci.plus_n(2)] ) ) > out_of_plane_tolerance ) &&
         ( std::abs( plane.signed_distance( points[ci.plus_n(4)] ) ) > out_of_plane_tolerance ) && ( std::abs( plane.signed_distance( points[ci.plus_n(5)] ) ) > out_of_plane_tolerance ) &&
         ( sign( plane.signed_distance( points[ci.plus_n(1)] ) ) == -sign( plane.signed_distance( points[ci.plus_n(2)] ) ) ) &&
         ( sign( plane.signed_distance( points[ci.plus_n(4)] ) ) == -sign( plane.signed_distance( points[ci.plus_n(5)] ) ) ) &&
         ( sign( plane.signed_distance( points[ci.plus_n(1)] ) ) == -sign( plane.signed_distance( points[ci.plus_n(5)] ) ) ) )
    {
        is_twisted_boat_ = true;
    }
    
    for ( size_t i( 0 ); i != 3; ++i )
    {
        std::vector< Vector3D > points2;
        CyclicInteger ci( 0, 5, i+1 );
        points2.push_back( points[ci.next_value()] );
        points2.push_back( points[ci.next_value()] );
        ++ci;
        points2.push_back( points[ci.next_value()] );
        points2.push_back( points[ci.next_value()] );
        Plane plane_1245( points2 );
        if ( root_mean_square_devation_from_mean_plane( points2, plane_1245 ) > 0.075 )
            continue;
        if ( std::abs( plane_1245.signed_distance( points[i] ) ) < 0.1 )
            continue;
        if ( std::abs( plane_1245.signed_distance( points[i+3] ) ) < 0.1 )
            continue;
        if ( sign( plane_1245.signed_distance( points[i] ) ) == -sign( plane_1245.signed_distance( points[i+3] ) ) )
            is_chair_ = true;
        else
            is_boat_ = true;
    }

    size_t nhits( 0 );
    if ( is_planar_ )
        ++nhits;
    if ( is_twisted_boat_ )
        ++nhits;
    if ( is_chair_ )
        ++nhits;
    if ( is_boat_ )
        ++nhits;
    if ( nhits == 0 )
        std::cout << "SixMemberedRingAnalyser::analyse(): ERROR: conformation type could not be assigned." << std::endl;
    if ( nhits > 1 )
        std::cout << "SixMemberedRingAnalyser::analyse(): ERROR: more than one type of conformation matched." << std::endl;

}

// ********************************************************************************

SixMemberedRingAnalyser::GeometryType SixMemberedRingAnalyser::axial_or_equatorial( const Vector3D & point ) const
{
    if ( ! is_chair() )
        return NONE;
    // Find the point among the six points that is closest to the new point
    size_t smallest_distance_index = 0;
    double smallest_distance = ( points_[0] - point ).length();
    for ( size_t i( 1 ); i != 6; ++i )
    {
        double distance = ( points_[i] - point ).length();
        if ( distance < smallest_distance )
        {
            smallest_distance = distance;
            smallest_distance_index = i;
        }
    }
    if ( smallest_distance < 0.1 )
        throw std::runtime_error( "SixMemberedRingAnalyser::axial_or_equatorial(): points too close together." );
    if ( smallest_distance > 3.0 )
        throw std::runtime_error( "SixMemberedRingAnalyser::axial_or_equatorial(): atoms not bound." );
    NormalisedVector3D r = normalised_vector( point - points_[smallest_distance_index] );
    std::vector< Vector3D > points2;
    CyclicInteger ci( 0, 5, smallest_distance_index+1 );
    points2.push_back( points_[ci.next_value()] );
    points2.push_back( points_[ci.next_value()] );
    ++ci;
    points2.push_back( points_[ci.next_value()] );
    points2.push_back( points_[ci.next_value()] );
    Plane plane_1245( points2 );
    if ( ( plane_1245.normal()*r ) < 0.0 )
        r = -r;
    if ( angle( plane_1245.normal(), r ) < Angle::from_degrees( 30.0 ) )
        return AXIAL;
    return EQUATORIAL;

    if ( angle( plane_1245.normal(), r ) > Angle::from_degrees( 60.0 ) )
        return EQUATORIAL;
    throw std::runtime_error( "SixMemberedRingAnalyser::axial_or_equatorial(): geometry is chair, but neighbour is neither axial nor equatorial." );
}

// ********************************************************************************

