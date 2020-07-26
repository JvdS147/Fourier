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

#include "VoidsFinder.h"
#include "3DCalculations.h"
#include "CrystalStructure.h"
#include "MathFunctions.h"
#include "Plane.h"
#include "RandomNumberGenerator.h"

#include <iostream>
#include <stdexcept>

namespace {

// Position is in Fractional coordinates
bool intersects_atoms( const CrystalStructure & crystal_structure, const Vector3D & position, const std::vector< double > & distances2 )
{
    for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
    {
        if ( crystal_structure.crystal_lattice().shortest_distance2( crystal_structure.atom( i ).position(), position ) < distances2[i] )
            return true;
    }
    return false;
}

double calculate_volume( const CrystalStructure & crystal_structure, const std::vector< Vector3D > & void_spheres, const double probe_radius )
{
    std::cout << "Now calculating the volume." << std::endl;
    if ( void_spheres.size() == 0 )
        return 0.0;
    RandomNumberGenerator_double random_number_generator;
    const size_t nprobes( ( crystal_structure.crystal_lattice().volume() * 1000.0 ) / crystal_structure.space_group().nsymmetry_operators() ); // 1000 sampling points per A^3
    size_t ninside_voids( 0 );
    double probe_radius2 = square( probe_radius );
    for ( size_t ix( 0 ); ix != nprobes; ++ix )
    {
        Vector3D random_vector( random_number_generator.next_number(), random_number_generator.next_number(), random_number_generator.next_number() );
        random_vector = crystal_structure.crystal_lattice().fractional_to_orthogonal( random_vector );
        for ( size_t i( 0 ); i != void_spheres.size(); ++i )
        {
            Vector3D difference_vector = random_vector - void_spheres[i];
            if ( difference_vector.norm2() < probe_radius2 )
            {
                ++ninside_voids;
                break;
            }
        }
    }
    return ( static_cast<double>( ninside_voids ) / static_cast<double>( nprobes ) ) * crystal_structure.crystal_lattice().volume();
}

} // namespace

// ********************************************************************************

// @@ Must use symmetry and we must somehow order the spheres to speed things up.

double find_voids( const CrystalStructure & crystal_structure, const double probe_radius )
{
    if ( crystal_structure.natoms() == 0 )
        return crystal_structure.crystal_lattice().volume();
    if ( ! crystal_structure.space_group_symmetry_has_been_applied() )
        throw std::runtime_error( "find_voids(): space-group symmetry has not been applied for input crystal structure." );
    double grid_spacing( 0.15 );
    double probe_reduction = probe_radius * 0.05;
    double probe_radius_intermediate = probe_radius - probe_reduction;
    std::vector< Vector3D > void_spheres_intermediate;
    std::vector< double > distances2;
    // For each atom, precalculate square( VdW radius + probe_radius )
    for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
        distances2.push_back( square( crystal_structure.atom( i ).element().Van_der_Waals_radius() + probe_radius_intermediate ) );
    // Our CrystalLattice class uses the a-along-x convention

    // We need the enclosing box / bounding box
    Vector3D min_min_min;
    Vector3D max_max_max;
    crystal_structure.crystal_lattice().enclosing_box( min_min_min, max_max_max );

    for ( int ix = round_to_int( min_min_min.x() / grid_spacing ) - 1; ix != round_to_int( max_max_max.x() / grid_spacing ) + 1; ++ix )
    {
        for ( int iy = round_to_int( min_min_min.y() / grid_spacing ) - 1; iy != round_to_int( max_max_max.y() / grid_spacing ) + 1; ++iy )
        {
            for ( int iz = round_to_int( min_min_min.z() / grid_spacing ) - 1; iz != round_to_int( max_max_max.z() / grid_spacing ) + 1; ++iz )
            {
                Vector3D position( ix * grid_spacing, iy * grid_spacing, iz * grid_spacing ); // Cartesian coordinates
                // Check that the probe point lies within the unit cell
                Vector3D fractional_position = crystal_structure.crystal_lattice().orthogonal_to_fractional( position );
                if ( ( fractional_position.x() < 0.0 ) || ( fractional_position.x() > 1.0 ) ||
                     ( fractional_position.y() < 0.0 ) || ( fractional_position.y() > 1.0 ) ||
                     ( fractional_position.z() < 0.0 ) || ( fractional_position.z() > 1.0 ) )
                    continue;
                if ( ! intersects_atoms( crystal_structure, fractional_position, distances2 ) )
                    void_spheres_intermediate.push_back( position ); // Cartesian coordinates
            }
        }
    }
    std::cout << "Number of intermediate spheres = " << void_spheres_intermediate.size() << std::endl;
    // Above was with smaller radius (-0.05). Now for every point try again, this time with the real radius, but on nine different points (vertices + centre of a cube).
    std::vector< Vector3D > void_spheres_final;
    // For each atom, precalculate square( VdW radius + probe_radius )
    distances2.clear();
    for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
        distances2.push_back( square( crystal_structure.atom( i ).element().Van_der_Waals_radius() + probe_radius ) );
    double delta = probe_reduction / sqrt(3.0); // 0.028868 Cartesian coordinates
    for ( size_t i( 0 ); i != void_spheres_intermediate.size(); ++i )
    {
        Vector3D position( void_spheres_intermediate[i] ); // Cartesian coordinates
        if ( ! intersects_atoms( crystal_structure, crystal_structure.crystal_lattice().orthogonal_to_fractional( position ), distances2 ) )
            void_spheres_final.push_back( position ); // Cartesian coordinates
        position = Vector3D( void_spheres_intermediate[i] + Vector3D( -delta, -delta, -delta ) ); // Cartesian coordinates
        if ( ! intersects_atoms( crystal_structure, crystal_structure.crystal_lattice().orthogonal_to_fractional( position ), distances2 ) )
            void_spheres_final.push_back( position ); // Cartesian coordinates
        position = Vector3D( void_spheres_intermediate[i] + Vector3D( -delta, -delta, +delta ) ); // Cartesian coordinates
        if ( ! intersects_atoms( crystal_structure, crystal_structure.crystal_lattice().orthogonal_to_fractional( position ), distances2 ) )
            void_spheres_final.push_back( position ); // Cartesian coordinates
        position = Vector3D( void_spheres_intermediate[i] + Vector3D( -delta, +delta, -delta ) ); // Cartesian coordinates
        if ( ! intersects_atoms( crystal_structure, crystal_structure.crystal_lattice().orthogonal_to_fractional( position ), distances2 ) )
            void_spheres_final.push_back( position ); // Cartesian coordinates
        position = Vector3D( void_spheres_intermediate[i] + Vector3D( -delta, +delta, +delta ) ); // Cartesian coordinates
        if ( ! intersects_atoms( crystal_structure, crystal_structure.crystal_lattice().orthogonal_to_fractional( position ), distances2 ) )
            void_spheres_final.push_back( position ); // Cartesian coordinates
        position = Vector3D( void_spheres_intermediate[i] + Vector3D( +delta, -delta, -delta ) ); // Cartesian coordinates
        if ( ! intersects_atoms( crystal_structure, crystal_structure.crystal_lattice().orthogonal_to_fractional( position ), distances2 ) )
            void_spheres_final.push_back( position ); // Cartesian coordinates
        position = Vector3D( void_spheres_intermediate[i] + Vector3D( +delta, -delta, +delta ) ); // Cartesian coordinates
        if ( ! intersects_atoms( crystal_structure, crystal_structure.crystal_lattice().orthogonal_to_fractional( position ), distances2 ) )
            void_spheres_final.push_back( position ); // Cartesian coordinates
        position = Vector3D( void_spheres_intermediate[i] + Vector3D( +delta, +delta, -delta ) ); // Cartesian coordinates
        if ( ! intersects_atoms( crystal_structure, crystal_structure.crystal_lattice().orthogonal_to_fractional( position ), distances2 ) )
            void_spheres_final.push_back( position ); // Cartesian coordinates
        position = Vector3D( void_spheres_intermediate[i] + Vector3D( +delta, +delta, +delta ) ); // Cartesian coordinates
        if ( ! intersects_atoms( crystal_structure, crystal_structure.crystal_lattice().orthogonal_to_fractional( position ), distances2 ) )
            void_spheres_final.push_back( position ); // Cartesian coordinates
    }
    std::cout << "Number of final spheres = " << void_spheres_final.size() << std::endl;

    // Calculate the volume of the void
    double volume = calculate_volume( crystal_structure, void_spheres_final, probe_radius );

//    for ( size_t i( 0 ); i != void_spheres_final.size(); ++i )
//        void_spheres_final[i].show();

    return volume;
}

// ********************************************************************************

// This gives the wrong answer
double find_voids_2( const CrystalStructure & crystal_structure )
{
    double probe_radius( 1.2 );
    // For each atom, precalculate square( VdW radius + probe_radius )
    std::vector< double > distances2;
//    for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
//        distances2.push_back( square( crystal_structure.atom( i ).element().Van_der_Waals_radius() + probe_radius ) );
    for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
        distances2.push_back( square( crystal_structure.atom( i ).element().Van_der_Waals_radius() ) );
    RandomNumberGenerator_double random_number_generator;
    const size_t nprobes( crystal_structure.crystal_lattice().volume() * 100 ); // 5000 sampling points per A^3
    size_t ninside_voids( 0 );
    for ( size_t ix( 0 ); ix != nprobes; ++ix )
    {
        Vector3D random_vector( random_number_generator.next_number(), random_number_generator.next_number(), random_number_generator.next_number() );
        bool intersection_found( false );
        for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
        {
            // Both arguments to CrystalLattice::shortest_distance2() must be in fractional coordinates
            if ( crystal_structure.crystal_lattice().shortest_distance2( crystal_structure.atom( i ).position(), random_vector ) < distances2[i] )
            {
                intersection_found = true;
                break;
            }
        }
        if ( ! intersection_found )
            ++ninside_voids;
    }
    return ( static_cast<double>( ninside_voids ) / static_cast<double>( nprobes ) ) * crystal_structure.crystal_lattice().volume();
}

// ********************************************************************************

double void_volume( const CrystalStructure & crystal_structure )
{
    if ( crystal_structure.natoms() == 0 )
        return crystal_structure.crystal_lattice().volume();
    if ( ! crystal_structure.space_group_symmetry_has_been_applied() )
        throw std::runtime_error( "molecular_volume(): space-group symmetry has not been applied for input crystal structure." );
    RandomNumberGenerator_double random_number_generator;
    const size_t nprobes( ( crystal_structure.crystal_lattice().volume() * 1000.0 ) / crystal_structure.space_group().nsymmetry_operators() ); // 1000 sampling points per A^3
    // For each atom, precalculate square( VdW radius )
    std::vector< double > distances2;
    for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
        distances2.push_back( square( crystal_structure.atom( i ).element().Van_der_Waals_radius() ) );
    size_t ninside_voids( 0 );
    for ( size_t ix( 0 ); ix != nprobes; ++ix )
    {
        Vector3D random_vector( random_number_generator.next_number(), random_number_generator.next_number(), random_number_generator.next_number() );
        random_vector = crystal_structure.crystal_lattice().fractional_to_orthogonal( random_vector );
        if ( ! intersects_atoms( crystal_structure, random_vector, distances2 ) )
                ++ninside_voids;
    }
    return ( static_cast<double>( ninside_voids ) / static_cast<double>( nprobes ) ) * crystal_structure.crystal_lattice().volume();
}

// ********************************************************************************

