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

#include "ModelBuilding.h"
#include "3DCalculations.h"
#include "CrystalStructure.h"
#include "NormalisedVector3D.h"
#include "Plane.h"
#include "Utilities.h"
#include "Vector3D.h"

#include <stdexcept>
#include <iostream>

// ********************************************************************************

// Everything in Cartesian coordinates
std::vector< Vector3D > add_methyl_group( const Vector3D & atom_1, const Vector3D & atom_2, const Angle angle )
{
    // The following generates an orthonormal basis based on one direction.
    Vector3D difference_vector_temp = atom_1 - atom_2;
    NormalisedVector3D basis_vector_1 = normalised_vector( difference_vector_temp );
    NormalisedVector3D basis_vector_2;
    NormalisedVector3D basis_vector_3;
    generate_basis_1( basis_vector_1, basis_vector_2, basis_vector_3 );

    double XH_bond_length( 1.089 );
    double bx = XH_bond_length * Angle::from_degrees( 110.5 - 90.0 ).sine();
    double projection = XH_bond_length * Angle::from_degrees( 110.5 - 90.0 ).cosine();
    std::vector< Vector3D > result;
    result.push_back( atom_1 + bx * basis_vector_1 + projection *                                    angle.sine() * basis_vector_2 + projection *                                    angle.cosine() * basis_vector_3 );
    result.push_back( atom_1 + bx * basis_vector_1 + projection * ( angle + Angle::from_degrees( 120.0 ) ).sine() * basis_vector_2 + projection * ( angle + Angle::from_degrees( 120.0 ) ).cosine() * basis_vector_3 );
    result.push_back( atom_1 + bx * basis_vector_1 + projection * ( angle + Angle::from_degrees( 240.0 ) ).sine() * basis_vector_2 + projection * ( angle + Angle::from_degrees( 240.0 ) ).cosine() * basis_vector_3 );
    return result;
}

// ********************************************************************************

// Everything in Cartesian coordinates
std::vector< Vector3D > add_methyl_group( const Vector3D & atom_1, const Vector3D & atom_2, const Vector3D & atom_3 )
{
    Vector3D difference_vector_temp = atom_1 - atom_2;
    NormalisedVector3D basis_vector_1 = normalised_vector( difference_vector_temp );
    Vector3D r = atom_3 - atom_2;
    NormalisedVector3D basis_vector_2 = orthogonalise( basis_vector_1, r );
    NormalisedVector3D basis_vector_3;
    generate_basis_2( basis_vector_1, basis_vector_2, basis_vector_3 );

    double XH_bond_length( 1.089 );
    double bx = XH_bond_length * Angle::from_degrees( 110.5 - 90.0 ).sine();
    double projection = XH_bond_length * Angle::from_degrees( 110.5 - 90.0 ).cosine();
    std::vector< Vector3D > result;
    result.push_back( atom_1 + bx * basis_vector_1 - projection *                                         basis_vector_2                                                                     );
    result.push_back( atom_1 + bx * basis_vector_1 - projection * Angle::from_degrees( 120.0 ).cosine() * basis_vector_2 - projection * Angle::from_degrees( 120.0 ).sine() * basis_vector_3 );
    result.push_back( atom_1 + bx * basis_vector_1 - projection * Angle::from_degrees( 240.0 ).cosine() * basis_vector_2 - projection * Angle::from_degrees( 240.0 ).sine() * basis_vector_3 );
    return result;
}

// ********************************************************************************

// Everything in Cartesian coordinates
std::vector< Vector3D > add_hydrogen_atoms( const Vector3D & atom_1, const Vector3D & atom_2, const size_t nhydrogens, const Angle angle )
{
    // The following generates an orthonormal basis based on one direction.
    Vector3D difference_vector_temp = atom_1 - atom_2;
    NormalisedVector3D basis_vector_1 = normalised_vector( difference_vector_temp );
    NormalisedVector3D basis_vector_2;
    NormalisedVector3D basis_vector_3;
    generate_basis_1( basis_vector_1, basis_vector_2, basis_vector_3 );
    double target_bond_length( 1.0 );
    double bx = target_bond_length * Angle::from_degrees( 110.5 - 90.0 ).sine();
    double projection = target_bond_length * Angle::from_degrees( 110.5 - 90.0 ).cosine();
    std::vector< Vector3D > result;
    if ( nhydrogens == 0 )
        throw std::runtime_error( "add_hydrogen_atoms() : zero hydrogen atoms requested." );
    if ( nhydrogens > 6 )
        std::cout << "add_hydrogen_atoms() : Warning: more than six hydrogen atoms requested." << std::endl;
    for ( size_t i( 0 ); i != nhydrogens; ++i )
        result.push_back( atom_1 + bx * basis_vector_1 + projection * ( angle + Angle::from_degrees( i * ( 360.0 / nhydrogens ) ) ).sine() * basis_vector_2 + projection * ( angle + Angle::from_degrees( i * ( 360.0 / nhydrogens ) ) ).cosine() * basis_vector_3 );
    return result;
}

// ********************************************************************************

// Everything in Cartesian coordinates
std::vector< Vector3D > add_2_hydrogen_atoms_to_sp3_nitrogen( const Vector3D & atom_N, const Vector3D & atom_2, const Vector3D & atom_H )
{
    Vector3D difference_vector_temp = atom_N - atom_2;
    NormalisedVector3D basis_vector_1 = normalised_vector( difference_vector_temp );
    NormalisedVector3D basis_vector_2;
    NormalisedVector3D basis_vector_3;
    generate_basis_1( basis_vector_1, basis_vector_2, basis_vector_3 );
    double bx = 1.015 * Angle::from_degrees( 111.0 - 90.0 ).sine();
    double projection = 1.015 * Angle::from_degrees( 111.0 - 90.0 ).cosine();
    std::vector< Vector3D > result;
    difference_vector_temp = atom_H - atom_N;
    if ( difference_vector_temp.length() < 0.001 )
        throw std::runtime_error( "add_2_hydrogen_atoms_to_sp3_nitrogen() : atoms are too close together." );
    NormalisedVector3D H_direction = normalised_vector( difference_vector_temp );
    double sine_phi = H_direction * basis_vector_2;
    double cosine_phi = H_direction * basis_vector_3;
    Angle angle = ATAN2( sine_phi, cosine_phi );
    result.push_back( atom_N + bx * basis_vector_1 + projection * ( angle + Angle::from_degrees( 120.0 ) ).sine() * basis_vector_2 + projection * ( angle + Angle::from_degrees( 120.0 ) ).cosine() * basis_vector_3 );
    result.push_back( atom_N + bx * basis_vector_1 + projection * ( angle + Angle::from_degrees( 240.0 ) ).sine() * basis_vector_2 + projection * ( angle + Angle::from_degrees( 240.0 ) ).cosine() * basis_vector_3 );
    return result;
}

// ********************************************************************************

// Add two hydrogen atoms to an sp2 nitrogen atom
// Returns the coordinates of the two hydrogen atoms
// atom_1 and atom_2 are required to define the plane in which the two atoms should lie.
// Everything in Cartesian coordinates
std::vector< Vector3D > add_2_hydrogen_atoms_to_sp2_nitrogen( const Vector3D & atom_N, const Vector3D & atom_bound_to_N, const Vector3D & atom_1, const Vector3D & atom_2 )
{
    double target_bond_length = 1.015;
    // The H--N--H angle is 120.0. That makes the H--N--X angle also 120.0.
    std::vector< Vector3D > points;
    points.push_back( atom_N );
    points.push_back( atom_bound_to_N );
    points.push_back( atom_1 );
    points.push_back( atom_2 );
    Plane plane( points );
    Vector3D NminC = atom_N - atom_bound_to_N;
    NminC /= NminC.length();
    Vector3D p1 =  cross_product( NminC, NormalisedVector3D2Vector3D( plane.normal() ) );
    p1 /= p1.length();
    Vector3D H_atom_1 = atom_N + ( NminC * 0.5 * target_bond_length ) + ( p1 * 0.5 * sqrt(3.0) * target_bond_length );
    Vector3D H_atom_2 = atom_N + ( NminC * 0.5 * target_bond_length ) - ( p1 * 0.5 * sqrt(3.0) * target_bond_length );
    std::vector< Vector3D > result;
    result.push_back( H_atom_1 );
    result.push_back( H_atom_2 );
    return result;
}

// ********************************************************************************

// Everything in Cartesian coordinates
Vector3D add_hydrogen_atom_to_sp2_atom( const Vector3D & central_atom, const Element element_central_atom, const Vector3D & neighbour_1, const Vector3D & neighbour_2 )
{
    Vector3D difference_vector_1 = neighbour_1 - central_atom;
    if ( difference_vector_1.length() < 0.001 )
        throw std::runtime_error( "add_hydrogen_atom_to_sp2_atom() : atoms are too close together 1." );
    difference_vector_1.set_length( 1.0 );
    Vector3D difference_vector_2 = neighbour_2 - central_atom;
    if ( difference_vector_2.length() < 0.001 )
        throw std::runtime_error( "add_hydrogen_atom_to_sp2_atom() : atoms are too close together 2." );
    difference_vector_2.set_length( 1.0 );
    Vector3D average_vector = ( difference_vector_1 + difference_vector_2 ) / 2.0;
    average_vector *= -1.0;
    double target_bond_length( 1.0 );
    if ( element_central_atom == Element( "B" ) )
        target_bond_length = 1.215;
    else if ( element_central_atom == Element( "C" ) )
        target_bond_length = 1.090;
    else if ( element_central_atom == Element( "N" ) )
        target_bond_length = 1.015;
    else if ( element_central_atom == Element( "O" ) )
        target_bond_length = 0.993;
    average_vector.set_length( target_bond_length );
    return central_atom + average_vector;
}

// ********************************************************************************

// Everything in Cartesian coordinates
Vector3D add_hydrogen_atom_to_sp3_atom( const Vector3D & central_atom, const Element element_central_atom, const Vector3D & neighbour_1, const Vector3D & neighbour_2, const Vector3D & neighbour_3 )
{
    Vector3D v1 = neighbour_1 - central_atom;
    v1.set_length( 1.0 );
    Vector3D v2 = neighbour_2 - central_atom;
    v2.set_length( 1.0 );
    Vector3D v3 = neighbour_3 - central_atom;
    v3.set_length( 1.0 );
    Vector3D average_vector = ( v1 + v2 + v3 ) / 3.0;
    average_vector *= -1.0;
    double target_bond_length( 1.0 );
    if ( element_central_atom == Element( "B" ) )
        target_bond_length = 1.215;
    else if ( element_central_atom == Element( "C" ) )
        target_bond_length = 1.090;
    else if ( element_central_atom == Element( "N" ) )
        target_bond_length = 1.015;
    else if ( element_central_atom == Element( "O" ) )
        target_bond_length = 0.993;
    average_vector.set_length( target_bond_length );
    return central_atom + average_vector;
}

// ********************************************************************************

// Everything in Cartesian coordinates
std::vector< Vector3D > add_2_hydrogen_atoms_to_sp3_atom( const Vector3D & central_atom, const Element element_central_atom, const Vector3D & neighbour_1, const Vector3D & neighbour_2 )
{
    Vector3D v1 = neighbour_1 - central_atom;
    v1.set_length( 1.0 );
    Vector3D v2 = neighbour_2 - central_atom;
    v2.set_length( 1.0 );
    Vector3D average_vector = ( v1 + v2 ) / 2.0;
    average_vector *= -1.0;
    NormalisedVector3D basis_vector_1( normalised_vector( average_vector ) );
    NormalisedVector3D basis_vector_2( normalised_vector( cross_product( v1, v2 ) ) );
    double target_bond_length( 1.0 );
    if ( element_central_atom == Element( "B" ) )
        target_bond_length = 1.215;
    else if ( element_central_atom == Element( "C" ) )
        target_bond_length = 1.090;
    else if ( element_central_atom == Element( "N" ) )
        target_bond_length = 1.015;
    else if ( element_central_atom == Element( "O" ) )
        target_bond_length = 0.993;
    Vector3D H_atom_1 = central_atom + target_bond_length * Angle::from_degrees( 110.5 / 2.0 ).cosine() * basis_vector_1 + target_bond_length * Angle::from_degrees( 110.5 / 2.0 ).sine() * basis_vector_2;
    Vector3D H_atom_2 = central_atom + target_bond_length * Angle::from_degrees( 110.5 / 2.0 ).cosine() * basis_vector_1 - target_bond_length * Angle::from_degrees( 110.5 / 2.0 ).sine() * basis_vector_2;
    std::vector< Vector3D > result;
    result.push_back( H_atom_1 );
    result.push_back( H_atom_2 );
    return result;
}

// ********************************************************************************

void normalise_X_H_bonds( CrystalStructure & crystal_structure )
{
    // Loop over all atoms.
    for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
    {
        // If not hydrogen, continue.
        if ( ! crystal_structure.atom( i ).element().is_H_or_D() )
            continue;
        // Find the non-hydrogen atom nearest to it.
        size_t smallest_distance_index = 0;
        double smallest_distance2 = 10.0;
        for ( size_t j( 0 ); j != crystal_structure.natoms(); ++j )
        {
            // If hydrogen, continue.
            if ( crystal_structure.atom( j ).element().is_H_or_D() )
                continue;
            double distance2 = crystal_structure.crystal_lattice().shortest_distance2( crystal_structure.atom( j ).position(), crystal_structure.atom( i ).position() );
            if ( distance2 < smallest_distance2 )
            {
                smallest_distance2 = distance2;
                smallest_distance_index = j;
            }
        }
        double smallest_distance = sqrt( smallest_distance2 );
        if ( smallest_distance < 0.001 )
            throw std::runtime_error( "normalise_X_H_bonds(): points too close together." );
        if ( smallest_distance > 3.0 )
            throw std::runtime_error( "normalise_X_H_bonds(): atoms not bound." );
        // We need to be able to set the length of the X-H bond, so we must work in Cartesian coordinates.
        Vector3D difference_frac;
        double distance;
        crystal_structure.shortest_distance( crystal_structure.atom( smallest_distance_index ).position(), crystal_structure.atom( i ).position(), distance, difference_frac );
        if ( ! nearly_equal( distance, smallest_distance ) )
            throw std::runtime_error( "normalise_X_H_bonds() : distances differ." );
        Vector3D difference_cart = crystal_structure.crystal_lattice().fractional_to_orthogonal( difference_frac );
        double target_bond_length( 1.0 );
        if ( crystal_structure.atom( smallest_distance_index ).element() == Element( "B" ) )
            target_bond_length = 1.215;
        else if ( crystal_structure.atom( smallest_distance_index ).element() == Element( "C" ) )
            target_bond_length = 1.090;
        else if ( crystal_structure.atom( smallest_distance_index ).element() == Element( "N" ) )
            target_bond_length = 1.015;
        else if ( crystal_structure.atom( smallest_distance_index ).element() == Element( "O" ) )
            target_bond_length = 0.993;
        difference_cart *= target_bond_length / distance;
        difference_frac = crystal_structure.crystal_lattice().orthogonal_to_fractional( difference_cart );
        Vector3D H_atom_frac = crystal_structure.atom( smallest_distance_index ).position() + difference_frac;
        Atom new_atom = crystal_structure.atom( i );
        new_atom.set_position( H_atom_frac );
        crystal_structure.set_atom( i, new_atom );
    }
}

// ********************************************************************************

void normalise_C_F_bonds( CrystalStructure & crystal_structure )
{
    // Loop over all atoms.
    for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
    {
        // If not fluorine, continue.
        if ( crystal_structure.atom( i ).element() != Element( "F" ) )
            continue;
        // Find the carbon atom nearest to it.
        size_t smallest_distance_index = 0;
        double smallest_distance2 = 10.0;
        for ( size_t j( 0 ); j != crystal_structure.natoms(); ++j )
        {
            // If not carbon, continue.
            if ( crystal_structure.atom( j ).element() != Element( "C" ) )
                continue;
            double distance2 = crystal_structure.crystal_lattice().shortest_distance2( crystal_structure.atom( j ).position(), crystal_structure.atom( i ).position() );
            if ( distance2 < smallest_distance2 )
            {
                smallest_distance2 = distance2;
                smallest_distance_index = j;
            }
        }
        double smallest_distance = sqrt( smallest_distance2 );
        if ( smallest_distance < 0.001 )
            throw std::runtime_error( "normalise_C_F_bonds(): points too close together." );
        if ( smallest_distance > 3.0 )
            throw std::runtime_error( "normalise_C_F_bonds(): atoms not bound." );
        // We need to be able to set the length of the C-F bond, so we must work in Cartesian coordinates.
        Vector3D difference_frac;
        double distance;
        crystal_structure.shortest_distance( crystal_structure.atom( smallest_distance_index ).position(), crystal_structure.atom( i ).position(), distance, difference_frac );
        if ( ! nearly_equal( distance, smallest_distance ) )
            throw std::runtime_error( "normalise_C_F_bonds() : distances differ." );
        Vector3D difference_cart = crystal_structure.crystal_lattice().fractional_to_orthogonal( difference_frac );
        double target_bond_length( 1.38 );
        difference_cart *= target_bond_length / distance;
        difference_frac = crystal_structure.crystal_lattice().orthogonal_to_fractional( difference_cart );
        Vector3D F_atom_frac = crystal_structure.atom( smallest_distance_index ).position() + difference_frac;
        Atom new_atom = crystal_structure.atom( i );
        new_atom.set_position( F_atom_frac );
        crystal_structure.set_atom( i, new_atom );
    }
}

// ********************************************************************************

void normalise_Cu_Cl_bonds( CrystalStructure & crystal_structure )
{
    // Loop over all atoms.
    for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
    {
        // If not chlorine, continue.
        if ( crystal_structure.atom( i ).element() != Element( "Cl" ) )
            continue;
        // Find the Cu atom nearest to it.
        size_t smallest_distance_index = 0;
        double smallest_distance2 = 10.0;
        for ( size_t j( 0 ); j != crystal_structure.natoms(); ++j )
        {
            // If not Cu, continue.
            if ( crystal_structure.atom( j ).element() != Element( "Cu" ) )
                continue;
            double distance2 = crystal_structure.crystal_lattice().shortest_distance2( crystal_structure.atom( j ).position(), crystal_structure.atom( i ).position() );
            if ( distance2 < smallest_distance2 )
            {
                smallest_distance2 = distance2;
                smallest_distance_index = j;
            }
        }
        double smallest_distance = sqrt( smallest_distance2 );
        if ( smallest_distance < 0.001 )
            throw std::runtime_error( "normalise_Cu_Cl_bonds(): points too close together." );
        if ( smallest_distance > 3.0 )
            throw std::runtime_error( "normalise_Cu_Cl_bonds(): atoms not bound." );
        // We need to be able to set the length of the Cu-Cl bond, so we must work in Cartesian coordinates.
        Vector3D difference_frac;
        double distance;
        crystal_structure.shortest_distance( crystal_structure.atom( smallest_distance_index ).position(), crystal_structure.atom( i ).position(), distance, difference_frac );
        if ( ! nearly_equal( distance, smallest_distance ) )
            throw std::runtime_error( "normalise_Cu_Cl_bonds() : distances differ." );
        Vector3D difference_cart = crystal_structure.crystal_lattice().fractional_to_orthogonal( difference_frac );
        double target_bond_length( 2.25 );
        difference_cart *= target_bond_length / distance;
        difference_frac = crystal_structure.crystal_lattice().orthogonal_to_fractional( difference_cart );
        Vector3D Cl_atom_frac = crystal_structure.atom( smallest_distance_index ).position() + difference_frac;
        Atom new_atom = crystal_structure.atom( i );
        new_atom.set_position( Cl_atom_frac );
        crystal_structure.set_atom( i, new_atom );
    }
}

// ********************************************************************************
