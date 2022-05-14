/* *********************************************
Copyright (c) 2013-2022, Cornelis Jan (Jacco) van de Streek
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

#include "SpaceGroup.h"
#include "3DCalculations.h"
#include "BasicMathsFunctions.h"
#include "PointGroup.h"
#include "Utilities.h"

#include <cmath>
#include <stdexcept>
#include <iostream> // for debugging

// ********************************************************************************

SpaceGroup::SpaceGroup():name_("P1")
{
    symmetry_operators_.push_back( SymmetryOperator() );
    decompose();
}

// ********************************************************************************

SpaceGroup::SpaceGroup( const std::vector< SymmetryOperator > & symmetry_operators, const std::string & name ):symmetry_operators_(symmetry_operators),name_(name)
{
    if ( symmetry_operators_.empty() )
        throw std::runtime_error( "SpaceGroup::SpaceGroup(): a space group must have at least one symmetry operator." );
    check_if_closed( symmetry_operators_ );
    // Make the identity the first symmetry operator.
    bool identity_found( false );
    for ( size_t i( 0 ); i != symmetry_operators_.size(); ++i )
    {
        if ( symmetry_operators_[i].is_nearly_the_identity() )
        {
            identity_found = true;
            std::swap( symmetry_operators_[0], symmetry_operators_[i] );
            break;
        }
    }
    if ( ! identity_found )
        throw std::runtime_error( "SpaceGroup::SpaceGroup(): identity not found." );
    decompose();
}

// ********************************************************************************

SpaceGroup SpaceGroup::P21c()
{
    std::vector< SymmetryOperator > symmetry_operators;
    symmetry_operators.push_back( SymmetryOperator( std::string( "x,y,z" ) ) );
    symmetry_operators.push_back( SymmetryOperator( std::string( "-x,1/2+y,1/2-z" ) ) );
    SpaceGroup result( symmetry_operators, "P21/c" );
    result.add_inversion_at_origin();
    return result;
}

// ********************************************************************************

SpaceGroup SpaceGroup::C2c()
{
    std::vector< SymmetryOperator > symmetry_operators;
    symmetry_operators.push_back( SymmetryOperator( std::string( "x,y,z" ) ) );
    symmetry_operators.push_back( SymmetryOperator( std::string( "-x,y,1/2-z" ) ) );
    SpaceGroup result( symmetry_operators, "C2/c" );
    std::vector< Vector3D > centring_vectors;
    centring_vectors.push_back( Vector3D() );
    centring_vectors.push_back( Vector3D( 0.5, 0.5, 0.0 ) );
    result.add_centring_vectors( Centring( centring_vectors ) );
    result.add_inversion_at_origin();
    return result;
}

// ********************************************************************************

void SpaceGroup::add_inversion_at_origin()
{
    if ( has_inversion_at_origin_ )
    {
        std::cout << "SpaceGroup::add_inversion_at_origin(): warning: space group already has an inversion at the origin, will not add it again." << std::endl;
        return;
    }
    if ( has_inversion_ )
        std::cout << "SpaceGroup::add_inversion_at_origin(): warning you are adding an inversion to a space group that already has an inversion." << std::endl;
    std::vector< SymmetryOperator > symmetry_operators;
    symmetry_operators.reserve( symmetry_operators_.size() * 2 );
    SymmetryOperator inversion( Matrix3D( -1.0 ), Vector3D() );
    for ( size_t i( 0 ); i != symmetry_operators_.size(); ++i )
    {
        symmetry_operators.push_back( symmetry_operators_[i] );
        symmetry_operators.push_back( inversion * symmetry_operators_[i] );
    }
    symmetry_operators_ = symmetry_operators;
    decompose();
}

// ********************************************************************************

void SpaceGroup::add_centring_vectors( const Centring & centring )
{
    if ( centring_.size() != 1 )
        std::cout << "SpaceGroup::add_centring_vectors(): warning you are adding centring vectors to a space group that already has centring vectors." << std::endl;
    if ( centring.size() == 1 )
    {
        std::cout << "SpaceGroup::add_centring_vectors(): warning you are adding the zero vector, will not add it again." << std::endl;
        return;
    }
    std::vector< Vector3D > centring_vectors = centring.centring_vectors();
    for ( size_t i( 1 ); i != centring_vectors.size(); ++i )
    {
        if ( centring_.contains( centring_vectors[i] ) )
            throw std::runtime_error( "SpaceGroup::add_centring_vectors(): the centring vector you are adding is already present." );
    }
    std::vector< SymmetryOperator > symmetry_operators;
    symmetry_operators.reserve( symmetry_operators_.size() * centring_vectors.size() );
    for ( size_t i( 0 ); i != symmetry_operators_.size(); ++i )
    {
        for ( size_t j( 0 ); j != centring_vectors.size(); ++j )
        {
            symmetry_operators.push_back( SymmetryOperator( Matrix3D(), centring_vectors[j] ) * symmetry_operators_[i] );
        }
    }
    symmetry_operators_ = symmetry_operators;
    decompose();
}

// ********************************************************************************

SymmetryOperator SpaceGroup::symmetry_operator( const size_t i ) const
{
    if ( i < symmetry_operators_.size() )
        return symmetry_operators_[i];
    throw std::runtime_error( "SpaceGroup::symmetry_operator( size_t ): index out of bounds." );
}

// ********************************************************************************

// All elements of the rotation matrix of a standard symmetry operator are -1, 0 or 1.
// All elements of the translation vector are 0, 1/6, 1/4, 1/3, 1/2, 2/3, 3/4 or 5/6.
bool SpaceGroup::contains_non_standard_symmetry_operator() const
{
    for ( size_t i( 0 ); i != symmetry_operators_.size(); ++i )
    {
        if ( symmetry_operators_[i].is_non_standard_symmetry_operator() )
            return true;
    }
    return false;
}

// ********************************************************************************

void SpaceGroup::apply_similarity_transformation( const Matrix3D & matrix )
{
    Matrix3D inverse = matrix;
    inverse.invert();
    std::vector< SymmetryOperator > new_symmetry_operators;
    new_symmetry_operators.reserve( symmetry_operators_.size() );
    for ( size_t i( 0 ); i != symmetry_operators_.size(); ++i )
    {
        SymmetryOperator new_symmetry_operator = matrix * symmetry_operators_[i] * inverse;
        new_symmetry_operators.push_back( new_symmetry_operator );
    }
    symmetry_operators_ = new_symmetry_operators;
    decompose();
}

// ********************************************************************************

void SpaceGroup::apply_similarity_transformation( const SymmetryOperator & symmetry_operator )
{
    SymmetryOperator inverse = symmetry_operator;
    inverse.invert();
    std::vector< SymmetryOperator > new_symmetry_operators;
    new_symmetry_operators.reserve( symmetry_operators_.size() );
    for ( size_t i( 0 ); i != symmetry_operators_.size(); ++i )
    {
        SymmetryOperator new_symmetry_operator = symmetry_operator * symmetry_operators_[i] * inverse;
        new_symmetry_operators.push_back( new_symmetry_operator );
    }
    symmetry_operators_ = new_symmetry_operators;
    decompose();
}

// ********************************************************************************

void SpaceGroup::remove_duplicate_symmetry_operators()
{
    std::vector< SymmetryOperator > new_symmetry_operators;
    new_symmetry_operators.reserve( symmetry_operators_.size() / 2 );
    new_symmetry_operators.push_back( symmetry_operators_[ 0 ] ); // The identity
    for ( size_t i( 1 ); i != symmetry_operators_.size(); ++i )
    {
        bool found( false );
        for ( size_t j( 0 ); j != new_symmetry_operators.size(); ++j )
        {
            if ( nearly_equal( symmetry_operators_[i], new_symmetry_operators[j] ) )
            {
                found = true;
                break;
            }
        }
        if ( ! found )
            new_symmetry_operators.push_back( symmetry_operators_[i] );
    }
    symmetry_operators_ = new_symmetry_operators;
    decompose();
}

// ********************************************************************************

// The point group of a space group is the sum of all symmetry operators with the translations removed.
PointGroup SpaceGroup::point_group() const
{
    std::vector< Matrix3D > rotations;
    if ( has_inversion_ )
        rotations.reserve( 2 * representative_symmetry_operators_.size() );
    else
        rotations.reserve( representative_symmetry_operators_.size() );
    for ( size_t i( 0 ); i != representative_symmetry_operators_.size(); ++i )
    {
        rotations.push_back( representative_symmetry_operators_[i].rotation() );
        if ( has_inversion_ )
            rotations.push_back( -1.0 * representative_symmetry_operators_[i].rotation() );
    }
    PointGroup result( rotations );
    return result;
}

// ********************************************************************************

// The point group augmented with the inversion.
PointGroup SpaceGroup::laue_class() const
{
    PointGroup result = point_group();
    if ( ! has_inversion_ )
        result.add_inversion();
    return result;
}

// ********************************************************************************

std::string SpaceGroup::crystal_system() const
{
    if ( representative_symmetry_operators_.size() == 1 )
        return "triclinic";
    size_t sum_2( 0 );
    size_t sum_3( 0 );
    size_t sum_4( 0 );
    size_t sum_6( 0 );
    for ( size_t i( 0 ); i != representative_symmetry_operators_.size(); ++i )
    {
        // @@ Ugly. std::abs() cannot return an integer type, abs() does not compile on all platforms.
        size_t absolute_rotation_part_type = representative_symmetry_operators_[i].rotation_part_type();
        if ( absolute_rotation_part_type < 0 )
            absolute_rotation_part_type = -absolute_rotation_part_type;
        switch( absolute_rotation_part_type )
        {
            case  2 : ++sum_2; break;
            case  3 : ++sum_3; break;
            case  4 : ++sum_4; break;
            case  6 : ++sum_6; break;
        }
    }
    if ( sum_3 == 8 )
        return "cubic";
    if ( sum_6 == 2 )
        return "hexagonal";
    if ( sum_3 == 2 )
        return "trigonal";
    if ( sum_4 == 2 )
        return "tetragonal";
    if ( sum_2 == 3 )
        return "orthorhombic";
    if ( sum_2 == 1 )
        return "monoclinic";
    throw std::runtime_error( "SpaceGroup::crystal_system(): we should never be here." );
}

// ********************************************************************************

void SpaceGroup::print_multiplication_table() const
{
    for ( size_t i( 0 ); i != symmetry_operators_.size(); ++i )
        std::cout << size_t2string( i, 3, ' ' ) << " " << symmetry_operators_[i].to_string() << std::endl;
    std::cout << "    |";
    for ( size_t i( 0 ); i != symmetry_operators_.size(); ++i )
        std::cout << size_t2string( i, 3, ' ' ) << " ";
    std::cout << std::endl;
    for ( size_t i( 0 ); i != symmetry_operators_.size(); ++i )
    {
        std::cout << size_t2string( i, 3, ' ' ) << " |";
        for ( size_t j( 0 ); j != symmetry_operators_.size(); ++j )
        {
            SymmetryOperator result = symmetry_operators_[i] * symmetry_operators_[j];
            bool found( false );
            for ( size_t k( 0 ); k != symmetry_operators_.size(); ++k )
            {
                if ( nearly_equal( result, symmetry_operators_[k] ) )
                {
                    found = true;
                    std::cout << size_t2string( k, 3, ' ' ) << " ";
                    break;
                }
            }
            if ( ! found )
                throw std::runtime_error( "SpaceGroup::print_multiplication_table(): operator not found." );
        }
        std::cout << std::endl;
    }
}

// ********************************************************************************

// Decomposes the space group into:
// 1. Presence of inversion yes / no and its position.
// 2. List of centring vectors.
// 3. Whatever is left such that multiplication with all of the above would generate the original set of space-group symmetry operators again.
void SpaceGroup::decompose()
{
    representative_symmetry_operators_.clear();
    has_inversion_ = false;
    has_inversion_at_origin_ = false;
    std::vector< Vector3D > centring_vectors;
    std::vector< Vector3D > inversion_translation_vectors;
    Matrix3D inversion( -1.0 );
    for ( size_t i( 0 ); i != symmetry_operators_.size(); ++i )
    {
        double determinant = symmetry_operators_[i].rotation().determinant();
        if ( nearly_equal( determinant, 1.0 ) )
        {
            if ( symmetry_operators_[i].rotation().is_nearly_the_identity() )
                centring_vectors.push_back( symmetry_operators_[i].translation() );
        }
        else if ( nearly_equal( determinant, -1.0 ) )
        {
            if ( nearly_equal( symmetry_operators_[i].rotation(), inversion ) )
            {
                has_inversion_ = true;
                inversion_translation_vectors.push_back( symmetry_operators_[i].translation() );
            }
        }
        else
            throw std::runtime_error( "SpaceGroup::decompose(): unexpected determinant = " + double2string( determinant ) );
    }
    if ( has_inversion_ )
    {
        // Add the three translations x, y, and z (they are all [0,1>) and find the smallest value.
        double smallest_value( 4.0 );
        for ( size_t i( 0 ); i != inversion_translation_vectors.size(); ++i )
        {
            // Basically a Manhattan distance.
            double value = inversion_translation_vectors[i].x() + inversion_translation_vectors[i].y() + inversion_translation_vectors[i].z();
            if ( value < smallest_value )
            {
                smallest_value = value;
                position_of_inversion_ = inversion_translation_vectors[i];
            }
        }
        if ( nearly_zero( smallest_value ) )
            has_inversion_at_origin_ = true;
    }
    // List of representative symmetry operators here.
    for ( size_t i( 0 ); i != symmetry_operators_.size(); ++i )
    {
        bool found( false );
        for ( size_t j( 0 ); j != representative_symmetry_operators_.size(); ++j )
        {
            if ( nearly_equal( symmetry_operators_[i].rotation(),      representative_symmetry_operators_[j].rotation() ) ||
                 nearly_equal( symmetry_operators_[i].rotation(), -1.0*representative_symmetry_operators_[j].rotation() ) )
            {
                 found = true;
                 break;
            }
        }
        if ( ! found )
            representative_symmetry_operators_.push_back( symmetry_operators_[i] );
    }
    // Now determine the centring.
    centring_ = Centring( centring_vectors );
    // Check for duplicates.
    size_t nexpected_symmetry_operators = representative_symmetry_operators_.size() * centring_.size();
    if ( has_inversion_ )
        nexpected_symmetry_operators *= 2;
    if ( nexpected_symmetry_operators != symmetry_operators_.size() )
        throw std::runtime_error( "SpaceGroup::decompose(): number of symmetry operators not consistent." );
}

// ********************************************************************************

void SpaceGroup::show() const
{
    std::cout << "Symmetry operators:" << std::endl;
    for ( size_t i( 0 ); i != symmetry_operators_.size(); ++i )
        std::cout << symmetry_operators_[i].to_string() << std::endl;
    std::cout << "Representative symmetry operators:" << std::endl;
    for ( size_t i( 0 ); i != representative_symmetry_operators_.size(); ++i )
        std::cout << representative_symmetry_operators_[i].to_string() << std::endl;
    std::cout << "Centring vectors:" << std::endl;
    centring_.show();
    if ( has_inversion_ )
    {
        if ( has_inversion_at_origin_ )
            std::cout << "Has inversion at origin" << std::endl;
        else
            std::cout << "Has inversion at: " << position_of_inversion_.to_string() << std::endl;
    }
    else
        std::cout << "No inversion" << std::endl;
    std::cout << "Name: " << name_ << std::endl;
}

// ********************************************************************************

std::ostream & operator<<( std::ostream & os, const SpaceGroup & space_group )
{
    for ( size_t i( 0 ); i != space_group.nsymmetry_operators(); ++i )
        os << space_group.symmetry_operator( i );
    return os;
}

// ********************************************************************************

// Technically, two space groups are the same if they correspond to the same entry in the IT,
// i.e. regardless of the unit-cell setting.
// This function checks if two space groups are the same except for the order of the symmetry operators.
bool same_symmetry_operators( const SpaceGroup & lhs, const SpaceGroup & rhs )
{
    if ( lhs.nsymmetry_operators() != rhs.nsymmetry_operators() )
        return false;
    // Technically, the following logic is incorrect if lhs contains the same symmetry operator twice,
    // but decompose() contains consistency checks that should prevent this.
    for ( size_t i( 0 ); i != lhs.nsymmetry_operators(); ++i )
    {
        bool found( false );
        for ( size_t j( 0 ); j != rhs.nsymmetry_operators(); ++j )
        {
            if ( nearly_equal( lhs.symmetry_operator( i ), rhs.symmetry_operator( j ) ) )
            {
                found = true;
                break;
            }
        }
        if ( ! found )
            return false;
    }
    return true;
}

// ********************************************************************************

void check_if_closed( const std::vector< SymmetryOperator > & symmetry_operators )
{
    for ( size_t i( 0 ); i != symmetry_operators.size(); ++i )
    {
        for ( size_t j( 0 ); j != symmetry_operators.size(); ++j )
        {
            SymmetryOperator result = symmetry_operators[i] * symmetry_operators[j];
            bool found( false );
            for ( size_t k( 0 ); k != symmetry_operators.size(); ++k )
            {
                if ( nearly_equal( result, symmetry_operators[k] ) )
                {
                    found = true;
                    break;
                }
            }
            if ( ! found )
            {
                throw std::runtime_error( "SpaceGroup::check_if_closed( std::vector< SymmetryOperator > ): operator not found." );
            }
        }
    }
}

// ********************************************************************************

