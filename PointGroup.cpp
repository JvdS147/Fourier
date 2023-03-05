/* *********************************************
Copyright (c) 2013-2023, Cornelis Jan (Jacco) van de Streek
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

#include "PointGroup.h"
#include "BasicMathsFunctions.h"
#include "Utilities.h"

#include <stdexcept>
#include <iostream> // for debugging

// ********************************************************************************

PointGroup::PointGroup()
{
    symmetry_operators_.push_back( Matrix3D() );
    name_ = "C1";
    has_inversion_ = false;
}

// ********************************************************************************

PointGroup::PointGroup( const std::vector< Matrix3D > & symmetry_operators, const std::string & name )
{
    check_if_closed( symmetry_operators );
    has_inversion_ = false;
    // Go through the list of symmetry operators and divide into proper and improper, check for identity and check if there is an inversion.
    bool identity_found( false );
    size_t identity_index;
    Matrix3D inversion( -1.0 );
    for ( size_t i( 0 ); i != symmetry_operators.size(); ++i )
    {
        double determinant = symmetry_operators[i].determinant();
        if ( nearly_equal( determinant, 1.0 ) )
        {
            if ( symmetry_operators[i].is_nearly_the_identity() )
            {
                identity_found = true;
                identity_index = i;
            }
        }
        else if ( nearly_equal( determinant, -1.0 ) )
        {
            if ( nearly_equal( symmetry_operators[i], inversion ) )
                has_inversion_ = true;
        }
        else
            throw std::runtime_error( "PointGroup::PointGroup(): unexpected determinant = " + double2string( determinant ) );
        symmetry_operators_.push_back( symmetry_operators[i] );
    }
    if ( ! identity_found )
        throw std::runtime_error( "PointGroup::PointGroup(): identity not found." );
    if ( identity_index != 0  )
        std::swap( symmetry_operators_[0], symmetry_operators_[identity_index] );
    name_ = name;
}

// ********************************************************************************

void PointGroup::add_inversion()
{
    if ( has_inversion_ )
        return;
    has_inversion_ = true;
    size_t old_size = symmetry_operators_.size();
    for ( size_t i( 0 ); i != old_size; ++i )
        symmetry_operators_.push_back( -1.0*symmetry_operators_[i] );
}

// ********************************************************************************

Matrix3D PointGroup::special_position_operator() const
{
    Matrix3D result;
    for ( size_t i( 1 ); i != symmetry_operators_.size(); ++i )
        result += symmetry_operators_[i];
    result /= symmetry_operators_.size();
    return result;
}

// ********************************************************************************

std::ostream & operator<<( std::ostream & os, const PointGroup & point_group )
{
    for ( size_t i( 0 ); i != point_group.nsymmetry_operators(); ++i )
        os << point_group.symmetry_operator( i );
    return os;
}

// ********************************************************************************

// This function checks if two point groups are the same except for the order of the symmetry operators.
bool same_symmetry_operators( const PointGroup & lhs, const PointGroup & rhs )
{
    if ( lhs.nsymmetry_operators() != rhs.nsymmetry_operators() )
        return false;
    // Technically, the following logic is incorrect if lhs contains the same symmetry operator twice, but that should not happen.
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

void check_if_closed( const std::vector< Matrix3D > & symmetry_operators )
{
    for ( size_t i( 0 ); i != symmetry_operators.size(); ++i )
    {
        for ( size_t j( 0 ); j != symmetry_operators.size(); ++j )
        {
            Matrix3D result = symmetry_operators[i] * symmetry_operators[j];
            bool found( false );
            for ( size_t i( 0 ); i != symmetry_operators.size(); ++i )
            {
                if ( nearly_equal( result, symmetry_operators[i] ) )
                {
                    found = true;
                    continue;
                }
            }
            if ( ! found )
                throw std::runtime_error( "check_if_closed( std::vector< Matrix3D > ): operator not found." );
        }
    }
}

// ********************************************************************************

