/* *********************************************
Copyright (c) 2013-2025, Cornelis Jan (Jacco) van de Streek
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

#include "SymmetryOperator.h"
#include "3DCalculations.h"
#include "BasicMathsFunctions.h"
#include "Fraction.h"
#include "StringFunctions.h"
#include "Utilities.h"

#include <iostream>
#include <stdexcept>

// ********************************************************************************

SymmetryOperator::SymmetryOperator()
{
}

// ********************************************************************************

SymmetryOperator::SymmetryOperator( const Matrix3D & rotation_matrix, const Vector3D & translation_vector ) :
rotation_matrix_( rotation_matrix ),
translation_vector_( translation_vector )
{
    canonicalise();
}

// ********************************************************************************

SymmetryOperator::SymmetryOperator( std::string input ) :
rotation_matrix_( 0.0 )
{
    std::string original_input( input );
    // Remove all spaces
    input = remove( input, ' ' );
    // Split into three parts separated by commas
    std::vector< std::string > parts;
    // Find the position of the first comma
    size_t iPos1 = input.find_first_of( "," );
    if ( iPos1 == std::string::npos )
        throw std::runtime_error( "SymmetryOperator::SymmetryOperator( std::string ): symmetry line must contain two commas: " + original_input );
    if ( iPos1 == 0 )
        throw std::runtime_error( "SymmetryOperator::SymmetryOperator( std::string ): symmetry line cannot contain empty field: " + original_input );
    parts.push_back( input.substr( 0, iPos1 ) );
    size_t iPos2 = input.find_first_of( ",", iPos1+1 );
    if ( iPos2 == std::string::npos )
        throw std::runtime_error( "SymmetryOperator::SymmetryOperator( std::string ): symmetry line must contain two commas: " + original_input );
    if ( ( iPos2 == iPos1+1 ) || ( iPos2 == input.length()-1 ) )
        throw std::runtime_error( "SymmetryOperator::SymmetryOperator( std::string ): symmetry line cannot contain empty field: " + original_input );
    parts.push_back( input.substr( iPos1+1, iPos2 - (iPos1+1) ) );
    parts.push_back( input.substr( iPos2+1 ) );
    // For each part, separate into parts separated by "+" or "-".
    // This can be made easier by replacing each "-" except the first by "+-" (this turns a binary minus into a unary minus).
    size_t iPos;
    for ( size_t i( 0 ); i != 3; ++i )
    {
        // Find reference to x.
        iPos = parts[i].find_first_of( "xX" );
        if ( iPos != std::string::npos )
        {
            parts[i].erase( iPos, 1 );
            double value = 1.0;
            if ( ( iPos != 0 ) && is_digit( parts[i][iPos-1] ) )
            {
                // We now have something like 2x,y,z , which can happen if we make angles close to 90 degrees after going from centred to primitive.
                value = string2double( parts[i].substr( iPos-1, 1 ) );
                // We only allow for a single digit, so 10x,y,z is still considered an error.
                --iPos;
                parts[i].erase( iPos, 1 );
            }
            if ( ( iPos == 0 ) || ( parts[i][iPos-1] == '+' ) )
                rotation_matrix_.set_value( i, 0, value );
            else if ( parts[i][iPos-1] == '-' )
                rotation_matrix_.set_value( i, 0, -value );
            else
                throw std::runtime_error( "SymmetryOperator::SymmetryOperator( std::string ): \"x\" cannot be preceded by character other than +, - or digit: " + original_input );
            if ( iPos != 0 )
                parts[i].erase( iPos-1, 1 );
        }
        // Find reference to y.
        iPos = parts[i].find_first_of( "yY" );
        if ( iPos != std::string::npos )
        {
            parts[i].erase( iPos, 1 );
            double value = 1.0;
            if ( ( iPos != 0 ) && is_digit( parts[i][iPos-1] ) )
            {
                // We now have something like x,2y,z , which can happen if we make angles close to 90 degrees after going from centred to primitive.
                value = string2double( parts[i].substr( iPos-1, 1 ) );
                // We only allow for a single digit, so x,10y,z is still considered an error.
                --iPos;
                parts[i].erase( iPos, 1 );
            }
            if ( ( iPos == 0 ) || ( parts[i][iPos-1] == '+' ) )
                rotation_matrix_.set_value( i, 1, value );
            else if ( parts[i][iPos-1] == '-' )
                rotation_matrix_.set_value( i, 1, -value );
            else
                throw std::runtime_error( "SymmetryOperator::SymmetryOperator( std::string ): \"y\" cannot be preceded by character other than +, - or digit: " + original_input );
            if ( iPos != 0 )
                parts[i].erase( iPos-1, 1 );
        }
        // Find reference to z.
        iPos = parts[i].find_first_of( "zZ" );
        if ( iPos != std::string::npos )
        {
            parts[i].erase( iPos, 1 );
            double value = 1.0;
            if ( ( iPos != 0 ) && is_digit( parts[i][iPos-1] ) )
            {
                // We now have something like x,y,2z , which can happen if we make angles close to 90 degrees after going from centred to primitive.
                value = string2double( parts[i].substr( iPos-1, 1 ) );
                // We only allow for a single digit, so x,y,10z is still considered an error.
                --iPos;
                parts[i].erase( iPos, 1 );
            }
            if ( ( iPos == 0 ) || ( parts[i][iPos-1] == '+' ) )
                rotation_matrix_.set_value( i, 2, +value );
            else if ( parts[i][iPos-1] == '-' )
                rotation_matrix_.set_value( i, 2, -value );
            else
                throw std::runtime_error( "SymmetryOperator::SymmetryOperator( std::string ): \"z\" cannot be preceded by character other than +, - or digit: " + original_input );
            if ( iPos != 0 )
                parts[i].erase( iPos-1, 1 );
        }
        // Whatever is left must be the translational part.
        double translational_part( 0.0 );
        if ( parts[i].length() != 0 )
        {
            // Is it expressed as a fraction?
            iPos = parts[i].find_first_of( "//" );
            if ( iPos == 0 )
                throw std::runtime_error( "SymmetryOperator::SymmetryOperator( std::string ): fractional part cannot start with // : " + original_input );
            if ( iPos == parts[i].length() )
                throw std::runtime_error( "SymmetryOperator::SymmetryOperator( std::string ): fractional part cannot end with // : " + original_input );
            if ( iPos != std::string::npos )
            {
                int numerator = string2integer( parts[i].substr( 0, iPos ) );
                if ( ( parts[i][iPos+1] == '+' ) || ( parts[i][iPos+1] == '-' ) )
                    throw std::runtime_error( "SymmetryOperator::SymmetryOperator( std::string ): could not parse //+ or //- : " + original_input );
                int denominator = string2integer( parts[i].substr( iPos+1 ) );
                translational_part = static_cast<double>(numerator) / static_cast<double>(denominator);
            }
            else
                translational_part = string2double( parts[i] );
        }
        translation_vector_.set_value( i, translational_part );
    }
    canonicalise();
}

// ********************************************************************************

// This is the "N" in Grosse-Kunstleve.
int SymmetryOperator::rotation_part_type() const
{
    switch( round_to_int( rotation_matrix_.trace() ) )
    {
        case -3 : return -1;
        case -2 : return -6;
        case -1 : return round_to_int( rotation_matrix_.determinant() ) == 1 ? 2 : -4;
        case  0 : return round_to_int( rotation_matrix_.determinant() ) == 1 ? 3 : -3;
        case  1 : return round_to_int( rotation_matrix_.determinant() ) == 1 ? 4 : -2;
        case  2 : return 6;
        case  3 : return 1;
    }
    throw std::runtime_error( "SymmetryOperator::rotation_part_type(): Error." );
}

// ********************************************************************************

// This is wi in Grosse-Kunstleve.
Vector3D SymmetryOperator::intrinsic_translation_part() const
{
    int N = rotation_part_type();
    size_t n = absolute( N );
    if ( N == -1 )
        n = 2;
    else if ( N == -3 )
        n = 6;
    SymmetryOperator result( *this );
    for ( size_t i( 1 ); i != n; ++i )
        result = result * ( *this );
    if ( ! result.rotation().is_nearly_the_identity() )
        throw std::runtime_error( "SymmetryOperator::intrinsic_translation_part(): result is not the identity." );
    return result.translation()/n;
}

// ********************************************************************************

bool SymmetryOperator::has_intrinsic_translation() const
{
    return ( ! nearly_zero( intrinsic_translation_part() ) );
}

// ********************************************************************************

// This is wl in Grosse-Kunstleve
Vector3D SymmetryOperator::location_translation_part() const
{
    return translation() - intrinsic_translation_part();
}

// ********************************************************************************

// All elements of the rotation matrix of a standard symmetry operator are -1, 0 or 1.
// All elements of the translation vector are multiples of 1/12.
bool SymmetryOperator::is_non_standard_symmetry_operator() const
{
    if ( is_non_standard_rotation( rotation_matrix_ ) )
        return true;
    if ( is_non_standard_translation( translation_vector_ ) )
        return true;
    return false;
}

// ********************************************************************************

bool SymmetryOperator::is_nearly_the_identity( const double tolerance ) const
{
    return ( rotation_matrix_.is_nearly_the_identity( tolerance ) && translation_vector_.nearly_zero( tolerance ) );
}

// ********************************************************************************

void SymmetryOperator::invert()
{
    rotation_matrix_.invert();
    translation_vector_ = -1.0 * rotation_matrix_ * translation_vector_;
    canonicalise();
}

// ********************************************************************************

SymmetryOperator & SymmetryOperator::operator*=( const SymmetryOperator & rhs )
{
    *this = SymmetryOperator( rotation() * rhs.rotation(), rotation() * rhs.translation() + translation() );
    return *this;
}

// ********************************************************************************

std::string SymmetryOperator::to_string() const
{
    std::string result;
    std::string xyz( "xyz" );
    for ( size_t row( 0 ); row != 3; ++row )
    {
        bool is_first_character( true );
        for ( size_t col( 0 ); col != 3; ++col )
        {
            int value = round_to_int( rotation_matrix_.value( row, col ) );
            if ( value == 0 )
                continue;
            if ( value < 0 )
                result += "-";
            else
            {
                if ( ! is_first_character )
                    result += "+";
            }
            if ( absolute( value ) > 1 )
                result += size_t2string( absolute( value ) );
            result += xyz[col];
            is_first_character = false;
        }
        // Translational part.
        double value = translation_vector_.value( row );
        // Convert to fraction.
        Fraction fraction = double2fraction( value, Fraction( 1, 12 ) );
        double error = absolute( fraction.to_double() - value );
        // If conversion to a fraction introduces too great an error, use original value.
        if ( error > 0.01 )
        {
            if ( ( value > 0.0 ) && ( ! is_first_character ) )
                result += "+";
            result += double2string( value );
        }
        else
        {
            if ( ! fraction.is_zero() )
            {
                std::string fraction_str = fraction.to_string();
                // In principle, fraction_str can now be "1 + 3/4", so we must remove the spaces. (Actually, the fraction is in the range [0,1> so this is *NOT* possible.).
                fraction_str = remove( fraction_str, ' ' );
                if ( ( fraction > Fraction( 0 ) ) && ( ! is_first_character ) )
                    result += "+";
                result += fraction_str;
            }
        }
        if ( row != 2 )
        result += ",";
    }
    return result;
}

// ********************************************************************************

std::ostream & operator<<( std::ostream & os, const SymmetryOperator & symmetry_operator )
{
    for ( size_t i( 0 ); i != 3; ++i )
        os << symmetry_operator.rotation().value( i, 0 ) << " " <<
              symmetry_operator.rotation().value( i, 1 ) << " " <<
              symmetry_operator.rotation().value( i, 2 ) << " " <<
              symmetry_operator.translation().value( i ) << std::endl;
    return os;
}

// ********************************************************************************

bool nearly_equal( const SymmetryOperator & lhs, const SymmetryOperator & rhs, const double tolerance )
{
    return ( nearly_equal( lhs.rotation(), rhs.rotation(), tolerance ) && nearly_equal( lhs.translation(), rhs.translation(), tolerance ) );
}

// ********************************************************************************

// All elements of the rotation matrix of a standard symmetry operator are -1, 0 or 1.
bool is_non_standard_rotation( const Matrix3D & rotation )
{
    for ( size_t i( 0 ); i != 3; ++i )
    {
        for ( size_t j( 0 ); j != 3; ++j )
        {
            if ( ! ( nearly_zero( rotation.value( i, j ) ) || nearly_equal( absolute ( rotation.value( i, j ) ), 1.0 ) ) )
                return true;
        }
    }
    return false;
}

// ********************************************************************************

// All elements of the translation vector are 0, 1/6, 1/4, 1/3, 1/2, 2/3, 3/4 or 5/6.
bool is_non_standard_translation( const Vector3D & translation )
{
    for ( size_t i( 0 ); i != 3; ++i )
    {
        Fraction fraction = double2fraction( translation.value( i ), Fraction( 1, 12 ) );
        if ( ! nearly_equal( fraction.to_double(), translation.value( i ), 0.01 ) )
            return true;
//        if ( fraction.integer_part() != 0 )
//            return true;
//        if ( fraction.numerator() < 0 )
//            return true;
        if ( fraction == Fraction( 0 ) )
            continue;
        if ( fraction == Fraction( 1, 6 ) )
            continue;
        if ( fraction == Fraction( 1, 4 ) )
            continue;
        if ( fraction == Fraction( 1, 3 ) )
            continue;
        if ( fraction == Fraction( 1, 2 ) )
            continue;
        if ( fraction == Fraction( 2, 3 ) )
            continue;
        if ( fraction == Fraction( 3, 4 ) )
            continue;
        if ( fraction == Fraction( 5, 6 ) )
            continue;
        return true;
    }
    return false;
}

// ********************************************************************************

void SymmetryOperator::canonicalise()
{
    // Check if all rotation entries are -1.0, 0.0 or 1.0
    if ( is_non_standard_symmetry_operator() )
    {
  //      std::cout << "SymmetryOperator::canonicalise(): warning: rotation matrix element is not -1, 0 or 1." << std::endl;
    }
    translation_vector_ = adjust_for_translations( translation_vector_ );
}

// ********************************************************************************

SymmetryOperator operator*( const Matrix3D & matrix, const SymmetryOperator & symmetry_operator )
{
    // Need to be careful with multiplying from the left or from the right.
    return SymmetryOperator( matrix * symmetry_operator.rotation(), matrix * symmetry_operator.translation() );
}

// ********************************************************************************

SymmetryOperator operator*( const SymmetryOperator & symmetry_operator, const Matrix3D & matrix )
{
    // Need to be careful with multiplying from the left or from the right.
    return SymmetryOperator( symmetry_operator.rotation() * matrix, symmetry_operator.translation() );
}

// ********************************************************************************

SymmetryOperator operator*( const SymmetryOperator & lhs, const SymmetryOperator & rhs )
{
    // Need to be careful with multiplying from the left or from the right.
    return SymmetryOperator( lhs.rotation() * rhs.rotation(), lhs.rotation() * rhs.translation() + lhs.translation() );
}

// ********************************************************************************

Vector3D operator*( const SymmetryOperator & symmetry_operator, const Vector3D & vector3D )
{
    return ( symmetry_operator.rotation() * vector3D ) + symmetry_operator.translation();
}

// ********************************************************************************

