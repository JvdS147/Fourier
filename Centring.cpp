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

#include "Centring.h"
#include "Matrix3D.h"
#include "StringFunctions.h"

#include <stdexcept>
#include <iostream> // for debugging

// ********************************************************************************

Centring::Centring()
{
    centring_vectors_.push_back( Vector3D() );
    centring_type_ = P;
}

// ********************************************************************************

Centring::Centring( const std::vector< Vector3D > & centring_vectors ):centring_vectors_(centring_vectors)
{
    if ( centring_vectors_.empty() )
        throw std::runtime_error( "Centring::Centring( std::vector< Vector3D > ): error: a centring must have at least one centring vector." );
    // Make the zero vector the first centring vector.
    bool zero_vector_found( false );
    for ( size_t i( 0 ); i != centring_vectors_.size(); ++i )
    {
        if ( centring_vectors_[i].nearly_zero() )
        {
            zero_vector_found = true;
            std::swap( centring_vectors_[0], centring_vectors_[i] );
            break;
        }
    }
    if ( ! zero_vector_found )
        throw std::runtime_error( "Centring::Centring(): error: zero vector not found." );
    // Check for duplicates
    for ( size_t i( 0 ); i != centring_vectors_.size(); ++i )
    {
        for ( size_t j( i+1 ); j != centring_vectors_.size(); ++j )
        {
            if ( nearly_equal( centring_vectors_[i], centring_vectors_[j] ) )
                throw std::runtime_error( "Centring::Centring( std::vector< Vector3D > ): error: duplicate vector found." );
        }
    }
    centring_type_ = U; // Unknown
    if ( centring_vectors_.size() == 1 )
        centring_type_ = P;
    else if ( centring_vectors_.size() == 2 ) // C-centred
    {
        if ( nearly_equal( centring_vectors_[1], Vector3D( 0.0, 0.5, 0.5 ) ) )
            centring_type_ = A;
        else if ( nearly_equal( centring_vectors_[1], Vector3D( 0.5, 0.0, 0.5 ) ) )
            centring_type_ = B;
        else if ( nearly_equal( centring_vectors_[1], Vector3D( 0.5, 0.5, 0.0 ) ) )
            centring_type_ = C;
        else if ( nearly_equal( centring_vectors_[1], Vector3D( 0.5, 0.5, 0.5 ) ) )
            centring_type_ = I;
    }
    else if ( centring_vectors_.size() == 3 ) // R-centred
    {
        if ( contains( Vector3D( 2.0/3.0, 1.0/3.0, 1.0/3.0 ) ) && contains( Vector3D( 1.0/3.0, 2.0/3.0, 2.0/3.0 ) ) )
            centring_type_ = R_OBVERSE;
        else if ( contains( Vector3D( 1.0/3.0, 2.0/3.0, 1.0/3.0 ) ) && contains( Vector3D( 2.0/3.0, 1.0/3.0, 2.0/3.0 ) ) )
            centring_type_ = R_REVERSE;
        else if ( contains( Vector3D( 1.0/3.0, 1.0/3.0, 1.0/3.0 ) ) && contains( Vector3D( 2.0/3.0, 2.0/3.0, 2.0/3.0 ) ) )
            centring_type_ = D;
    }
    else if ( centring_vectors_.size() == 4 ) // F-centred
    {
        if ( contains( Vector3D( 0.0, 0.5, 0.5 ) ) && contains( Vector3D( 0.5, 0.0, 0.5 ) ) && contains( Vector3D( 0.5, 0.5, 0.0 ) ) )
            centring_type_ = F;
    }
    else if ( centring_vectors_.size() == 1000 ) // J-centred
        centring_type_ = J;
}

// ********************************************************************************

Centring::Centring( std::string centring_name )
{
    if ( to_upper( centring_name ) != centring_name )
    {
        std::cout << "Centring::Centring( std::string ): Warning: centring name must be uppercase." << std::endl;
        centring_name = to_upper( centring_name );
    }
    if ( centring_name == "U" )
        throw std::runtime_error( "Centring::Centring( std::string ): error: cannot construct unknown centring with this constructor." );
    if ( centring_name == "J" )
    {
        centring_vectors_.reserve( 1000 );
        for ( size_t i1( 0 ); i1 != 10; ++i1 )
        {
            for ( size_t i2( 0 ); i2 != 10; ++i2 )
            {
                for ( size_t i3( 0 ); i3 != 10; ++i3 )
                {
                    centring_vectors_.push_back( Vector3D( i1/10.0, i2/10.0, i3/10.0 ) );
                }
            }
        }
        centring_type_ = J;
        return;
    }
    centring_vectors_.push_back( Vector3D() );
    if ( centring_name == "P" )
        centring_type_ = P;
    else if ( centring_name == "A" )
    {
        centring_vectors_.push_back( Vector3D( 0.0, 0.5, 0.5 ) );
        centring_type_ = A;
    }
    else if ( centring_name == "B" )
    {
        centring_vectors_.push_back( Vector3D( 0.5, 0.0, 0.5 ) );
        centring_type_ = B;
    }
    else if ( centring_name == "C" )
    {
        centring_vectors_.push_back( Vector3D( 0.5, 0.5, 0.0 ) );
        centring_type_ = C;
    }
    else if ( centring_name == "D" )
    {
        centring_vectors_.push_back( Vector3D( 1.0/3.0, 1.0/3.0, 1.0/3.0 ) );
        centring_vectors_.push_back( Vector3D( 2.0/3.0, 2.0/3.0, 2.0/3.0 ) );
        centring_type_ = D;
    }
    else if ( ( centring_name == "R_OBVERSE" ) || ( centring_name == "R" ) )
    {
        centring_vectors_.push_back( Vector3D( 2.0/3.0, 1.0/3.0, 1.0/3.0 ) );
        centring_vectors_.push_back( Vector3D( 1.0/3.0, 2.0/3.0, 2.0/3.0 ) );
        centring_type_ = R_OBVERSE;
    }
    else if ( centring_name == "R_REVERSE" )
    {
        centring_vectors_.push_back( Vector3D( 1.0/3.0, 2.0/3.0, 1.0/3.0 ) );
        centring_vectors_.push_back( Vector3D( 2.0/3.0, 1.0/3.0, 2.0/3.0 ) );
        centring_type_ = R_REVERSE;
    }
    else if ( centring_name == "F" )
    {
        centring_vectors_.push_back( Vector3D( 0.0, 0.5, 0.5 ) );
        centring_vectors_.push_back( Vector3D( 0.5, 0.0, 0.5 ) );
        centring_vectors_.push_back( Vector3D( 0.5, 0.5, 0.0 ) );
        centring_type_ = F;
    }
    else
        throw std::runtime_error( "Centring::Centring( std::string ): error: centring name not recognised." );
}

// ********************************************************************************

Vector3D Centring::centring_vector( const size_t i ) const
{
    if ( i < centring_vectors_.size() )
        return centring_vectors_[i];
    throw std::runtime_error( "Centring::centring_vector( size_t ): error: index out of bounds." );
}

// ********************************************************************************

Matrix3D Centring::to_primitive() const
{
    if ( is_primitive() )
    {
        std::cout << "Centring::to_primitive(): warning: centring is primitive." << std::endl;
        return Matrix3D();
    }
    if ( centring_type() == A )
        return Matrix3D(  1.0,  0.0,  0.0,
                          0.0,  0.5,  0.5,
                          0.0, -0.5,  0.5 );
    if ( centring_type() == B )
        return Matrix3D(  0.5,  0.0,  0.5,
                          0.0,  1.0,  0.0,
                         -0.5,  0.0,  0.5 );
    if ( centring_type() == C )
        return Matrix3D(  0.5,  0.5,  0.0,
                         -0.5,  0.5,  0.0,
                          0.0,  0.0,  1.0 );
    if ( centring_type() == D )
        throw std::runtime_error( "Centring::to_primitive(): centring D not yet implemented." );
    if ( centring_type() == F )
        return Matrix3D(  0.0,  0.5,  0.5,
                          0.5,  0.0,  0.5,
                          0.5,  0.5,  0.0 );
    if ( centring_type() == I )
        return Matrix3D( -0.5,  0.5,  0.5,
                          0.5, -0.5,  0.5,
                          0.5,  0.5, -0.5 );
    if ( centring_type() == R_OBVERSE )
        return Matrix3D( 2.0/3.0, 1.0/3.0, 1.0/3.0,
                         1.0/3.0, 2.0/3.0, 2.0/3.0,
                           0.0,     0.0,     1.0 );
    if ( centring_type() == R_REVERSE )
        return Matrix3D( 1.0/3.0, 2.0/3.0, 2.0/3.0,
                         2.0/3.0, 1.0/3.0, 1.0/3.0,
                           0.0,     0.0,     1.0 );
    if ( centring_type() == U )
        throw std::runtime_error( "Centring::to_primitive(): error: no transformation matrix for centring U." );
    if ( centring_type() == J )
        return Matrix3D(  0.1,  0.0,  0.0,
                          0.0,  0.1,  0.0,
                          0.0,  0.0,  0.1 );
    throw std::runtime_error( "Centring::to_primitive(): centring not yet implemented." );
}

// ********************************************************************************

bool Centring::contains( const Vector3D & centring_vector, const double tolerance ) const
{
    for ( size_t i( 0 ); i != centring_vectors_.size(); ++i )
    {
        if ( nearly_equal( centring_vector, centring_vectors_[i], tolerance ) )
            return true;
    }
    return false;
}

// ********************************************************************************

std::string Centring::centring_name() const
{
    return centring_type_to_string( centring_type_ );
}

// ********************************************************************************

void Centring::show() const
{
    std::cout << centring_name() << std::endl;
    for ( size_t i( 0 ); i != centring_vectors_.size(); ++i )
        std::cout << centring_vectors_[i].to_string() << std::endl;
}

// ********************************************************************************

std::string centring_type_to_string( const Centring::CentringType centring_type )
{
    switch ( centring_type )
    {
        case Centring::P : return "P";
        case Centring::A : return "A";
        case Centring::B : return "B";
        case Centring::C : return "C";
        case Centring::D : return "D";
        case Centring::F : return "F";
        case Centring::I : return "I";
        case Centring::R_OBVERSE : return "R_OBVERSE";
        case Centring::R_REVERSE : return "R_REVERSE";
        case Centring::U : return "U";
        case Centring::J : return "J";
        default : return "Error";
    }
}

// ********************************************************************************

