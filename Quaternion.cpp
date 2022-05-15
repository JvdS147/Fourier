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

#include "Quaternion.h"
#include "BasicMathsFunctions.h"
#include "Matrix3D.h"
#include "Utilities.h"

#include <cmath>
#include <iostream> // For debugging
#include <stdexcept>


// ********************************************************************************

Quaternion::Quaternion():a_(1.0),b_(0.0),c_(0.0),d_(0.0)
{
}

// ********************************************************************************

Quaternion::Quaternion( const double a, const double b, const double c, const double d ):a_(a),b_(b),c_(c),d_(d)
{
    normalise();
}

// ********************************************************************************

Quaternion::Quaternion( const Matrix3D & rotation_matrix )
{
}

// ********************************************************************************

Matrix3D Quaternion::rotation_matrix() const
{
    return Matrix3D();
}

// ********************************************************************************

// The reciprocal value is taken.
void Quaternion::reciprocal()
{
    b_ = -b_;
    c_ = -c_;
    d_ = -d_;    
}

// ********************************************************************************

// Multiplies the quaternion by itself
void Quaternion::square()
{
    *this = Quaternion( a_*a_ - b_*b_ - c_*c_ - d_*d_,
                        2.0*a_*b_,
                        2.0*a_*c_,
                        2.0*a_*d_ );
}

// ********************************************************************************

// Raises the quaternion to n.
void Quaternion::power( const int n )
{
    if ( n == 0 )
    {
        *this = Quaternion();
        return;
    }
    if ( n < 0 )
        reciprocal();
    size_t abs_n = absolute( n );
    if ( abs_n == 1 )
        return;
    if ( abs_n == 2 )
    {
        square();
        return;
    }
    // Many clever things can be done here, e.g. using the binary power method or recursion.
    // Note that we gain something from using the binary method because our square() function is optimised
    // for Quaternion objects and faster than simply x * x (with x a Quaternion).
    // The current implementation is extremely simple.
    Quaternion quaternion( *this );
    for ( size_t i( 1 ); i != abs_n; ++i )
        *this *= quaternion;
}

// ********************************************************************************

Quaternion & Quaternion::operator*=( const Quaternion & rhs )
{
    *this = Quaternion( a()*rhs.a() - b()*rhs.b() - c()*rhs.c() - d()*rhs.d(),
                        a()*rhs.b() + b()*rhs.a() + c()*rhs.d() - d()*rhs.c(),
                        a()*rhs.c() - b()*rhs.d() + c()*rhs.a() + d()*rhs.b(),
                        a()*rhs.d() + b()*rhs.c() - c()*rhs.b() + d()*rhs.a() );
    return *this;
}

// ********************************************************************************

Quaternion & Quaternion::operator/=( const Quaternion & rhs )
{
    Quaternion rhs_2( rhs );
    rhs_2.reciprocal();
    *this *= rhs_2;
    return *this;
}

// ********************************************************************************

bool Quaternion::operator<( const Quaternion & rhs ) const
{
    if ( rhs.a() < a_ )
        return true;
    if ( a_ < rhs.a() )
        return false;
    if ( rhs.b() < b_ )
        return true;
    if ( b_ < rhs.b() )
        return false;
    if ( rhs.c() < c_ )
        return true;
    return false;
}

// ********************************************************************************

std::string Quaternion::to_string() const
{
    return std::string( double2string( a_ ) + " " + double2string( b_ ) + "*i + " + double2string( c_ ) + "*j + " + double2string( d_ ) + "*k" );
}

// ********************************************************************************

void Quaternion::show() const
{
    std::cout << to_string() << std::endl;
}

// ********************************************************************************

void Quaternion::normalise()
{
    double l = a_*a_ + b_*b_ + c_*c_ + d_*d_;
    if ( l < TOLERANCE )
        throw std::runtime_error( "Quaternion::normalise(): zero vector" );
    l = sqrt( l );
    a_ /= l;
    b_ /= l;
    c_ /= l;
    d_ /= l;
}

// ********************************************************************************

bool nearly_equal( const Quaternion lhs, const Quaternion rhs, const double tolerance )
{
    return ( nearly_equal( lhs.a(), rhs.a(), tolerance ) &&
             nearly_equal( lhs.b(), rhs.b(), tolerance ) &&
             nearly_equal( lhs.c(), rhs.c(), tolerance ) &&
             nearly_equal( lhs.d(), rhs.d(), tolerance ) );
}

// ********************************************************************************

