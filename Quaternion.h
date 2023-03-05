#ifndef QUATERNION_H
#define QUATERNION_H

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

class Matrix3D; // It would be cleaner if Quaternion did not have to know about Matrix3D, but because 
                // this Quaternion class is not really a general quaternion class but a rotation class,
                // it makes sense to add the conversion to and from a rotation matrix in here.

#include "BasicMathsFunctions.h"

#include <string>

/*
  This is a "unit quaternion" or a "normalised quaternion" or "versor" used for representing rotations, not a general quaternion.
  THere is an ambiguity in that q and -q represent the same rotation (as on a 3D unit sphere). This can be taken care of
  by forcing a_ to be positive (if a_<0 then all components change sign), but that is not done in this class.
*/
class Quaternion
{
public:

    // Default constructor. Returns the identity quaternion.
    Quaternion();
    
    Quaternion( const double a, const double b, const double c, const double d );

    Quaternion( const Matrix3D & rotation_matrix );

    Matrix3D rotation_matrix() const;

    double a() const { return a_; }
    double b() const { return b_; }
    double c() const { return c_; }
    double d() const { return d_; }

    // The reciprocal value is taken.
    void reciprocal();

    // Multiplies the quaternion by itself
    void square();

    // Raises the quaternion to n.
    void power( const int n );

    Quaternion operator*( const Quaternion & rhs ) const { return Quaternion(*this) *= rhs; }
    Quaternion operator/( const Quaternion & rhs ) const { return Quaternion(*this) /= rhs; }

    Quaternion & operator*=( const Quaternion & rhs );
    Quaternion & operator/=( const Quaternion & rhs );

    bool operator< ( const Quaternion & rhs ) const;
    bool operator> ( const Quaternion & rhs ) const { return ( rhs < *this ); }
    bool operator>=( const Quaternion & rhs ) const { return ! ( *this < rhs ); }
    bool operator<=( const Quaternion & rhs ) const { return ! ( rhs < *this ); }

    std::string to_string() const;

    void show() const; // For debugging

private:
    double a_;
    double b_;
    double c_;
    double d_;
    
    void normalise();

};

bool nearly_equal( const Quaternion lhs, const Quaternion rhs, const double tolerance = TOLERANCE );

#endif // QUATERNION_H

