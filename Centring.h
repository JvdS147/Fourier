#ifndef CENTRING_H
#define CENTRING_H

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

class Matrix3D;

#include "Vector3D.h"

#include <cstddef> // For definition of size_t
#include <vector>
#include <string>
//#include <iosfwd>

/*
     "P", "A", "B", "C", "I", "R", "F" or "U" for Unknown
    
    "J" for Joke or Jacco, a P1 unit cell with 10 centrings in each direction
    
    [ 0.0, 0.0, 0.0 ] is guaranteed to be the first centring vector. We could have chosen not to add it at all,
    but there is the minor advantage that if you decide that [ 0.0, 0.0, 0.0 ] should not be included,
    if it is not present at construction time, you do not know if that is intended or if that is by mistake.
    Also, the SpaceGroup class does contain the identity, so this class should also do that.
    And if [ 0.0, 0.0, 0.0 ] is included, you can simply do nsymmetry_operators = ncentrings * nrepresentative_symmetry_operators

    R comes in two flavours: obverse and reverse.
    
    D (for Diagonal) is rare.
*/
class Centring
{
public:

    enum CentringType { P, A, B, C, D, F, I, R_OBVERSE, R_REVERSE, U, J };

    // Default constructor
    Centring();

    explicit Centring( const std::vector< Vector3D > & centring_vectors );

    explicit Centring( std::string centring_name );

    size_t size() const { return centring_vectors_.size(); }

    Vector3D centring_vector( const size_t i ) const;

    bool is_primitive() const { return size() == 1; }

    Matrix3D to_primitive() const;

    bool contains( const Vector3D & centring_vector, const double tolerance = 0.000001 ) const;

    std::vector< Vector3D > centring_vectors() const { return centring_vectors_; }

    CentringType centring_type() const { return centring_type_; }

    std::string centring_name() const;

    void show() const;

private:
    std::vector< Vector3D > centring_vectors_;
    CentringType centring_type_;
};

#endif // CENTRING_H

std::string centring_type_to_string( const Centring::CentringType centring_type );

