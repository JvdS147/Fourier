#ifndef NORMALISEDVECTOR3D_H
#define NORMALISEDVECTOR3D_H

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

#include <cstddef> // For definition of size_t
#include <string>

/*
  A vector of dimension 3 with a length normalised to 1.0.
*/
class NormalisedVector3D
{
public:

    // Default constructor: [ 1.0, 0.0, 0.0 ]
    NormalisedVector3D();

    // The vector is normalised upon construction
    // throws if length = 0.0
    NormalisedVector3D( const double x, const double y, const double z );

    // This constructor does not exist because it would make it necessary for the NormalisedVector3D to know about the Vector3D class.
    // This has all been moved to the file 3DCalculations
    //explicit NormalisedVector3D( const Vector3D & vector_3d );

    double x() const { return data_[0]; }
    double y() const { return data_[1]; }
    double z() const { return data_[2]; }
    
    double value( const size_t i ) const { return data_[i]; }

    NormalisedVector3D operator-() const;
    NormalisedVector3D operator+() const;

    double norm2() const { return 1.0; }

    double length() const { return 1.0; }

    // Could be const, decided not to make it const to make the signature identical to the signature of Vector3D
    // Alternatively, since this member function makes no sense, we could make it private
    void normalise() {}

    std::string to_string() const;

    void show() const; // For debugging

private:
    double data_[3];

    void normalise_2();
};

// The transposition is implied
double operator*( const NormalisedVector3D & lhs, const NormalisedVector3D & rhs );

NormalisedVector3D orthogonalise( const NormalisedVector3D & n, const NormalisedVector3D & r );

#endif // NORMALISEDVECTOR3D_H

