#ifndef SYMMETRYOPERATOR_H
#define SYMMETRYOPERATOR_H

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

#include "Matrix3D.h"
#include "Vector3D.h"

//class Fraction;

#include <vector>
#include <iosfwd>

/*
  A SymmetryOperator.

  Symmetry operators are canonicalised, so if two symmetry operators are equivalent they are guaranteed to have the same representation.
  In practice this means that the three elements of the translation vector are adjusted to be in the range [ 0, 1 >.

  All of the symmetry must be independent of the existence of molecules.
  
  The translation vector must be three fractions

*/
class SymmetryOperator
{
public:

    // Default constructor: identity operator
    SymmetryOperator();

    // Symmetry operator consisting of a rotation matrix (3x3) and a translation vector (3 elements)
    // The elements of the translational part are normalised to be in the range [0,1>,
    // so this class is only useful for crystal structures
    SymmetryOperator( const Matrix3D & rotation_matrix, const Vector3D & translation_vector );

    // cif style symmetry-operator string, e.g. "0.5+x, -y -1/2, z"
    explicit SymmetryOperator( std::string input );

    Matrix3D rotation() const { return rotation_matrix_; }

    // The elements of the translational part are normalised to be in the range [0,1>.
    Vector3D translation() const { return translation_vector_; }

    // This is the "N" in Grosse-Kunstleve
    int rotation_part_type() const;

    void invert();

    // cif format: "x,y,-z+1/2"
    std::string to_string() const;

private:
    Matrix3D rotation_matrix_;
    Vector3D translation_vector_;

    void canonicalise();
};

std::ostream & operator<<( std::ostream & os, const SymmetryOperator & symmetry_operator );

bool nearly_equal( const SymmetryOperator & lhs, const SymmetryOperator & rhs, const double tolerance = 0.0000001 );

SymmetryOperator operator*( const SymmetryOperator & lhs, const SymmetryOperator & rhs );

// Careful: multiplication from the left and from the right is very different
Vector3D operator*( const SymmetryOperator & symmetry_operator, const Vector3D & vector3D );

#endif // SYMMETRYOPERATOR_H

