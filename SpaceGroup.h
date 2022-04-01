#ifndef SPACEGROUP_H
#define SPACEGROUP_H

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

#include "SymmetryOperator.h"

#include <vector>
#include <string>
#include <iosfwd>

class PointGroup;

/*
  A space group.

  For space group P21/c, use the static constructor:

  SpaceGroup space_group = SpaceGroup::P21c();

  For e.g. space group P-1, use:

  SpaceGroup space_group;
  space_group.add_inversion_at_origin(); // Space group is now P-1
  
  It is guaranteed that the identity is the first symmetry operator.

*/
class SpaceGroup
{
public:

    // Default constructor: P1
    SpaceGroup();

    explicit SpaceGroup( const std::vector< SymmetryOperator > & symmetry_operators, const std::string & name = "" );

    static SpaceGroup P21c();

    // Allows quick building of space groups
    void add_inversion_at_origin();

    // Allows quick building of space groups
    // Should not include [ 0, 0, 0 ].
    void add_centring_vectors( const std::vector< Vector3D > & centring_vectors );

    size_t nsymmetry_operators() const { return symmetry_operators_.size(); }

    SymmetryOperator symmetry_operator( const size_t i ) const { return symmetry_operators_[i]; }
    
    std::vector< SymmetryOperator > symmetry_operators() const { return symmetry_operators_; }

// The following would avoid the unnecessary copying of symmetry operators as forced by the member function SymmetryOperator symmetry_operator( const size_t i ) const; .
//    Vector3D apply_symmetry_operator( const size_t i, const Vector3D & vector ) const { return symmetry_operators_[i] * vector; }

    std::string name() const { return name_; }
    void set_name( const std::string & name ) { name_ = name; }

    bool has_inversion() const { return has_inversion_; }
    
    Vector3D position_of_inversion() const { return position_of_inversion_; }
    
    bool has_inversion_at_origin() const { return has_inversion_at_origin_; }
    
    void apply_similarity_transformation( const SymmetryOperator & symmetry_operator );

    // If you change from C-centred to primitive using the standard transformation [ 0.5, 0.5, 0.0, -0.5, 0.5, 0.0, 0.0, 0.0, 1.0 ],
    // the unit cell volume is halved and so is the number of symmetry operators, because half of the
    // symmetry operators transform to the other half. This function removes the duplicate symmetry operators.
    void remove_duplicate_symmetry_operators();

    // The point group of a space group is the sum of all symmetry operators with the translations removed
    PointGroup point_group() const;
    
    // The point group augmented with the inversion.
    PointGroup laue_class() const;

    std::string crystal_system() const;

    // Does not include [ 0.0, 0.0, 0.0 ]
    std::vector< Vector3D > centring_vectors() const { return centring_vectors_; }

    // Returns "P", "A", "B", "C", "I", "R", "F" or "U" for Unknown
    std::string centring() const { return centring_; }

    void print_multiplication_table() const;
    
    void show() const;

private:
    std::vector< SymmetryOperator > symmetry_operators_;
    std::vector< SymmetryOperator > representative_symmetry_operators_;
    std::vector< Vector3D > centring_vectors_;
    std::string centring_; // "P", "A", "B", "C", "I", "R", "F" or "U" for Unknown
    bool has_inversion_;
    bool has_inversion_at_origin_;
    Vector3D position_of_inversion_;
    std::string name_;

    void decompose();

};

std::ostream & operator<<( std::ostream & os, const SpaceGroup & space_group );

// Technically, two space groups are the same if they correspond to the same entry in the IT,
// i.e. regardless of the unit-cell setting.
// This function checks if two space groups are the same except for the order of the symmetry operators.
bool same_symmetry_operators( const SpaceGroup & lhs, const SpaceGroup & rhs );

// Multiplies all combinations (both ways) and checks that the result is in the symmetry operators
// Currently throws if not successful, should probably return bool.
// Should not only check if the result if a symmetry operator that is in the set, but should
// calculate the entire multiplication table which has as an additional condition that each
// element must occur exactly once in each row and in each column.
void check_if_closed( const std::vector< SymmetryOperator > & symmetry_operators );

#endif // SPACEGROUP_H

