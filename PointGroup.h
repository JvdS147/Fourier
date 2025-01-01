#ifndef POINTGROUP_H
#define POINTGROUP_H

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

#include "Matrix3D.h"

#include <vector>
#include <string>
#include <iosfwd>

/*
  A crystallographic point group.

  Not a general point group, but the point group of a crystal. Five-fold rotation or rotation higher than six-fold are not possible.

  The identity is guaranteed to be the first operator.

  The number of symmetry operators is at most 48.
*/
class PointGroup
{
public:

    // Default constructor: C1.
    PointGroup();

    PointGroup( const std::vector< Matrix3D > & symmetry_operators, const std::string & name = "" );
    
    // Allows quick building of point groups.
    void add_inversion();
    
    // Alternatively, we can add a general member function add_symmetry_element() that checks for the presence of the element and that adds it such that the group remains closed.

    // If we had an add_symmetry_element() member function we might end up in a situation where while building the point group it is inconsistent
    // so adding an operator should always complete the multiplication table.

    size_t nsymmetry_operators() const { return symmetry_operators_.size(); }

    // @@ Boundary should be checked
    Matrix3D symmetry_operator( const size_t i ) const { return symmetry_operators_[ i ]; }

    std::string name() const { return name_; }

    bool has_inversion() const { return has_inversion_; }
    
    // Average of the matrices.
    Matrix3D special_position_operator() const;

private:
    std::vector< Matrix3D > symmetry_operators_;
    std::string name_;
    bool has_inversion_;

};

std::ostream & operator<<( std::ostream & os, const PointGroup & point_group );

// This function checks if two point groups are the same except for the order of the symmetry operators.
bool same_symmetry_operators( const PointGroup & lhs, const PointGroup & rhs );

// Multiplies all combinations (both ways) and checks that the result is in the symmetry operators
// Currently throws if not successful, should probably return bool.
// Should not only check if the result is a symmetry operator that is in the set, but should
// calculate the entire multiplication table which has as an additional condition that each
// element must occur exactly once in each row and in each column.
void check_if_closed( const std::vector< Matrix3D > & symmetry_operators );

#endif // POINTGROUP_H

