#ifndef COLLECTIONOFPOINTS_H
#define COLLECTIONOFPOINTS_H

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

#include "Vector3D.h"

#include <vector>

/*
  Main purpose is two things:

  - Allow for caching of results of calculations in the future (currently not done)
  - Make it easy and efficient to refer to the points with respect to their centre of mass.

  This is ambiguous: does this class store the points or the points w.r.t. their c.o.m.? In other words,
  when calculating a Plane from a CollectionOfPoints object, which of the two sets of points does it use?

  I think we should get rid of the ambiguity and instead let the user call move_to_centre_of_mass(); when that is
  desired. That has one disadvantage: adding more points after a call to move_to_centre_of_mass(); will probably
  not have the intended effect.

  @@ Maybe should be called SetOfPoints
  @@ We then also need SetOfDoubles
  @@ Maybe these should be selected to be low-level classes, so that e.g. Plane is allowed to include SetOfPoints
     (but not the opposite way round).

@@ We could make this more general. The SetOfDoubles class and the SetOfPoints class could have a "reference" or "offset" defined
which is *always* subtracted when using Vector3D point( const size_t i ) const { return points_[i]; }. I.e., we always store all
original points and an offset, but the getters only return the original point minus the offset.
For doubles, a natural offset is e.g. the lowest value, or the average value.

*/
class CollectionOfPoints
{
public:

    // Default constructor
    CollectionOfPoints();

    explicit CollectionOfPoints( const std::vector< Vector3D > & points );

    void add_point( const Vector3D point );
    void add_points( const std::vector< Vector3D > & points );

    // Until I have made this class a low-level class, this member function
    // helps reduce interclass dependencies.
    std::vector< Vector3D > points() const { return points_; }
    std::vector< Vector3D > points_wrt_com() const { return points_wrt_com_; }

    size_t size() const { return points_.size(); }
    void reserve( const size_t value ) { points_.reserve( value ); }

    Vector3D point( const size_t i ) const { return points_[i]; }
    Vector3D point_wrt_com( const size_t i ) const { return points_wrt_com_[i]; }

    void move_to_centre_of_mass();

    Vector3D average() const { return average_; }
    Vector3D centre_of_mass() const { return average_; }

    // To reduce the number of dependencies, the CollectionOfPoints class does not know
    // about the Plane class. Use the corresponding function in 3DCalculations.

private:
    std::vector< Vector3D > points_;
    std::vector< Vector3D > points_wrt_com_;
    Vector3D average_;

    void update();
};

#endif // COLLECTIONOFPOINTS_H

