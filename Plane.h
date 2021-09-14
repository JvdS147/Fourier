#ifndef PLANE_H
#define PLANE_H

/* *********************************************
Copyright (c) 2013-2021, Cornelis Jan (Jacco) van de Streek
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

class Vector2D;
class Vector3D;

#include "CoordinateFrame.h"
#include "NormalisedVector3D.h"

#include <string>
#include <vector>

/*
  A plane.

  Internally represented by Hesse's normal representation of a plane. This implementation detail is relevant
  for signed_distance() which can be used to characterise the inversion geometry of nitrogen atoms.

*/
class Plane
{
public:

    // Default constructor.
    Plane();

    Plane( const NormalisedVector3D & n, const double c );

    // The order of the points determines the sign of the normal n.
    Plane( const Vector3D & r1, const Vector3D & r2, const Vector3D & r3 );
    
    // Least squares plane.
    // To reduce the number of dependencies, the Plane class does not know
    // about the CollectionOfPoints class. Use the corresponding function in
    // 3DCalculations.
    explicit Plane( const std::vector< Vector3D > & points );

    NormalisedVector3D normal() const { return n_; }

    double constant() const { return c_; }

    double signed_distance( const Vector3D & r ) const;

    double distance( const Vector3D & r ) const;
    
    // Projection along the normal onto the plane, turning a 3D point into a 2D point.
    Vector3D projection3D( const Vector3D & point ) const;

    // Projection along the normal onto the plane, turning a 3D point into a 2D point.
//    Vector2D projection2D( const Vector3D & point ) const;

    // The normal is the z-axis, x-axis and y-axis arbitrary but reproducible.
    CoordinateFrame coordinate_frame() const { return coordinate_frame_; }

    std::string to_string() const;

    void show() const; // For debugging.

// @@ Equality operator: N.B. special case if plane1.normal() == -plane2.normal() and plane1.constant() == -plane2.constant()

private:
    NormalisedVector3D n_;
    double c_;
    CoordinateFrame coordinate_frame_; // Could be cached / mutable, I simply precalculate and store it.
    
    void calculate_coordinate_frame();
};

#endif // PLANE_H

