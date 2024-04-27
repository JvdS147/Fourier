#ifndef COORDINATEFRAME_H
#define COORDINATEFRAME_H

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

#include "NormalisedVector3D.h"

// A CrystalLattice is a crystallographic unit cell, a CoordinateFrame is a Cartesian coordinate frame
class CoordinateFrame
{
public:

    CoordinateFrame();

    // @@ Because I expect this overload to get called with the normal to a plane, the direction that is provided would most naturally be the z-axis
    // Generates three orthonormal basis vectors based on a single direction. The direction can be e.g. the normal to a plane,
    // the basis within the plane is then chosen arbitrarily but reproducibly.    
    explicit CoordinateFrame( const NormalisedVector3D & x_axis );
    
    CoordinateFrame( const NormalisedVector3D & x_axis, const NormalisedVector3D & y_axis );
    
    CoordinateFrame( const NormalisedVector3D & x_axis, const NormalisedVector3D & y_axis, const NormalisedVector3D & z_axis );

    NormalisedVector3D x_axis() const { return x_axis_; }
    NormalisedVector3D y_axis() const { return y_axis_; }
    NormalisedVector3D z_axis() const { return z_axis_; }

    void print() const;

    void show() const; // For debugging

private:
    NormalisedVector3D x_axis_;
    NormalisedVector3D y_axis_;
    NormalisedVector3D z_axis_;
};

//bool nearly_equal( const CoordinateFrame & lhs, const CoordinateFrame & rhs, double tolerance = TOLERANCE );

#endif // COORDINATEFRAME_H

