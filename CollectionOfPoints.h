#ifndef COLLECTIONOFPOINTS_H
#define COLLECTIONOFPOINTS_H

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

#include "Vector3D.h"

#include <vector>

/*
  Main purpose is two things:
  
  - Allow for caching of results of calculations in the future (currently not done)
  - Make it easy and efficient to refer to the points with respect to their centre of mass.
*/
class CollectionOfPoints
{
public:

    // Default constructor
    CollectionOfPoints();

    explicit CollectionOfPoints( const std::vector< Vector3D > & points );

    void add_point( const Vector3D point );
    void add_points( const std::vector< Vector3D > & points );

    size_t size() const { return points_.size(); }
    void reserve( const size_t value ) { points_.reserve( value ); }

    Vector3D point( const size_t i ) const { return points_[i]; }
    Vector3D point_wrt_com( const size_t i ) const { return points_wrt_com_[i]; }
    
    Vector3D average() const { return average_; }
    Vector3D centre_of_mass() const { return average_; }

private:
    std::vector< Vector3D > points_;
    std::vector< Vector3D > points_wrt_com_;
    Vector3D average_;

    void update();
};

#endif // COLLECTIONOFPOINTS_H

