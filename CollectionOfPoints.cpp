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

#include "CollectionOfPoints.h"
#include "Vector3DCalculations.h"

//#include <iostream>

// ********************************************************************************

CollectionOfPoints::CollectionOfPoints()
{
}

// ********************************************************************************

CollectionOfPoints::CollectionOfPoints( const std::vector< Vector3D > & points ):
points_(points)
{
    update();
}

// ********************************************************************************

void CollectionOfPoints::add_point( const Vector3D point )
{
    points_.push_back( point );
    update();
}

// ********************************************************************************

void CollectionOfPoints::add_points( const std::vector< Vector3D > & points )
{
    points_.reserve( points_.size() + points.size() );
    for ( std::vector< Vector3D >::const_iterator it(points.begin()); it != points.end(); ++it )
        points_.push_back( *it );
    update();
}

// ********************************************************************************

void CollectionOfPoints::update()
{
    average_ = ::average( points_ );
    points_wrt_com_.clear();
    points_wrt_com_.reserve( points_.size() );
    for ( std::vector< Vector3D >::const_iterator it(points_.begin()); it != points_.end(); ++it )
        points_wrt_com_.push_back( (*it) - average_ );
}

// ********************************************************************************

