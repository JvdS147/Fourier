#ifndef DRUNKARDSWALK_H
#define DRUNKARDSWALK_H

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

#include "BagOfNumbers.h"

#include <cmath>

struct GridPoint2D
{
    
    GridPoint2D() { x_ = 0; y_ = 0; }

    GridPoint2D( const int x, const int y ) { x_ = x; y_ = y; }

    int x_;
    int y_;
};

inline size_t Manhattan_distance( const GridPoint2D lhs, const GridPoint2D rhs )
{
    return ( std::abs( rhs.x_ - lhs.x_ ) + std::abs( rhs.y_ - lhs.y_ ) );
}

/*

*/
class DrunkardsWalk
{
public:

    // Default constructor
    // Infinite grid, start position ( 0, 0 )
    explicit DrunkardsWalk( const int idum = 1539 );

    // Infinite grid
    explicit DrunkardsWalk( const GridPoint2D start_position, const int idum = 1539 );

    // Finite grid
    DrunkardsWalk( const GridPoint2D lower_left, const GridPoint2D upper_right, const GridPoint2D start_position, const int idum = 1539 );

    GridPoint2D next_step();
    
    GridPoint2D current_position() const { return current_position_; }

    double displacement() const;

private:
    GridPoint2D start_position_;
    GridPoint2D current_position_;
    bool length_x_is_infinite_;
    bool length_y_is_infinite_;
    GridPoint2D lower_left_;
    GridPoint2D upper_right_;
    BagOfNumbers bag_of_numbers_;
};

#endif // DRUNKARDSWALK_H

