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

#include "DrunkardsWalk.h"
#include "BasicMathsFunctions.h"

#include <cmath>
#include <stdexcept>

// ********************************************************************************

DrunkardsWalk::DrunkardsWalk( const int idum ):bag_of_numbers_( 4, idum )
{
    current_position_ = start_position_;
    length_x_is_infinite_ = true;
    length_y_is_infinite_ = true;
}

// ********************************************************************************

// Infinite grid
DrunkardsWalk::DrunkardsWalk( const GridPoint2D start_position, const int idum ):bag_of_numbers_( 4, idum )
{
    start_position_ = start_position;
    current_position_ = start_position_;
    length_x_is_infinite_ = true;
    length_y_is_infinite_ = true;
}

// ********************************************************************************

// Finite grid
DrunkardsWalk::DrunkardsWalk( const GridPoint2D lower_left, const GridPoint2D upper_right, const GridPoint2D start_position, const int idum ):bag_of_numbers_(4, idum)
{
    start_position_ = start_position;
    current_position_ = start_position_;
    length_x_is_infinite_ = false;
    length_y_is_infinite_ = false;
    lower_left_ = lower_left;
    upper_right_ = upper_right;
    if ( current_position_.x_ < lower_left_.x_ )
        throw std::runtime_error( "DrunkardsWalk::DrunkardsWalk( GridPoint2D, GridPoint2D, GridPoint2D ): boundaries inconsistent" );
    if ( current_position_.x_ > upper_right_.x_ )
        throw std::runtime_error( "DrunkardsWalk::DrunkardsWalk( GridPoint2D, GridPoint2D, GridPoint2D ): boundaries inconsistent" );
    if ( current_position_.y_ < lower_left_.y_ )
        throw std::runtime_error( "DrunkardsWalk::DrunkardsWalk( GridPoint2D, GridPoint2D, GridPoint2D ): boundaries inconsistent" );
    if ( current_position_.y_ > upper_right_.y_ )
        throw std::runtime_error( "DrunkardsWalk::DrunkardsWalk( GridPoint2D, GridPoint2D, GridPoint2D ): boundaries inconsistent" );
    if ( lower_left_.x_ > upper_right_.x_ )
        throw std::runtime_error( "DrunkardsWalk::DrunkardsWalk( GridPoint2D, GridPoint2D, GridPoint2D ): boundaries inconsistent" );
    if ( lower_left_.y_ > upper_right_.y_ )
        throw std::runtime_error( "DrunkardsWalk::DrunkardsWalk( GridPoint2D, GridPoint2D, GridPoint2D ): boundaries inconsistent" );
    if ( ( lower_left_.x_ == upper_right_.x_ ) && ( lower_left_.y_ == upper_right_.y_ ) )
        throw std::runtime_error( "DrunkardsWalk::DrunkardsWalk( GridPoint2D, GridPoint2D, GridPoint2D ): boundaries inconsistent" );
}

// ********************************************************************************

GridPoint2D DrunkardsWalk::next_step()
{
    bool valid_move( false );
    do
    {
        size_t random_direction = bag_of_numbers_.draw_with_replace();
        // 0 = -x, 1 = -y, 2 = +x, 3 = +y
        if ( random_direction == 0 )
        {
            if ( length_x_is_infinite_ || ( current_position_.x_ > lower_left_.x_ ) )
            {
                --current_position_.x_;
                valid_move = true;
            }
        }
        else if ( random_direction == 1 )
        {
            if ( length_y_is_infinite_ || ( current_position_.y_ > lower_left_.y_ ) )
            {
                --current_position_.y_;
                valid_move = true;
            }
        }
        else if ( random_direction == 2 )
        {
            if ( length_x_is_infinite_ || ( current_position_.x_ < upper_right_.x_ ) )
            {
                ++current_position_.x_;
                valid_move = true;
            }
        }
        else if ( random_direction == 3 )
        {
            if ( length_y_is_infinite_ || ( current_position_.y_ < upper_right_.y_ ) )
            {
                ++current_position_.y_;
                valid_move = true;
            }
        }
    } while ( ! valid_move );
    return current_position_;
}

// ********************************************************************************

double DrunkardsWalk::displacement() const
{
    return sqrt( square( static_cast<double>(current_position_.x_) - static_cast<double>(start_position_.x_) ) +
                 square( static_cast<double>(current_position_.y_) - static_cast<double>(start_position_.y_) ) );
}

// ********************************************************************************

