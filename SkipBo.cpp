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

#include "SkipBo.h"

#include <stdexcept>

// ********************************************************************************

SkipBoGame::SkipBoGame( const size_t nplayers, const size_t stock_pile_size ):
draw_pile_(1359),
nplayers_(nplayers)
{
    if ( nplayers_ < 2 )
        throw std::runtime_error( "SkipBoGame::SkipBoGame(): number of players must be > 1." );
    // Initialise
    std::vector< size_t > empty_build_pile( 12, 0 );
    build_piles_.push_back( empty_build_pile );
    build_piles_.push_back( empty_build_pile );
    build_piles_.push_back( empty_build_pile );
    build_piles_.push_back( empty_build_pile );
    std::vector< size_t > empty_hand( 5, 0 );
    for ( size_t i(0); i != nplayers_; ++i )
        hands_.push_back( empty_hand );
    // Populate draw pile
// @@    draw_pile_.reserve( 162 );
    for ( size_t i( 0 ); i != 18; ++i )
        draw_pile_.add( 13 );
    for ( size_t i( 0 ); i != 12; ++i )
    {
        for ( size_t j( 1 ); j != 13; ++j )
        {
            draw_pile_.add( j );
        }
    }
    if ( draw_pile_.size() != 162 )
        throw std::runtime_error( "SkipBoGame::SkipBoGame(): draw_pile_.size() != 162." );
    for ( size_t i(0); i != nplayers_; ++i )
    {
        std::vector< size_t > stock_pile;
        stock_pile.reserve( 20 );
        for ( size_t j(0); j != stock_pile_size; ++j )
        {
            stock_pile.push_back( draw_pile_.draw() );
        }
        stock_piles_.push_back( stock_pile );
    }
}

// ********************************************************************************

size_t SkipBoGame::deal_from_draw_pile()
{
    return draw_pile_.draw();
}

// ********************************************************************************

SkipBoPlayer::SkipBoPlayer()
{
}

// ********************************************************************************

void SkipBoPlayer::play( SkipBoGame & skip_bo_game )
{
    // Fill hand up to five from draw pile.
    for ( size_t i(0); i != hand_.size(); ++i )
    {
        if ( hand_[i] == 0 )
        {
            hand_[i] = skip_bo_game.deal_from_draw_pile();
        }
    }
    // Check if stock card fits on a build pile.
    // Check if it is a Skip-Bo card
//    if ( stock )
    {
    }
//    else
    {

    }

    // If no stock card can be played, check if build cards can be used to play stock card.


    //

}

// ********************************************************************************

