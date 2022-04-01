#ifndef SKIPBO_H
#define SKIPBO_H

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

#include "BagOfNumbers.h"

#include <vector>

/*
 * Skip-Bo game class.
 *
 * A deck of cards is currently a std::vector< size_t >. A 0 is no card, a 13 is a Skip-Bo card.
 * There are 162 cards, 12 of each + 18 Skip-Bo cards.
 *
 */
class SkipBoGame
{
public:

    SkipBoGame( const size_t nplayers, const size_t stock_pile_size );

    size_t deal_from_draw_pile();

    // Only returns the number of the top card (if Skip-Bo card, the card it represents is returned).
    // i = 0, 1, 2, 3
    size_t get_build_pile( const size_t i ) const;

    void add_to_build_pile( const size_t i );

private:
    BagOfNumbers draw_pile_;
    std::vector< std::vector< size_t > > stock_piles_; // The players' decks
    std::vector< std::vector< size_t > > build_piles_; // Can't store top card only because must keep track of Skip-Bo cards
    std::vector< std::vector< size_t > > hands_;
    std::vector< std::vector< std::vector< size_t > > > discard_piles_;

    size_t nplayers_;
};

class SkipBoPlayer
{
public:

    SkipBoPlayer();

    void play( SkipBoGame & skip_bo_game );

    bool finished() const { return false; } // @@

private:
    std::vector< size_t > hand_;
    std::vector< std::vector< size_t > > discard_piles_;

};


#endif // SKIPBO_H

