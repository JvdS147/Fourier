#ifndef MILLERINDICES_H
#define MILLERINDICES_H

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

#include <iosfwd>

/*
  Miller indices. True Miller indices are relative prime (co-prime), if they are not then technically they are
  called Laue indices. These Miller indices are NOT relative prime, so 002 is a perfectly acceptable
  value for this class. However, there is a member function to make all indices relative prime.

  We mainly need this class for sorting.
*/
class MillerIndices
{
public:

    MillerIndices( const int h, const int k, const int l );

    int h() const { return h_; }
    int k() const { return k_; }
    int l() const { return l_; }
    
    void make_relative_prime();

    // if (hkl) = (000), there can be surprises.
    bool is_000() const;

    MillerIndices operator-() const
    {
        return MillerIndices( -this->h(), -this->k(), -this->k() );
    }

    std::string to_string() const;

private:
    int h_;
    int k_;
    int l_;
};

std::ostream & operator<<( std::ostream & os, const MillerIndices & miller_indices );

bool operator<( const MillerIndices & lhs, const MillerIndices & rhs );
bool operator==( const MillerIndices & lhs, const MillerIndices & rhs );
// The transposition is implied
int operator*( const MillerIndices & lhs, const MillerIndices & rhs );

#endif // MILLERINDICES_H

