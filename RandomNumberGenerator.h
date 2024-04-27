#ifndef RANDOMNUMBERGENERATOR_H
#define RANDOMNUMBERGENERATOR_H

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

#include <vector>

// The classes for integers and for doubles are virtually identical, the decision to split them
// was deliberate: if both are combined into one class, the integer random number generator is used
// for both the integers and the doubles and this gives a correlation between the two. The way it is
// currently implemented, the two are entirely independent, yet both give reproducible results.

class RandomNumberGenerator_integer
{
public:
    
    // idum should be positive and not zero
    RandomNumberGenerator_integer( const int idum = 1539 );

    int next_number( const int start, const int end );

private:
    int idum_;
    int idum2_;
    int iy_;
    std::vector< int > iv_;
    
    int IM1;
    int IM2;
    int IA1;
    int IA2;
    int IQ1;
    int IQ2;
    int IR1;
    int IR2;
    int NTAB;
    int IMM1;
    int NDIV;
};

// Returns a random number in the range [0.0, 1.0].
class RandomNumberGenerator_double
{
public:
    
    // idum should be positive and not zero
    explicit RandomNumberGenerator_double( const int idum = 7381 );

    double next_number();

private:
    int idum_;
    int idum2_;
    int iy_;
    std::vector< int > iv_;
    
    int IM1;
    int IM2;
    int IA1;
    int IA2;
    int IQ1;
    int IQ2;
    int IR1;
    int IR2;
    int NTAB;
    int IMM1;
    int NDIV;
    double AM;
};

#endif // RANDOMNUMBERGENERATOR_H

