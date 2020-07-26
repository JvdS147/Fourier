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

#include "RandomQuaternionGenerator.h"
#include "MathConstants.h"
#include "Quaternion.h"

#include <cmath>

// ********************************************************************************

RandomQuaternionGenerator::RandomQuaternionGenerator( const int idum ):RNG_double_(idum)
{
}

// ********************************************************************************

Quaternion RandomQuaternionGenerator::next_quaternion()
{
    double u1 = RNG_double_.next_number();
    double u2 = RNG_double_.next_number();
    double u3 = RNG_double_.next_number();
    // @@ We should check that it does not have a length near 0.0
    return Quaternion( sqrt(1.0-u1)*sin(2.0*CONSTANT_PI*u2), sqrt(1.0-u1)*cos(2.0*CONSTANT_PI*u2), sqrt(u1)*sin(2.0*CONSTANT_PI*u3), sqrt(u1)*cos(2.0*CONSTANT_PI*u3) );
}

// ********************************************************************************

