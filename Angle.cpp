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

#include "Angle.h"

#include <iostream>

// The Angle class is supposed to be lightweight, so there should not that much code related to the class in here.
// But there are some helper functions in here.

// ********************************************************************************

std::ostream & operator<<( std::ostream & os, const Angle angle )
{
    os << angle.value_in_degrees();
    return os;
}

// ********************************************************************************

void sincos( Angle angle, double & sine, double & cosine )
{
//    sine = sin( x );

// Sine
    double x = angle.value_in_radians();
    while ( x > CONSTANT_PI )
        x -= 2.0*CONSTANT_PI;
    while ( x < -CONSTANT_PI )
        x += 2.0*CONSTANT_PI;
    double B = 4.0/CONSTANT_PI;
    double C = -4.0/(CONSTANT_PI*CONSTANT_PI);
    double y = B * x + C * x * fabs(x);
    double P = 0.225;
    sine = P * ( y * fabs(y) - y ) + y; // = Q * y + P * y * fabs(y), Q = 0.775 ( P + Q = 1.0 )

//    cosine = cos( x );
    x += CONSTANT_PI/2.0; // cos(x) = sin(x + pi/2)
    if ( x > CONSTANT_PI ) // Original x > pi/2
        x -= 2.0*CONSTANT_PI; // Wrap: cos(x) = cos(x - 2 pi)
    y = B * x + C * x * fabs(x);
    cosine = P * ( y * fabs(y) - y ) + y; // = Q * y + P * y * fabs(y), Q = 0.775 ( P + Q = 1.0 )
}

// ********************************************************************************

bool triquality( const Angle x1, const Angle x2, const Angle x3, const Angle tolerance )
{
    Angle average = ( x1 + x2 + x3 ) / 3.0;
    if ( ! nearly_equal( x1, average, tolerance ) )
        return false;
    if ( ! nearly_equal( x2, average, tolerance ) )
        return false;
    if ( ! nearly_equal( x3, average, tolerance ) )
        return false;
    return true;
}

// ********************************************************************************

