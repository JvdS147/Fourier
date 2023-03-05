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

#include "ChebyshevBackground.h"

#include <cmath>
#include <stdexcept>

// ********************************************************************************

// Returns a Chebyshev polynomial of the first kind.
double Chebyshev( const size_t order, const double x )
{
    switch( order )
    {
        case  0 : return 1;
        case  1 : return x;
        case  2 : return 2.0*std::pow( x, 2 ) - 1.0;
        case  3 : return 4.0*std::pow( x, 3 ) - 3.0*x;
        case  4 : return 8.0*std::pow( x, 4 ) - 8.0*std::pow( x, 2 ) + 1.0;
        case  5 : return 16.0*std::pow( x, 5 ) - 20.0*std::pow( x, 3 ) + 5.0*x;
        case  6 : return 32.0*std::pow( x, 6 ) - 48.0*std::pow( x, 4 ) + 18.0*std::pow( x, 2 ) - 1.0;
        case  7 : return 64.0*std::pow( x, 7 ) - 112.0*std::pow( x, 5 ) + 56.0*std::pow( x, 3 ) - 7.0*x;
        case  8 : return 128.0*std::pow( x, 8 ) - 256.0*std::pow( x, 6 ) + 160.0*std::pow( x, 4 ) - 32.0*std::pow( x, 2 ) + 1.0;
        case  9 : return 256.0*std::pow( x, 9 ) - 576.0*std::pow( x, 7 ) + 432.0*std::pow( x, 5 ) - 120.0*std::pow( x, 3 ) + 9.0*x;
        case 10 : return 512.0*std::pow( x, 10 ) - 1280.0*std::pow( x, 8 ) + 1120.0*std::pow( x, 6 ) - 400.0*std::pow( x, 4 ) + 50.0*std::pow( x, 2 ) - 1.0;
        default : throw std::runtime_error( "Chebyshev(): order not yet implemented." );
    }
}

// ********************************************************************************

