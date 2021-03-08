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

#include "TLS.h"

#include "TestSuite.h"

#include <iostream>

void test_TLS( TestSuite & test_suite )
{
    std::cout << "Now running tests for TLS." << std::endl;

    {
        double from_nm = 10.0;
        double from_rad = 1.0;

        double T11 = from_nm * from_nm * 0.00024;
        double T12 = from_nm * from_nm * 0.0;
        double T13 = from_nm * from_nm * -0.00026;
        double T23 = T13; // Why ??? Because it is on a centre of symmetry???
        double T22 = from_nm * from_nm * 0.00035;
        double T33 = from_nm * from_nm * -0.0005;

        SymmetricMatrix3D T( T11, T22, T33, T12, T13, T23 );
        
        double L11 = from_rad * 0.05;
        double L12 = from_rad * 0.0;
        double L13 = from_rad * -0.002;
        double L23 = L13; // Why ??? Because it is on a centre of symmetry???
        double L22 = from_rad * 0.06;
        double L33 = from_rad * 0.004;

        SymmetricMatrix3D L( L11, L22, L33, L12, L13, L23 );

    // The following is not possible: S33 = -S11 - S22, so if S11 is 0.0, S33 and S22 must be of opposite sign.
    // I think that they chose as constraint either that S11 is 0.0 or that S22 = S33

        double S11 = from_nm * from_rad * 0.0;
        double S12 = from_nm * from_rad * 0.0;
        double S13 = from_nm * from_rad * 0.0;
        double S21 = S12;
        double S22 = from_nm * from_rad * 0.0002;
        double S23 = from_nm * from_rad * 0.0;
        double S31 = from_nm * from_rad * 0.0;
        double S32 = S23;
        double S33 = from_nm * from_rad * 0.0002;

        Matrix3D S( S11, S12, S13,
                    S21, S22, S23,
                    S31, S32, S33 );
        
        TLS tls( T, L, S );

        // test_suite.test_equality( tls, , "TLS()" );
    }

}

