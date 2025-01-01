/* *********************************************
Copyright (c) 2013-2025, Cornelis Jan (Jacco) van de Streek
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

#include "SphericalHarmonics.h"

#include "TestSuite.h"

#include <iostream>

void test_SphericalHarmonics( TestSuite & test_suite )
{
    std::cout << "Now running tests for SphericalHarmonics." << std::endl;
    {
        Angle alpha = Angle::from_degrees( 45.0 );
        Angle beta = Angle::from_degrees( 90.0 );
        Complex result;
        Complex target;
        result = spherical_harmonics( 2, -2, alpha, beta );
        target = 0.25 * sqrt(15.0/(2.0*CONSTANT_PI))*square(beta.sine())*exponential( Complex( 0.0, -2.0*alpha.value_in_radians() ) );
        test_suite.test_equality_Complex( result, target, "SphericalHarmonics() 01" );
        result = spherical_harmonics( 2, -1, alpha, beta );
        target = sqrt(15.0/(8.0*CONSTANT_PI))*beta.sine()*beta.cosine()*exponential( Complex( 0.0, -alpha.value_in_radians() ) );
        test_suite.test_equality_Complex( result, target, "SphericalHarmonics() 02" );        
        result = spherical_harmonics( 2,  0, alpha, beta );
        target = Complex( sqrt(5.0/(4.0*CONSTANT_PI))*( (3.0/2.0)*square(beta.cosine()) - 0.5 ) );
        test_suite.test_equality_Complex( result, target, "SphericalHarmonics() 03" );
        result = spherical_harmonics( 2,  1, alpha, beta );
        target = -sqrt(15.0/(8.0*CONSTANT_PI))*beta.sine()*beta.cosine()*exponential( Complex( 0.0, alpha.value_in_radians() ) );
        test_suite.test_equality_Complex( result, target, "SphericalHarmonics() 04" );
        result = spherical_harmonics( 2,  2, alpha, beta );
        target = 0.25 * sqrt(15.0/(2.0*CONSTANT_PI))*square(beta.sine())*exponential( Complex( 0.0, 2.0*alpha.value_in_radians() ) );
        test_suite.test_equality_Complex( result, target, "SphericalHarmonics() 05" );
    }

}

