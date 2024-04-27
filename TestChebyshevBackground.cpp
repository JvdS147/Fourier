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

#include "ChebyshevBackground.h"
#include "PowderPattern.h"
#include "Utilities.h"

#include "TestSuite.h"

#include <iostream>

void test_Chebyshev_background( TestSuite & test_suite )
{
    std::cout << "Now running tests for ChebyshevBackground." << std::endl;

    PowderPattern powder_pattern; // Dummy input pattern. Only the 2theta values are used.
    powder_pattern.push_back( Angle::from_degrees(  2.00000 ), 20.0 );
    powder_pattern.push_back( Angle::from_degrees(  2.01638 ), 20.0 );
    powder_pattern.push_back( Angle::from_degrees(  2.03276 ), 20.0 );
    powder_pattern.push_back( Angle::from_degrees(  2.04914 ), 20.0 );
    powder_pattern.push_back( Angle::from_degrees(  2.06552 ), 20.0 );
    powder_pattern.push_back( Angle::from_degrees(  2.08190 ), 20.0 );
    powder_pattern.push_back( Angle::from_degrees(  2.09828 ), 20.0 );
    powder_pattern.push_back( Angle::from_degrees(  2.11466 ), 20.0 );
    powder_pattern.push_back( Angle::from_degrees( 40.7379 ), 20.0 );
    powder_pattern.push_back( Angle::from_degrees( 40.7543 ), 20.0 );
    powder_pattern.push_back( Angle::from_degrees( 40.7707 ), 20.0 );
    powder_pattern.push_back( Angle::from_degrees( 40.7871 ), 20.0 );
    powder_pattern.push_back( Angle::from_degrees( 40.8034 ), 20.0 );
    powder_pattern.push_back( Angle::from_degrees( 40.8198 ), 20.0 );
    powder_pattern.push_back( Angle::from_degrees( 40.8362 ), 20.0 );
    powder_pattern.push_back( Angle::from_degrees( 40.8526 ), 20.0 );
    powder_pattern.push_back( Angle::from_degrees( 40.8690 ), 20.0 );
    powder_pattern.push_back( Angle::from_degrees( 40.8853 ), 20.0 );
    powder_pattern.push_back( Angle::from_degrees( 40.9017 ), 20.0 );

    PowderPattern back_ground_calculated; // The calculated background with coefficients from a TOPAS refinement
    for ( size_t i( 0 ); i != powder_pattern.size(); ++i )
    {
        double yb = 0.0;
        std::vector< double > B; // The refineable parameters.
        B.push_back( 173.697413 );    //  0
        B.push_back(  -3.44529792 );  //  1
        B.push_back(  12.3258236 );   //  2
        B.push_back( -18.1777706 );   //  3
        B.push_back(  29.4247509 );   //  4
        B.push_back( -20.9786348 );   //  5
        B.push_back(  16.7415553 );   //  6
        B.push_back( -14.0710459 );   //  7
        B.push_back(  10.4740022 );   //  8
        B.push_back(  -4.10834492 );  //  9
        B.push_back(   4.86598724 );  // 10
        // The factor "2.0" and the term "-1.0" are there to "shift" the input from the interval [0,1] to [-1,1].
        double x = ( ( 2.0 * ( powder_pattern.two_theta( i ) - powder_pattern.two_theta_start() ) ) / ( powder_pattern.two_theta_end() - powder_pattern.two_theta_start() ) ) - 1.0;
        for ( size_t order( 0 ); order != 11; ++order )
        {
            yb += B[order] * Chebyshev( order, x );
        }
        back_ground_calculated.push_back( powder_pattern.two_theta( i ), yb );
    }
        
    PowderPattern back_ground_target; // The known correct values
    back_ground_target.push_back( Angle::from_degrees(  2.00000 ), 308.311 );
    back_ground_target.push_back( Angle::from_degrees(  2.01638 ), 304.970 );
    back_ground_target.push_back( Angle::from_degrees(  2.03276 ), 301.675 );
    back_ground_target.push_back( Angle::from_degrees(  2.04914 ), 298.424 );
    back_ground_target.push_back( Angle::from_degrees(  2.06552 ), 295.217 );
    back_ground_target.push_back( Angle::from_degrees(  2.08190 ), 292.054 );
    back_ground_target.push_back( Angle::from_degrees(  2.09828 ), 288.934 );
    back_ground_target.push_back( Angle::from_degrees(  2.11466 ), 285.857 );
    back_ground_target.push_back( Angle::from_degrees( 40.7379 ), 182.489 );
    back_ground_target.push_back( Angle::from_degrees( 40.7543 ), 182.872 );
    back_ground_target.push_back( Angle::from_degrees( 40.7707 ), 183.264 );
    back_ground_target.push_back( Angle::from_degrees( 40.7871 ), 183.666 );
    back_ground_target.push_back( Angle::from_degrees( 40.8034 ), 184.074 );
    back_ground_target.push_back( Angle::from_degrees( 40.8198 ), 184.495 );
    back_ground_target.push_back( Angle::from_degrees( 40.8362 ), 184.926 );
    back_ground_target.push_back( Angle::from_degrees( 40.8526 ), 185.367 );
    back_ground_target.push_back( Angle::from_degrees( 40.8690 ), 185.818 );
    back_ground_target.push_back( Angle::from_degrees( 40.8853 ), 186.277 );
    back_ground_target.push_back( Angle::from_degrees( 40.9017 ), 186.748 );
    
    for ( size_t i( 0 ); i != powder_pattern.size(); ++i )
    {
        test_suite.test_equality_double( back_ground_calculated.intensity( i ), back_ground_target.intensity( i ), "ChebyshevBackground " + size_t2string(i), 0.001 );
    }

}

