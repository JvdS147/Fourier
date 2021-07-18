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

#include "TestSuite.h"

#include "3DCalculations.h"
#include "AnisotropicDisplacementParameters.h"
#include "CrystalLattice.h"
#include "BasicMathsFunctions.h"

#include <iostream>

void test_TLS_ADPs( TestSuite & test_suite )
{
    std::cout << "Now running tests for TLS and ADPs." << std::endl;

    // R. T. Downs tests. Example page 10.
    {
        CrystalLattice crystal_lattice(  4.9137,
                                         4.9137,
                                         5.4047,
                                        Angle::from_degrees( 90.00000 ),
                                        Angle::from_degrees( 90.00000 ),
                                        Angle::from_degrees( 120.00000 ) );
        AnisotropicDisplacementParameters ADPs_C1( U_star_2_U_cart( SymmetricMatrix3D( 0.0179, 0.0130, 0.0085, 0.0102, -0.0026, -0.0041 ) / ( 2.0 * square( CONSTANT_PI ) ), crystal_lattice ) );
        Vector3D v( 0.0564, -0.2672, -0.1188 );
        double mu2_v = ( v * ( transpose( crystal_lattice.Downs_G() ) * ( 2.0 * square( CONSTANT_PI ) * ADPs_C1.U_star( crystal_lattice ) ) * crystal_lattice.Downs_G() ) * v ) / ( 2.0 * square( CONSTANT_PI ) * v * crystal_lattice.Downs_G() * v );
        test_suite.test_equality_double( mu2_v, 0.00693478, "mu2_v" );
    }
    {
        CrystalLattice crystal_lattice(  4.9137,
                                         4.9137,
                                         5.4047,
                                        Angle::from_degrees( 90.00000 ),
                                        Angle::from_degrees( 90.00000 ),
                                        Angle::from_degrees( 120.00000 ) );
        AnisotropicDisplacementParameters ADPs_C1( U_star_2_U_cart( SymmetricMatrix3D( 0.0179, 0.0130, 0.0085, 0.0102, -0.0026, -0.0041 ) / ( 2.0 * square( CONSTANT_PI ) ), crystal_lattice ) );
        Vector3D v( -0.0564, 0.2672, 0.1188 );
        double mu2_v = ( v * ( transpose( crystal_lattice.Downs_G() ) * ( 2.0 * square( CONSTANT_PI ) * ADPs_C1.U_star( crystal_lattice ) ) * crystal_lattice.Downs_G() ) * v ) / ( 2.0 * square( CONSTANT_PI ) * v * crystal_lattice.Downs_G() * v );
        test_suite.test_equality_double( mu2_v, 0.00693478, "mu2_v" );
    }

}

