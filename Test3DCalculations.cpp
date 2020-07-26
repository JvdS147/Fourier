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

#include "TestSuite.h"
#include "3DCalculations.h"
#include "AnisotropicDisplacementParameters.h"
#include "CrystalLattice.h"
#include "Eigenvalue.h"
#include "NormalisedVector3D.h"

#include <iostream>

void test_3D_calculations( TestSuite & test_suite )
{
    std::cout << "Now running tests for 3DCalculations." << std::endl;
    {
        SymmetricMatrix3D a = SymmetricMatrix3D( 0.01038, 0.01333, 0.01294, -0.00184, 0.00248, -0.00623 ); // KOFHOC01
        std::vector< double > eigenvalues;
        std::vector< NormalisedVector3D > eigenvectors;
        calculate_eigenvalues( a, eigenvalues, eigenvectors );
        test_suite.test_equality_double( eigenvalues[0], 0.00681342, "calculate_eigenvalues() value 0" );
        test_suite.test_equality_double( eigenvalues[1], 0.00953122, "calculate_eigenvalues() value 1" );
        test_suite.test_equality_double( eigenvalues[2], 0.0203054 , "calculate_eigenvalues() value 2" );
        test_suite.test_equality_double( eigenvectors[0].value(0), -0.174105, "calculate_eigenvalues() vector 0" );
        test_suite.test_equality_double( eigenvectors[0].value(1),  0.654352, "calculate_eigenvalues() vector 1" );
        test_suite.test_equality_double( eigenvectors[0].value(2),  0.735874, "calculate_eigenvalues() vector 2" );
    }

    {
        SymmetricMatrix3D U_cif( 0.04682, 0.02399, 0.05938, 0.00258, -0.00735, -0.00165 );

        Matrix3D transformation_matrix(  1.0,  0.0,  0.0,
                                         0.0, -1.0,  0.0,
                                        -1.0,  0.0, -1.0 );

        CrystalLattice crystal_lattice( 13.8593,
                                        10.5242,
                                        14.8044,
                                        Angle::from_degrees( 90.0 ),
                                        Angle::from_degrees( 92.0014 ),
                                        Angle::from_degrees( 90.0 ) );

        SymmetricMatrix3D U_cif_new = transform_adps( U_cif, transformation_matrix, crystal_lattice );
        test_suite.test_equality_double( U_cif_new.value( 0, 0 ),  0.06219, "transform_adps() element 0 0", 0.00001 );
        test_suite.test_equality_double( U_cif_new.value( 1, 1 ),  0.02399, "transform_adps() element 1 1", 0.00001 );
        test_suite.test_equality_double( U_cif_new.value( 2, 2 ),  0.05938, "transform_adps() element 2 2", 0.00001 );
        test_suite.test_equality_double( U_cif_new.value( 0, 1 ), -0.00306, "transform_adps() element 0 1", 0.00001 );
        test_suite.test_equality_double( U_cif_new.value( 0, 2 ),  0.04677, "transform_adps() element 0 2", 0.00001 );
        test_suite.test_equality_double( U_cif_new.value( 1, 2 ), -0.00165, "transform_adps() element 1 2", 0.00001 );
    }

}

