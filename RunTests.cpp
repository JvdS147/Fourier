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

#include "TestSuite.h"
#include "RunTests.h"

#include <iostream>

void run_tests()
{
    TestSuite test_suite;
    try
    {
        test_angle( test_suite );
        test_Chebyshev_background( test_suite );
        test_chemical_formula( test_suite );
        test_Complex( test_suite );
        test_Constraints( test_suite );
        test_ConvexPolygon( test_suite );
        test_correlation_matrix( test_suite );
        test_crystal_lattice( test_suite );
        test_crystal_structure( test_suite );
        test_file_name( test_suite );
        test_fraction( test_suite );
        test_linear_regression( test_suite );
        test_mapping( test_suite );
        test_matrix3D( test_suite );
        test_MatrixFraction3D( test_suite );
        test_maths( test_suite );
        test_ModelBuilding( test_suite );
        test_OrientationalOrderParameters( test_suite );
        test_PowderPattern( test_suite );
        test_quaternion( test_suite );
        test_ReadCell( test_suite );
        test_sort( test_suite );
        test_space_group( test_suite );
        test_SphericalHarmonics( test_suite );
        test_StringFunctions( test_suite );
        test_StringConversions( test_suite );
        test_SudokuSolver( test_suite );
        test_TextFileReader_2( test_suite );
        test_TLS_ADPs( test_suite );
        test_utilities( test_suite );
        test_3D_calculations( test_suite );
    }
    catch ( std::exception& e )
    {
        std::cout << "An exception was thrown" << std::endl;
        std::cout << e.what() << std::endl;
    }
    test_suite.report();
    std::cout << "Test suite done" << std::endl;
}

