#ifndef RUNTESTS_H
#define RUNTESTS_H

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

class TestSuite;

void test_angle( TestSuite & test_suite );
void test_Chebyshev_background( TestSuite & test_suite );
void test_chemical_formula( TestSuite & test_suite );
void test_Complex( TestSuite & test_suite );
void test_Constraints( TestSuite & test_suite );
void test_ConvexPolygon( TestSuite & test_suite );
void test_correlation_matrix( TestSuite & test_suite );
void test_crystal_lattice( TestSuite & test_suite );
void test_crystal_structure( TestSuite & test_suite );
void test_file_name( TestSuite & test_suite );
void test_fraction( TestSuite & test_suite );
void test_mapping( TestSuite & test_suite );
void test_matrix3D( TestSuite & test_suite );
void test_MatrixFraction3D( TestSuite & test_suite );
void test_maths( TestSuite & test_suite );
void test_ModelBuilding( TestSuite & test_suite );
void test_PowderMatchTable( TestSuite & test_suite );
void test_PowderPattern( TestSuite & test_suite );
void test_quaternion( TestSuite & test_suite );
void test_ReadCell( TestSuite & test_suite );
void test_sort( TestSuite & test_suite );
void test_TextFileReader_2( TestSuite & test_suite );
void test_TLS_ADPs( TestSuite & test_suite );
void test_utilities( TestSuite & test_suite );
void test_3D_calculations( TestSuite & test_suite );

void run_tests();

#endif // RUNTESTS_H

