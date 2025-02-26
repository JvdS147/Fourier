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

#include "CrystalLattice.h"

#include "TestSuite.h"

#include <iostream>

void test_crystal_lattice( TestSuite & test_suite )
{
    std::cout << "Now running tests for CrystalLattice." << std::endl;

    {
    CrystalLattice crystal_lattice( 1.0, 1.0, 100.0, Angle::from_degrees( 90.0 ), Angle::from_degrees( 100.0 ), Angle::from_degrees( 90.0 ) );
    test_suite.test_equality_double( crystal_lattice.beta_star().value_in_degrees(), 80.0, "beta_star()" );
    }
    {
    CrystalLattice crystal_lattice( 4.56, 10.2, 12.34, Angle::from_degrees( 89.0 ), Angle::from_degrees( 92.0 ), Angle::from_degrees( 75.6 ) );
    test_suite.test_equality( crystal_lattice.lattice_system(), CrystalLattice::TRICLINIC, "deduce_lattice_system() triclinic" );
    }
    {
    CrystalLattice crystal_lattice( 4.56, 10.2, 12.34, Angle::from_degrees( 92.0 ), Angle::angle_90_degrees(), Angle::angle_90_degrees() );
    test_suite.test_equality( crystal_lattice.lattice_system(), CrystalLattice::MONOCLINIC_A, "deduce_lattice_system() monoclinic 1" );
    }
    {
    CrystalLattice crystal_lattice( 4.56, 10.2, 12.34, Angle::angle_90_degrees(), Angle::from_degrees( 92.0 ), Angle::angle_90_degrees() );
    test_suite.test_equality( crystal_lattice.lattice_system(), CrystalLattice::MONOCLINIC_B, "deduce_lattice_system() monoclinic 2" );
    }
    {
    CrystalLattice crystal_lattice( 4.56, 10.2, 12.34, Angle::angle_90_degrees(), Angle::angle_90_degrees(), Angle::from_degrees( 92.0 ) );
    test_suite.test_equality( crystal_lattice.lattice_system(), CrystalLattice::MONOCLINIC_C, "deduce_lattice_system() monoclinic 3" );
    }
    {
    CrystalLattice crystal_lattice( 4.56, 10.2, 12.34, Angle::angle_90_degrees(), Angle::angle_90_degrees(), Angle::angle_90_degrees() );
    test_suite.test_equality( crystal_lattice.lattice_system(), CrystalLattice::ORTHORHOMBIC, "deduce_lattice_system() orthorhombic" );
    }
    {
    CrystalLattice crystal_lattice( 4.56, 4.56, 12.34, Angle::angle_90_degrees(), Angle::angle_90_degrees(), Angle::angle_90_degrees() );
    test_suite.test_equality( crystal_lattice.lattice_system(), CrystalLattice::TETRAGONAL, "deduce_lattice_system() tetragonal" );
    }
    {
    CrystalLattice crystal_lattice( 4.56, 4.56, 12.34, Angle::angle_90_degrees(), Angle::angle_90_degrees(), Angle::angle_120_degrees() );
    test_suite.test_equality( crystal_lattice.lattice_system(), CrystalLattice::HEXAGONAL, "deduce_lattice_system() hexagonal" );
    }
    {
    CrystalLattice crystal_lattice( 10.2, 10.2, 10.2, Angle::from_degrees( 92.0 ), Angle::from_degrees( 92.0 ), Angle::from_degrees( 92.0 ) );
    test_suite.test_equality( crystal_lattice.lattice_system(), CrystalLattice::RHOMBOHEDRAL, "deduce_lattice_system() rhombohedral" );
    }
    {
    CrystalLattice crystal_lattice( 10.2, 10.2, 10.2, Angle::angle_90_degrees(), Angle::angle_90_degrees(), Angle::angle_90_degrees() );
    test_suite.test_equality( crystal_lattice.lattice_system(), CrystalLattice::CUBIC, "deduce_lattice_system() cubic" );
    }
}

