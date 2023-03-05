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

#include "CrystalStructure.h"

#include "TestSuite.h"

#include <iostream>

void test_crystal_structure( TestSuite & test_suite )
{
    std::cout << "Now running tests for CrystalStructure." << std::endl;
    {
    CrystalStructure crystal_structure_1;
    CrystalStructure crystal_structure_2;
    CrystalLattice crystal_lattice_1( 10.0, 10.0, 10.0, Angle::angle_90_degrees(), Angle::angle_90_degrees(), Angle::angle_90_degrees() );
    CrystalLattice crystal_lattice_2( 10.0, 10.0, 10.0, Angle::angle_90_degrees(), Angle::angle_90_degrees(), Angle::angle_90_degrees() );
    crystal_structure_1.set_crystal_lattice( crystal_lattice_1 );
    crystal_structure_2.set_crystal_lattice( crystal_lattice_2 );
    Atom atom_1( Element( "C" ), Vector3D( 0.1, 0.0, 0.0 ), "C1" );
    crystal_structure_1.add_atom( atom_1 );
    Atom atom_2( Element( "C" ), Vector3D( 0.2, 0.0, 0.0 ), "C2" );
    crystal_structure_2.add_atom( atom_2 );
    double result = root_mean_square_Cartesian_displacement( crystal_structure_1, crystal_structure_2, false );
    test_suite.test_equality_double( result, 1.0, "CrystalStructure::root_mean_square_Cartesian_displacement() : maximum error has changed." );
    }
    {
    CrystalStructure crystal_structure;
    CrystalLattice crystal_lattice( 10.0, 10.0, 10.0, Angle::angle_90_degrees(), Angle::angle_90_degrees(), Angle::angle_90_degrees() );
    crystal_structure.set_crystal_lattice( crystal_lattice );
    SpaceGroup space_group;
    space_group.add_inversion_at_origin();
    crystal_structure.set_space_group( space_group );
    Atom atom_1( Element( "C" ), Vector3D( 0.1, 0.2, 0.3 ), "C1" );
    crystal_structure.add_atom( atom_1 );
    Atom atom_2( Element( "C" ), Vector3D( -0.1, 0.8, -1.3 ), "C2" );
    crystal_structure.add_atom( atom_2 );
    test_suite.test_equality( crystal_structure.natoms(), static_cast<size_t>(2), "CrystalStructure::reduce_to_asymmetric_unit() : Error 1." );
    crystal_structure.reduce_to_asymmetric_unit();
    test_suite.test_equality( crystal_structure.natoms(), static_cast<size_t>(1), "CrystalStructure::reduce_to_asymmetric_unit() : Error 2." );
    }
    {
    CrystalStructure crystal_structure;
    CrystalLattice crystal_lattice( 13.173, 51.417, 6.4553, Angle::angle_90_degrees(), Angle::angle_90_degrees(), Angle::angle_90_degrees() );
    crystal_structure.set_crystal_lattice( crystal_lattice );
    SpaceGroup space_group;
    Centring centring( "F" );
    space_group.add_centring( centring );
    crystal_structure.set_space_group( space_group );
    crystal_structure.reduce_to_primitive();
    test_suite.test_equality_double( crystal_structure.crystal_lattice().orthogonality_defect(), 4.6142, "CrystalStructure::choose_angles_close_to_90() : Error 1.", 0.0001 );
    Matrix3D transformation_matrix = crystal_structure.crystal_lattice().choose_angles_close_to_90();
    crystal_structure.transform( transformation_matrix );
    test_suite.test_equality_double( crystal_structure.crystal_lattice().orthogonality_defect(), 1.12236, "CrystalStructure::choose_angles_close_to_90() : Error 2.", 0.0001 );
    }

}

