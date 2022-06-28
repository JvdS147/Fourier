/* *********************************************
Copyright (c) 2013-2022, Cornelis Jan (Jacco) van de Streek
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

#include "CrystalStructures.h"

#include "CrystalStructure.h"

// ********************************************************************************

CrystalStructure NaCl()
{
    CrystalStructure result;
    result.set_name( "NaCl" );
    result.add_atom( Atom( Element( "Na" ), Vector3D( 0.0, 0.0, 0.0 ), "Na1", 1.0, AnisotropicDisplacementParameters( 0.01 ) ) );
    result.add_atom( Atom( Element( "Cl" ), Vector3D( 0.0, 0.5, 0.0 ), "Cl1", -1.0, AnisotropicDisplacementParameters( 0.01 ) ) );
    result.set_crystal_lattice( CrystalLattice( 5.6393, 5.6393, 5.6393, Angle::angle_90_degrees(), Angle::angle_90_degrees(), Angle::angle_90_degrees() ) );
    std::vector< SymmetryOperator > symmetry_operators;
    symmetry_operators.push_back( SymmetryOperator( "x,y,z" ) );
    symmetry_operators.push_back( SymmetryOperator( "z,x,y" ) );
    symmetry_operators.push_back( SymmetryOperator( "y,z,x" ) );
    symmetry_operators.push_back( SymmetryOperator( "-x,-y,z" ) );
    symmetry_operators.push_back( SymmetryOperator( "z,-x,-y" ) );
    symmetry_operators.push_back( SymmetryOperator( "-y,z,-x" ) );
    symmetry_operators.push_back( SymmetryOperator( "-z,x,-y" ) );
    symmetry_operators.push_back( SymmetryOperator( "-y,-z,x" ) );
    symmetry_operators.push_back( SymmetryOperator( "y,-z,-x" ) );
    symmetry_operators.push_back( SymmetryOperator( "-x,y,-z" ) );
    symmetry_operators.push_back( SymmetryOperator( "-z,-x,y" ) );
    symmetry_operators.push_back( SymmetryOperator( "x,-y,-z" ) );
    symmetry_operators.push_back( SymmetryOperator( "-y,-x,-z" ) );
    symmetry_operators.push_back( SymmetryOperator( "-z,-y,-x" ) );
    symmetry_operators.push_back( SymmetryOperator( "-x,-z,-y" ) );
    symmetry_operators.push_back( SymmetryOperator( "y,x,-z" ) );
    symmetry_operators.push_back( SymmetryOperator( "-z,y,x" ) );
    symmetry_operators.push_back( SymmetryOperator( "x,-z,y" ) );
    symmetry_operators.push_back( SymmetryOperator( "z,-y,x" ) );
    symmetry_operators.push_back( SymmetryOperator( "x,z,-y" ) );
    symmetry_operators.push_back( SymmetryOperator( "-x,z,y" ) );
    symmetry_operators.push_back( SymmetryOperator( "y,-x,z" ) );
    symmetry_operators.push_back( SymmetryOperator( "z,y,-x" ) );
    symmetry_operators.push_back( SymmetryOperator( "-y,x,z" ) );
    SpaceGroup space_group( symmetry_operators );
    space_group.add_inversion_at_origin();
    space_group.add_centring( Centring( "F" ) );
    result.set_space_group( space_group );
    return result;
}

// ********************************************************************************


