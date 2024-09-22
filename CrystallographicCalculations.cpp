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

#include "CrystallographicCalculations.h"
#include "3DCalculations.h"
#include "Centring.h"
#include "CrystalStructure.h"
#include "Matrix3D.h"
#include "MillerIndices.h"
#include "NormalisedVector3D.h"
#include "PointGroup.h"
#include "SpaceGroup.h"
#include "SymmetryOperator.h"
#include "Utilities.h"
#include "Vector3D.h"

#include <stdexcept>


// ********************************************************************************

std::vector< Matrix3D > orthorhombic_unit_cell_axes_permutations()
{
    std::vector< Matrix3D > result;
    result.push_back( Matrix3D(  1.0,  0.0,  0.0,
                                 0.0,  1.0,  0.0,
                                 0.0,  0.0,  1.0 ) );
    result.push_back( Matrix3D(  0.0,  0.0,  1.0,
                                 1.0,  0.0,  0.0,
                                 0.0,  1.0,  0.0 ) );
    result.push_back( Matrix3D(  0.0,  1.0,  0.0,
                                 0.0,  0.0,  1.0,
                                 1.0,  0.0,  0.0 ) );
    result.push_back( Matrix3D(  0.0,  1.0,  0.0,
                                 1.0,  0.0,  0.0,
                                 0.0,  0.0, -1.0 ) );
    result.push_back( Matrix3D(  0.0,  0.0,  1.0,
                                 0.0, -1.0,  0.0,
                                 1.0,  0.0,  0.0 ) );
    result.push_back( Matrix3D( -1.0,  0.0,  0.0,
                                 0.0,  0.0,  1.0,
                                 0.0,  1.0,  0.0 ) );
    return result;
}

// ********************************************************************************

std::vector< std::string > orthorhombic_unit_cell_axes_permutation_labels()
{
    std::vector< std::string > result;
    result.push_back( "abc" );
    result.push_back( "cab" );
    result.push_back( "bca" );
    result.push_back( "ba-c" );
    result.push_back( "c-ba" );
    result.push_back( "-acb" );
    return result;
}

// ********************************************************************************

Vector3D reciprocal_lattice_point( const MillerIndices miller_indices, const CrystalLattice & crystal_lattice )
{
    return ( miller_indices.h() * crystal_lattice.a_star_vector() +
             miller_indices.k() * crystal_lattice.b_star_vector() +
             miller_indices.l() * crystal_lattice.c_star_vector() );
}

// ********************************************************************************

NormalisedVector3D reciprocal_lattice_direction( const MillerIndices miller_indices, const CrystalLattice & crystal_lattice )
{
    return normalised_vector( miller_indices.h() * crystal_lattice.a_star_vector() +
                              miller_indices.k() * crystal_lattice.b_star_vector() +
                              miller_indices.l() * crystal_lattice.c_star_vector() );
}

// ********************************************************************************

//MillerIndices operator*( const Matrix3D & matrix, const MillerIndices & Miller_indices )
//{
//    return MillerIndices( round_to_int( matrix.value( 0, 0 ) * Miller_indices.h() + matrix.value( 0, 1 ) * Miller_indices.k() + matrix.value( 0, 2 ) * Miller_indices.l() ),
//                          round_to_int( matrix.value( 1, 0 ) * Miller_indices.h() + matrix.value( 1, 1 ) * Miller_indices.k() + matrix.value( 1, 2 ) * Miller_indices.l() ),
//                          round_to_int( matrix.value( 2, 0 ) * Miller_indices.h() + matrix.value( 2, 1 ) * Miller_indices.k() + matrix.value( 2, 2 ) * Miller_indices.l() )
//                        );
//}

// ********************************************************************************

MillerIndices operator*( const MillerIndices & Miller_indices, const Matrix3D & matrix )
{
    return MillerIndices( round_to_int( Miller_indices.h() * matrix.value( 0, 0 ) + Miller_indices.k() * matrix.value( 1, 0 ) + Miller_indices.l() * matrix.value( 2, 0 ) ),
                          round_to_int( Miller_indices.h() * matrix.value( 0, 1 ) + Miller_indices.k() * matrix.value( 1, 1 ) + Miller_indices.l() * matrix.value( 2, 1 ) ),
                          round_to_int( Miller_indices.h() * matrix.value( 0, 2 ) + Miller_indices.k() * matrix.value( 1, 2 ) + Miller_indices.l() * matrix.value( 2, 2 ) )
                        );
}

// ********************************************************************************

MillerIndices operator*( const MillerIndices & Miller_indices, const SymmetricMatrix3D & matrix )
{
    return MillerIndices( round_to_int( Miller_indices.h() * matrix.value( 0, 0 ) + Miller_indices.k() * matrix.value( 1, 0 ) + Miller_indices.l() * matrix.value( 2, 0 ) ),
                          round_to_int( Miller_indices.h() * matrix.value( 0, 1 ) + Miller_indices.k() * matrix.value( 1, 1 ) + Miller_indices.l() * matrix.value( 2, 1 ) ),
                          round_to_int( Miller_indices.h() * matrix.value( 0, 2 ) + Miller_indices.k() * matrix.value( 1, 2 ) + Miller_indices.l() * matrix.value( 2, 2 ) )
                        );
}

// ********************************************************************************

double operator*( const MillerIndices & miller_indices, const Vector3D & vector_3D )
{
    return miller_indices.h() * vector_3D.x() + miller_indices.k() * vector_3D.y() + miller_indices.l() * vector_3D.z();
}

// ********************************************************************************

void add_centring_to_space_group_after_transformation( Matrix3D tranformation_matrix, SpaceGroup & space_group )
{
    double d = tranformation_matrix.determinant();
    if ( ! nearly_integer( d ) )
        throw std::runtime_error( "add_centring_to_space_group_after_transformation() : determinant is not an integer." );
    if ( nearly_zero( d ) )
        throw std::runtime_error( "add_centring_to_space_group_after_transformation() : determinant is zero." );
    if ( d < 0.0 )
        throw std::runtime_error( "add_centring_to_space_group_after_transformation() : determinant is negative." );
    size_t D = round_to_int( tranformation_matrix.determinant() );
    if ( D == 1 )
        return;
    // I have not been able to find a smart way to extract the possible additional lattice points
    // from the transformation matrix, so we simply try them all.
    // We want to find the points [ f/D, g/D, h/D ], with D the determinant, that lie within the unit cell
    // but that is not one of the current lattice points. So f/D, g/D and h/D are not allowed all to be integers at once.
    // f/D must be in the range [ 0, 1 >. 
    tranformation_matrix.transpose();
    // tranformation_matrix /= d;
    std::vector< Vector3D > centring_vectors;
    for ( size_t f( 0 ); f != D; ++f )
    {
        for ( size_t g( 0 ); g != D; ++g )
        {
            for ( size_t h( 0 ); h != D; ++h )
            {
                Vector3D trial_vector( f/d, g/d, h/d );
                // I guess it would be more efficient to divide the transformation matrix by d and
                // to construct the trial vector as [ f, g, h ].
                Vector3D lp = tranformation_matrix * trial_vector; // lp = (trial) lattice point in the old coordinate frame.
                if ( ! nearly_integer( lp.x() ) )
                    continue;
                if ( ! nearly_integer( lp.y() ) )
                    continue;
                if ( ! nearly_integer( lp.z() ) )
                    continue;
                centring_vectors.push_back( trial_vector );
            }
        }
    }
    if ( centring_vectors.size() != D )
        throw std::runtime_error( "add_centring_to_space_group_after_transformation() : centring_vectors.size() != D." );
    space_group.add_centring( Centring( centring_vectors ) );
}

// ********************************************************************************

std::vector< SymmetryOperator > centring_generators( const Centring & centring )
{
    if ( centring.centring_type() == Centring::U )
        throw std::runtime_error( "centring_generators() : error: unknown centring." );
    std::vector< SymmetryOperator > result;
    if ( centring.centring_type() == Centring::P )
        return result;
    Matrix3D identity;
    if ( centring.centring_type() == Centring::J )
    {
        result.push_back( SymmetryOperator( identity, Vector3D( 0.1, 0.0, 0.0 ) ) );
        result.push_back( SymmetryOperator( identity, Vector3D( 0.0, 0.1, 0.0 ) ) );
        result.push_back( SymmetryOperator( identity, Vector3D( 0.0, 0.0, 0.1 ) ) );
        return result;
    }
    if ( centring.centring_type() == Centring::F )
    {
        result.push_back( SymmetryOperator( identity, Vector3D( 0.0, 0.5, 0.5 ) ) );
        result.push_back( SymmetryOperator( identity, Vector3D( 0.5, 0.0, 0.5 ) ) );
        return result;
    }
    std::vector< Vector3D > centring_vectors = centring.centring_vectors();
    for ( size_t i( 1 ); i != centring_vectors.size(); ++i )
        result.push_back( SymmetryOperator( identity, centring_vectors[i] ) );
    return result;
}

// ********************************************************************************

Centring expand_centring_generators( const std::vector< SymmetryOperator > & generators )
{
    // @@ Currently we can only cope with whatever was generated by std::vector< SymmetryOperator > centring_generators( const Centring & centring ).
    if ( generators.size() == 0 )
        return Centring( "P" );
    if ( generators.size() == 3 )
        return Centring( "J" );
    if ( ( generators.size() == 2 ) && nearly_contains( generators, SymmetryOperator( Matrix3D(), Vector3D( 0.0, 0.5, 0.5 ) ) ) && nearly_contains( generators, SymmetryOperator( Matrix3D(), Vector3D( 0.5, 0.0, 0.5 ) ) ) )
        return Centring( "F" );
    std::vector< Vector3D > centring_vectors;
    centring_vectors.push_back( Vector3D() );
    for ( size_t i( 0 ); i != generators.size(); ++i )
        centring_vectors.push_back( generators[i].translation() );
    return Centring( centring_vectors );
}

// ********************************************************************************

bool nearly_equal( const CrystalStructure & lhs, const CrystalStructure & rhs )
{
    return true;
}

// ********************************************************************************

AnisotropicDisplacementParameters adjust_to_site_symmetry( const AnisotropicDisplacementParameters & adps, const PointGroup & point_group, const CrystalLattice & crystal_lattice )
{
    // We need U_star.
    SymmetricMatrix3D U_star = adps.U_star( crystal_lattice );
    SymmetricMatrix3D sum( U_star );
    for ( size_t i( 1 ); i != point_group.nsymmetry_operators(); ++i )
        sum += Matrix3D2SymmetricMatrix3D( point_group.symmetry_operator( i ) * U_star * transpose( point_group.symmetry_operator( i ) ) );
    U_star = sum / point_group.nsymmetry_operators();
    SymmetricMatrix3D U_cart = U_star_2_U_cart( U_star, crystal_lattice );
    // Some values are now e.g. 2.23E-34, round them to 0.0, but not for the diagonal.
//    for ( size_t i( 0 ); i != 3; ++i )
//    {
//        for ( size_t j( i+1 ); j != 3; ++j )
//        {
//            if ( nearly_zero( U_cart.value( i, j ) ) )
//                U_cart.set_value( i, j, 0.0 );
//        }
//    }
    return AnisotropicDisplacementParameters( U_cart );
}

// ********************************************************************************

MillerIndices select_realistic_preferred_orientation_direction( const CrystalLattice & crystal_lattice )
{
    switch ( crystal_lattice.lattice_system() )
    {
        case CrystalLattice::TRICLINIC    :
        case CrystalLattice::ORTHORHOMBIC : {
                                                if ( crystal_lattice.a() < crystal_lattice.b() )
                                                {
                                                    if ( crystal_lattice.a() < crystal_lattice.c() )
                                                        return MillerIndices( 1, 0, 0 );
                                                }
                                                else // b < a
                                                {
                                                    if ( crystal_lattice.b() < crystal_lattice.c() )
                                                        return MillerIndices( 0, 1, 0 );
                                                }
                                                return MillerIndices( 0, 0, 1 );
                                            }
        case CrystalLattice::MONOCLINIC_A : return MillerIndices( 1, 0, 0 );
        case CrystalLattice::MONOCLINIC_B : return MillerIndices( 0, 1, 0 );
        case CrystalLattice::MONOCLINIC_C :
        case CrystalLattice::TETRAGONAL   :
        case CrystalLattice::HEXAGONAL    : return MillerIndices( 0, 0, 1 );
        case CrystalLattice::RHOMBOHEDRAL :
        case CrystalLattice::CUBIC        : throw std::runtime_error( "select_realistic_preferred_orientation_direction() : error: no PO possible." );
    }
}

// ********************************************************************************

MillerIndices select_realistic_preferred_orientation_direction( const CrystalLattice & crystal_lattice, const SpaceGroup & space_group )
{
    MillerIndices result = select_realistic_preferred_orientation_direction( crystal_lattice );
//    // Check that the PO direction is commensurate with the space-group symmetry.
    MillerIndices reflection( 37, -23, 3 );
    Vector3D PO_vector = reciprocal_lattice_point( result, crystal_lattice );
    Vector3D H = reciprocal_lattice_point( reflection, crystal_lattice );
    double reference_dot_product = absolute( PO_vector * H );
//    for ( size_t i( 0 ); i != Laue_class_.nsymmetry_operators(); ++i )
//    {
//        MillerIndices equivalent_reflection = reflection * Laue_class_.symmetry_operator( i );
//        // Now check that the March-Dollase PO corrections are the same for all of them.
//        Vector3D H = reciprocal_lattice_point( equivalent_reflection, crystal_lattice);
//        double current_dot_product = absolute( PO_vector * H );
//        if ( ! nearly_equal( current_dot_product, reference_dot_product ) )
//        {
//            std::cout << "select_realistic_preferred_orientation_direction(): Warning: PO direction is not commensurate with space-group symmetry." << std::endl;
//            return;
//        }
//    }
    return result;
}

// ********************************************************************************

