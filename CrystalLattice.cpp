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

#include "CrystalLattice.h"
#include "3DCalculations.h"
#include "BasicMathsFunctions.h"
#include "Sort.h"
#include "Utilities.h"

#include <iostream>
#include <stdexcept>

// ********************************************************************************

CrystalLattice::CrystalLattice()
{
    *this = CrystalLattice( 1.0, 1.0, 1.0, Angle::angle_90_degrees(), Angle::angle_90_degrees(), Angle::angle_90_degrees() );
}

// ********************************************************************************

CrystalLattice::CrystalLattice( const double a,
                                const double b,
                                const double c,
                                const Angle alpha,
                                const Angle beta,
                                const Angle gamma ):
a_(a), b_(b), c_(c),
alpha_(alpha), beta_(beta), gamma_(gamma)
{
    double x;
    double y;
    double z;
    x = a;
    y = 0.0;
    z = 0.0;
    a_vector_ = Vector3D( x, y, z );
    x = b * gamma.cosine();
    y = b * gamma.sine();
    z = 0.0;
    b_vector_ = Vector3D( x, y, z );
    x = c * beta.cosine();
    y = ( b*c*alpha.cosine() - b_vector_.x()*x ) / b_vector_.y();
    z = sqrt( square( c ) - square( x ) - square( y ) );
    c_vector_ = Vector3D( x, y, z );
    // Build a matrix and invert it to get the reciprocal axes
    fractional_to_orthogonal_matrix_ = Matrix3D( a_vector_.x(), b_vector_.x(), c_vector_.x(),
                                                 a_vector_.y(), b_vector_.y(), c_vector_.y(),
                                                 a_vector_.z(), b_vector_.z(), c_vector_.z()
                                               );
    volume_ = fractional_to_orthogonal_matrix_.determinant();
    orthogonal_to_fractional_matrix_ = fractional_to_orthogonal_matrix_;
    orthogonal_to_fractional_matrix_.invert();
    a_star_vector_ = Vector3D( orthogonal_to_fractional_matrix_.value( 0, 0 ), orthogonal_to_fractional_matrix_.value( 0, 1 ), orthogonal_to_fractional_matrix_.value( 0, 2 ) );
    b_star_vector_ = Vector3D( orthogonal_to_fractional_matrix_.value( 1, 0 ), orthogonal_to_fractional_matrix_.value( 1, 1 ), orthogonal_to_fractional_matrix_.value( 1, 2 ) );
    c_star_vector_ = Vector3D( orthogonal_to_fractional_matrix_.value( 2, 0 ), orthogonal_to_fractional_matrix_.value( 2, 1 ), orthogonal_to_fractional_matrix_.value( 2, 2 ) );
    a_star_ = a_star_vector_.length();
    b_star_ = b_star_vector_.length();
    c_star_ = c_star_vector_.length();
    N_ = SymmetricMatrix3D( a_star_, b_star_, c_star_, 0.0, 0.0, 0.0 );
    N_inverse_ = SymmetricMatrix3D( 1.0/a_star_, 1.0/b_star_, 1.0/c_star_, 0.0, 0.0, 0.0 );
    alpha_star_ = angle( b_star_vector_, c_star_vector_ );
    beta_star_  = angle( a_star_vector_, c_star_vector_ );
    gamma_star_ = angle( a_star_vector_, b_star_vector_ );
    deduce_lattice_system();
    constraints_ = ::constraints( lattice_system_ );
    b_is_constrained_ = ( ( constraints_ & 1 ) == 1 );
    c_is_constrained_ = ( ( constraints_ & 2 ) == 2 );
}

// ********************************************************************************

// Initialises from a CASTEP LATTICE_CART block in a .cell file.
void CrystalLattice::from_CASTEP( const Matrix3D & matrix )
{
    Vector3D a_vector( matrix.value( 0, 0 ), matrix.value( 0, 1 ), matrix.value( 0, 2 ) );
    Vector3D b_vector( matrix.value( 1, 0 ), matrix.value( 1, 1 ), matrix.value( 1, 2 ) );
    Vector3D c_vector( matrix.value( 2, 0 ), matrix.value( 2, 1 ), matrix.value( 2, 2 ) );
    double a = a_vector.length();
    double b = b_vector.length();
    double c = c_vector.length();
    Angle alpha = arccosine( (b_vector*c_vector) / (b*c) );
    Angle beta  = arccosine( (a_vector*c_vector) / (a*c) );
    Angle gamma = arccosine( (a_vector*b_vector) / (a*b) );
    *this = CrystalLattice( a, b, c, alpha, beta, gamma );
}

// ********************************************************************************

double CrystalLattice::orthogonality_defect() const
{
    return ( a_ * b_ * c_ ) / volume_;
}

// ********************************************************************************

Vector3D CrystalLattice::lattice_vector( const size_t index ) const
{
    switch( index )
    {
        case  0 : return a_vector_;
        case  1 : return b_vector_;
        case  2 : return c_vector_;
        default : throw std::runtime_error( "CrystalLattice::lattice_vector( size_t ): error: index out of range." );
    }
}

// ********************************************************************************

void CrystalLattice::enclosing_box( Vector3D & min_min_min, Vector3D & max_max_max ) const
{
    min_min_min = Vector3D();
    max_max_max = Vector3D();
    Vector3D current_point;
//    current_point = 0.0 * a_vector() + 0.0 * b_vector() + 0.0 * c_vector();
//    min_min_min.set_x( std::min( min_min_min.x(), current_point.x() ) );
//    min_min_min.set_y( std::min( min_min_min.y(), current_point.y() ) );
//    min_min_min.set_z( std::min( min_min_min.z(), current_point.z() ) );
//    max_max_max.set_x( std::max( max_max_max.x(), current_point.x() ) );
//    max_max_max.set_y( std::max( max_max_max.y(), current_point.y() ) );
//    max_max_max.set_z( std::max( max_max_max.z(), current_point.z() ) );

    current_point = c_vector();
    min_min_min.set_x( std::min( min_min_min.x(), current_point.x() ) );
    min_min_min.set_y( std::min( min_min_min.y(), current_point.y() ) );
    min_min_min.set_z( std::min( min_min_min.z(), current_point.z() ) );
    max_max_max.set_x( std::max( max_max_max.x(), current_point.x() ) );
    max_max_max.set_y( std::max( max_max_max.y(), current_point.y() ) );
    max_max_max.set_z( std::max( max_max_max.z(), current_point.z() ) );

    current_point = b_vector();
    min_min_min.set_x( std::min( min_min_min.x(), current_point.x() ) );
    min_min_min.set_y( std::min( min_min_min.y(), current_point.y() ) );
    min_min_min.set_z( std::min( min_min_min.z(), current_point.z() ) );
    max_max_max.set_x( std::max( max_max_max.x(), current_point.x() ) );
    max_max_max.set_y( std::max( max_max_max.y(), current_point.y() ) );
    max_max_max.set_z( std::max( max_max_max.z(), current_point.z() ) );

    current_point = b_vector() + c_vector();
    min_min_min.set_x( std::min( min_min_min.x(), current_point.x() ) );
    min_min_min.set_y( std::min( min_min_min.y(), current_point.y() ) );
    min_min_min.set_z( std::min( min_min_min.z(), current_point.z() ) );
    max_max_max.set_x( std::max( max_max_max.x(), current_point.x() ) );
    max_max_max.set_y( std::max( max_max_max.y(), current_point.y() ) );
    max_max_max.set_z( std::max( max_max_max.z(), current_point.z() ) );

    current_point = a_vector();
    min_min_min.set_x( std::min( min_min_min.x(), current_point.x() ) );
    min_min_min.set_y( std::min( min_min_min.y(), current_point.y() ) );
    min_min_min.set_z( std::min( min_min_min.z(), current_point.z() ) );
    max_max_max.set_x( std::max( max_max_max.x(), current_point.x() ) );
    max_max_max.set_y( std::max( max_max_max.y(), current_point.y() ) );
    max_max_max.set_z( std::max( max_max_max.z(), current_point.z() ) );

    current_point = a_vector() + c_vector();
    min_min_min.set_x( std::min( min_min_min.x(), current_point.x() ) );
    min_min_min.set_y( std::min( min_min_min.y(), current_point.y() ) );
    min_min_min.set_z( std::min( min_min_min.z(), current_point.z() ) );
    max_max_max.set_x( std::max( max_max_max.x(), current_point.x() ) );
    max_max_max.set_y( std::max( max_max_max.y(), current_point.y() ) );
    max_max_max.set_z( std::max( max_max_max.z(), current_point.z() ) );

    current_point = a_vector() + b_vector();
    min_min_min.set_x( std::min( min_min_min.x(), current_point.x() ) );
    min_min_min.set_y( std::min( min_min_min.y(), current_point.y() ) );
    min_min_min.set_z( std::min( min_min_min.z(), current_point.z() ) );
    max_max_max.set_x( std::max( max_max_max.x(), current_point.x() ) );
    max_max_max.set_y( std::max( max_max_max.y(), current_point.y() ) );
    max_max_max.set_z( std::max( max_max_max.z(), current_point.z() ) );

    current_point = a_vector() + b_vector() + c_vector();
    min_min_min.set_x( std::min( min_min_min.x(), current_point.x() ) );
    min_min_min.set_y( std::min( min_min_min.y(), current_point.y() ) );
    min_min_min.set_z( std::min( min_min_min.z(), current_point.z() ) );
    max_max_max.set_x( std::max( max_max_max.x(), current_point.x() ) );
    max_max_max.set_y( std::max( max_max_max.y(), current_point.y() ) );
    max_max_max.set_z( std::max( max_max_max.z(), current_point.z() ) );

}

// ********************************************************************************

Matrix3D CrystalLattice::metric_matrix() const
{
    return Matrix3D( a_ * a_                  , a_ * b_ * gamma_.cosine() , a_ * c_ * beta_.cosine() ,
                     a_ * b_ * gamma_.cosine(), b_ * b_                   , b_ * c_ * alpha_.cosine(),
                     a_ * c_ * beta_.cosine() , b_ * c_ * alpha_.cosine() , c_ * c_ );
}

// ********************************************************************************

Matrix3D CrystalLattice::for_CASTEP() const
{
    double x;
    double y;
    double z;
    z = c_;
    y = 0.0;
    x = 0.0;
    Vector3D c_vector = Vector3D( x, y, z );
    
    z = b_ * alpha_.cosine();
    y = b_ * alpha_.sine();
    x = 0.0;
    Vector3D b_vector = Vector3D( x, y, z );
    
    z = a_ * beta_.cosine();
    y = ( b_*a_*gamma_.cosine() - b_vector.z()*z ) / b_vector.y();
    x = sqrt( square( a_ ) - square( z ) - square( y ) );
    Vector3D a_vector = Vector3D( x, y, z );
    
    return Matrix3D( a_vector.x(), a_vector.y(), a_vector.z(),
                     b_vector.x(), b_vector.y(), b_vector.z(),
                     c_vector.x(), c_vector.y(), c_vector.z()
                   );
}

// ********************************************************************************

Vector3D CrystalLattice::orthogonal_to_fractional( const Vector3D & input ) const
{
    return orthogonal_to_fractional_matrix_ * input;
}

// ********************************************************************************

Vector3D CrystalLattice::fractional_to_orthogonal( const Vector3D & input ) const
{
    return fractional_to_orthogonal_matrix_ * input;
}

// ********************************************************************************

void CrystalLattice::rescale_volume( const double target_volume, size_t Z )
{
    size_t current_Z(1);
    if ( Z == 0 )
        Z = 1;
    else
        current_Z = round_to_int( ( volume() / target_volume ) * Z );
    double k = std::pow( (target_volume/Z) / (volume()/current_Z), 1.0/3.0 );
    *this = CrystalLattice( a()*k, b()*k, c()*k, alpha(), beta(), gamma() );
}

// ********************************************************************************

// Finds shortest distance, in Angstrom, between two positions given in fractional coordinates.
double CrystalLattice::shortest_distance( const Vector3D & lhs, const Vector3D & rhs ) const
{
    return sqrt( shortest_distance2( lhs, rhs ) );
}

// ********************************************************************************

// Finds shortest distance, in Angstrom^2, between two positions given in fractional coordinates.
double CrystalLattice::shortest_distance2( const Vector3D & lhs, const Vector3D & rhs ) const
{
    Vector3D difference_vector = adjust_for_translations( rhs - lhs ); // In fractional coordinates
    // Now we must find the shortest distance.
    double shortest_distance2 = fractional_to_orthogonal( difference_vector ).norm2();
    // "shortest_distance" may now be something like 0.95, which clearly should have been 0.05. Likewise,
    // with very acute unit-cell angles, it may be necessary to add or subtract +/- 1 (fractional coordinates).
    Vector3D original_difference_vector( difference_vector );
    bool shortest_distance_changed( false );
    do
    {
        shortest_distance_changed = false;
        for ( int i( -1 ); i != 2; ++i )
        {
            for ( int j( -1 ); j != 2; ++j )
            {
                for ( int k( -1 ); k != 2; ++k )
                {
                    Vector3D new_difference_vector = original_difference_vector + Vector3D( i, j, k );
                    double distance2 = fractional_to_orthogonal( new_difference_vector ).norm2();
                    if ( distance2 < shortest_distance2 )
                    {
                        difference_vector = new_difference_vector;
                        shortest_distance2 = distance2;
                        shortest_distance_changed = true;
                    }
                }
            }
        }
    }
    while ( shortest_distance_changed );
    return shortest_distance2;
}

// ********************************************************************************

// Finds shortest distance, in Angstrom, between two positions given in fractional coordinates.
// Returns the shortest distance and the shortest difference vector (in fractional coordinates).
void CrystalLattice::shortest_distance( const Vector3D & lhs, const Vector3D & rhs, double & output_distance, Vector3D & output_difference_vector ) const
{
    Vector3D difference_vector = adjust_for_translations( rhs - lhs ); // In fractional coordinates
    // Now we must find the shortest distance.
    double shortest_distance = fractional_to_orthogonal( difference_vector ).norm2();
    // "shortest_distance" may now be something like 0.95, which clearly should have been 0.05. Likewise,
    // with very acute unit-cell angles, it may be necessary to add or subtract +/- 1 (fractional coordinates).
    Vector3D original_difference_vector( difference_vector );
    bool shortest_distance_changed( false );
    do
    {
        shortest_distance_changed = false;
        for ( int i( -1 ); i != 2; ++i )
        {
            for ( int j( -1 ); j != 2; ++j )
            {
                for ( int k( -1 ); k != 2; ++k )
                {
                    Vector3D new_difference_vector = original_difference_vector + Vector3D( i, j, k );
                    double distance = fractional_to_orthogonal( new_difference_vector ).norm2();
                    if ( distance < shortest_distance )
                    {
                        difference_vector = new_difference_vector;
                        shortest_distance = distance;
                        shortest_distance_changed = true;
                    }
                }
            }
        }
    }
    while ( shortest_distance_changed );
    output_distance = sqrt( shortest_distance );
    output_difference_vector = difference_vector;
}

// ********************************************************************************

void CrystalLattice::set_lattice_system( const LatticeSystem lattice_system )
{
    // It should be possible to start with the unit-cell parameters of a cubic unit cell,
    // then manually set it to triclinic and then manually set it to cubic again.
    char constraints_from_unit_cell_parameters = ::constraints( a_, b_, c_, alpha_, beta_, gamma_ );
    char constraints_from_lattice_system = ::constraints( lattice_system );
    if ( ( constraints_from_lattice_system & constraints_from_unit_cell_parameters ) != constraints_from_lattice_system )
        throw std::runtime_error( "CrystalLattice::set_lattice_system(): error: unit-cell parameters incommensurate with lattice system." );
    lattice_system_ = lattice_system;
    constraints_ = constraints_from_lattice_system;
    b_is_constrained_ = ( ( constraints_ & 1 ) == 1 );
    c_is_constrained_ = ( ( constraints_ & 2 ) == 2 );
}

// ********************************************************************************

Matrix3D CrystalLattice::choose_angles_close_to_90() const
{
    Matrix3D result;
    // TRICLINIC, MONOCLINIC, ORTHORHOMBIC, TETRAGONAL, HEXAGONAL, RHOMBOHEDRAL, CUBIC
    if ( lattice_system_ == CUBIC || lattice_system_ == HEXAGONAL || lattice_system_ == TETRAGONAL || lattice_system_ == ORTHORHOMBIC )
    {
        std::cout << "CrystalLattice::choose_angles_close_to_90(): warning: lattice was already known to be optimal." << std::endl;
        return result;
    }
    // TRICLINIC, MONOCLINIC, RHOMBOHEDRAL
    if ( lattice_system_ == MONOCLINIC_A )
    {
        // Just use Gram-Schmidt and round to the nearest integer.
        // Add the shortest vector to the longest.
        if ( b_ < c_ )
        {
            double omega = ( c_vector_ * b_vector_ ) / ( b_ * b_ );
            int nomegas = round_to_int( omega );
            if ( nomegas != 0 )
                result.set_value( 1, 2, -nomegas );
        }
        else
        {
            double omega = ( c_vector_ * b_vector_ ) / ( c_ * c_ );
            int nomegas = round_to_int( omega );
            if ( nomegas != 0 )
                result.set_value( 2, 1, -nomegas );
        }
        return result;
    }
    if ( lattice_system_ == MONOCLINIC_B )
    {
        if ( a_ < c_ )
        {
            double omega = ( c_vector_ * a_vector_ ) / ( a_ * a_ );
            int nomegas = round_to_int( omega );
            if ( nomegas != 0 )
                result.set_value( 0, 2, -nomegas );
        }
        else
        {
            double omega = ( c_vector_ * a_vector_ ) / ( c_ * c_ );
            int nomegas = round_to_int( omega );
            if ( nomegas != 0 )
                result.set_value( 2, 0, -nomegas );
        }
        return result;
    }
    if ( lattice_system_ == MONOCLINIC_C )
    {
        if ( a_ < b_ )
        {
            double omega = ( b_vector_ * a_vector_ ) / ( a_ * a_ );
            int nomegas = round_to_int( omega );
            if ( nomegas != 0 )
                result.set_value( 0, 1, -nomegas );
        }
        else
        {
            double omega = ( b_vector_ * a_vector_ ) / ( b_ * b_ );
            int nomegas = round_to_int( omega );
            if ( nomegas != 0 )
                result.set_value( 1, 0, -nomegas );
        }
        return result;
    }
    // TRICLINIC, RHOMBOHEDRAL
    CrystalLattice new_lattice( *this );
  //  {
        double best_orthogonality_defect = orthogonality_defect();
        std::cout << "orthogonality defect before = " << best_orthogonality_defect << std::endl;
        int limit = 3;
        for ( int i1( -limit ); i1 != limit+1; ++i1 )
        {
        for ( int i2( -limit ); i2 != limit+1; ++i2 )
        {
        for ( int i3( -limit ); i3 != limit+1; ++i3 )
        {
            for ( int j1( -limit ); j1 != limit+1; ++j1 )
            {
            for ( int j2( -limit ); j2 != limit+1; ++j2 )
            {
            for ( int j3( -limit ); j3 != limit+1; ++j3 )
            {
                for ( int k1( -limit ); k1 != limit+1; ++k1 )
                {
                for ( int k2( -limit ); k2 != limit+1; ++k2 )
                {
                for ( int k3( -limit ); k3 != limit+1; ++k3 )
                {
                    // Make a copy.
                    CrystalLattice new_lattice( *this );
                    // Transform.
                    Matrix3D transformation_matrix( i1, i2, i3, j1, j2, j3, k1, k2, k3 );
                    if ( ! nearly_equal( transformation_matrix.determinant(), 1.0 ) )
                        continue;
                    new_lattice.transform( transformation_matrix );
                    if ( new_lattice.orthogonality_defect() < best_orthogonality_defect )
                    {
                        result = transformation_matrix;
                        best_orthogonality_defect = new_lattice.orthogonality_defect();
                    }
                }
                }
                }
            }
            }
            }
        }
        }
        }

if ( false )
{
        std::vector< double > unit_cell_lengths;
        unit_cell_lengths.push_back( a_ );
        unit_cell_lengths.push_back( b_ );
        unit_cell_lengths.push_back( c_ );
        Mapping sorted_map = sort( unit_cell_lengths );
        Vector3D r_s = new_lattice.lattice_vector( sorted_map[0] ); // s = short
        Vector3D r_l = new_lattice.lattice_vector( sorted_map[2] ); // l = long
        double omega = ( r_l * r_s ) / ( r_s * r_s ); // r_s*r_s is just a^2, b^2 or c^2.
        int nomegas = round_to_int( omega );
        if ( nomegas != 0 )
            result.set_value( sorted_map[2], sorted_map[0], -nomegas );
}
    std::cout << "orthogonality defect after = " << best_orthogonality_defect << std::endl;
    return result;
}

// ********************************************************************************

void CrystalLattice::transform( const Matrix3D & m )
{
    // In practice, the elements of the transformation matrix will be integers like 0, -1, 1 and it would be cleaner
    // to test for closeness to 0 and then to discard the contribution.
    if ( ! nearly_equal( m.determinant(), 1.0 ) )
        std::cout << "CrystalLattice::transform(): determinant = " + double2string( m.determinant() ) << std::endl;
    Vector3D new_a = m.value( 0, 0 ) * a_vector_ + m.value( 0, 1 ) * b_vector_ + m.value( 0, 2 ) * c_vector_;
    Vector3D new_b = m.value( 1, 0 ) * a_vector_ + m.value( 1, 1 ) * b_vector_ + m.value( 1, 2 ) * c_vector_;
    Vector3D new_c = m.value( 2, 0 ) * a_vector_ + m.value( 2, 1 ) * b_vector_ + m.value( 2, 2 ) * c_vector_;
    *this = CrystalLattice( new_a.length(), new_b.length(), new_c.length(), angle( new_b, new_c ), angle( new_a, new_c ), angle( new_a, new_b ) );
}

// ********************************************************************************

void CrystalLattice::print() const
{
    std::cout << "a = " << a() << ", " <<
                 "b = " << b() << ", " <<
                 "c = " << c() << ", " <<
                 "al = " << alpha() << ", " <<
                 "be = " << beta()  << ", " <<
                 "ga = " << gamma() << std::endl;
}

// ********************************************************************************

void CrystalLattice::show() const
{
    std::cout << "a = " << a_vector_;
    std::cout << "b = " << b_vector_;
    std::cout << "c = " << c_vector_;
    std::cout << "a* = " << a_star_vector_;
    std::cout << "b* = " << b_star_vector_;
    std::cout << "c* = " << c_star_vector_;
}

// ********************************************************************************

Matrix3D CrystalLattice::Downs_D() const
{
    return Matrix3D( a_vector_.x(), b_vector_.x(), c_vector_.x(),
                     a_vector_.y(), b_vector_.y(), c_vector_.y(),
                     a_vector_.z(), b_vector_.z(), c_vector_.z() );
}

// ********************************************************************************
Matrix3D CrystalLattice::Downs_D_star() const
{
    return Matrix3D( a_star_vector_.x(), b_star_vector_.x(), c_star_vector_.x(),
                     a_star_vector_.y(), b_star_vector_.y(), c_star_vector_.y(),
                     a_star_vector_.z(), b_star_vector_.z(), c_star_vector_.z() );
}

// ********************************************************************************
Matrix3D CrystalLattice::Downs_G() const
{
    return Matrix3D( a_ * a_                  , a_ * b_ * gamma_.cosine() , a_ * c_ * beta_.cosine() ,
                     a_ * b_ * gamma_.cosine(), b_ * b_                   , b_ * c_ * alpha_.cosine(),
                     a_ * c_ * beta_.cosine() , b_ * c_ * alpha_.cosine() , c_ * c_ );
}

// ********************************************************************************

Matrix3D CrystalLattice::Downs_G_star() const
{
    return Matrix3D( a_star_ * a_star_                       , a_star_ * b_star_ * gamma_star_.cosine() , a_star_ * c_star_ * beta_star_.cosine() ,
                     a_star_ * b_star_ * gamma_star_.cosine(), b_star_ * b_star_                        , b_star_ * c_star_ * alpha_star_.cosine(),
                     a_star_ * c_star_ * beta_star_.cosine() , b_star_ * c_star_ * alpha_star_.cosine() , c_star_ * c_star_ );
}

// ********************************************************************************

void CrystalLattice::deduce_lattice_system()
{
    bool angles_equal = triquality( alpha(), beta(), gamma() );
    bool ab_equal = nearly_equal( a(), b() );
    bool alpha_is_90 = alpha().nearly_90();
// { TRICLINIC, MONOCLINIC, ORTHORHOMBIC, TETRAGONAL, HEXAGONAL, RHOMBOHEDRAL, CUBIC }
    if ( angles_equal )
    {
        bool lengths_equal = triquality( a(), b(), c() );
        if ( alpha_is_90 )
        {
            if ( lengths_equal )
            {
                lattice_system_ = CUBIC;
                return;
            }
            if ( ab_equal )
                lattice_system_ = TETRAGONAL;
            else
                lattice_system_ = ORTHORHOMBIC;
            return;
        }
        if ( lengths_equal )
        {
                lattice_system_ = RHOMBOHEDRAL;
                return;
        }
        std::cout << "CrystalLattice::deduce_lattice_system(): Warning: angles are all equal, but system is triclinic." << std::endl;
        lattice_system_ = TRICLINIC;
        return;
    }
// { TRICLINIC, MONOCLINIC, HEXAGONAL }
    bool beta_is_90 = beta().nearly_90();
    if ( alpha_is_90 && beta_is_90 )
    {
        if ( ab_equal && nearly_equal( gamma(), Angle::angle_120_degrees() ) )
            lattice_system_ = HEXAGONAL;
        else
            lattice_system_ = MONOCLINIC_C;
        return;
    }
// { TRICLINIC, MONOCLINIC }
    if ( gamma().nearly_90() )
    {
        if ( beta_is_90  )
        {
            lattice_system_ = MONOCLINIC_A;
            return;
        }
        if ( alpha_is_90  )
        {
            lattice_system_ = MONOCLINIC_B;
            return;
        }
    }
    lattice_system_ = TRICLINIC;
}

// ********************************************************************************

std::string LatticeSystem2string( const CrystalLattice::LatticeSystem lattice_system )
{
    switch ( lattice_system )
    {
        case CrystalLattice::TRICLINIC    : return "Triclinic";
        case CrystalLattice::MONOCLINIC_A : return "Monoclinic";
        case CrystalLattice::MONOCLINIC_B : return "Monoclinic";
        case CrystalLattice::MONOCLINIC_C : return "Monoclinic";
        case CrystalLattice::ORTHORHOMBIC : return "Orthorhombic";
        case CrystalLattice::TETRAGONAL   : return "Tetragonal";
        case CrystalLattice::HEXAGONAL    : return "Hexagonal";
        case CrystalLattice::RHOMBOHEDRAL : return "Rhombohedral";
        case CrystalLattice::CUBIC        : return "Cubic";
        default                           : throw std::runtime_error( "CrystalLattice::set_lattice_system(): error: unit-cell parameters incommensurate with lattice system." );
    }
}

// ********************************************************************************

// 1.  1 b = a.
// 2.  2 All lengths are equal.
// 3.  4 All angles are equal.
// 4.  8 alpha = 90.
// 5. 16 beta = 90.
// 6. 32 gamma = 90.
// 7. 64 gamma = 120.
char constraints( const CrystalLattice::LatticeSystem lattice_system )
{
    switch ( lattice_system )
    {
        case CrystalLattice::TRICLINIC    : return 0;
        case CrystalLattice::MONOCLINIC_A : return 16+32;
        case CrystalLattice::MONOCLINIC_B : return 8+32;
        case CrystalLattice::MONOCLINIC_C : return 8+16;
        case CrystalLattice::ORTHORHOMBIC : return 4+8+16+32;
        case CrystalLattice::TETRAGONAL   : return 1+4+8+16+32;
        case CrystalLattice::HEXAGONAL    : return 1+4+8+16+32;
        case CrystalLattice::RHOMBOHEDRAL : return 1+2+4;
        case CrystalLattice::CUBIC        : return 1+2+4+8+16+32;
        default                           : throw std::runtime_error( "CrystalLattice::set_lattice_system(): error: unit-cell parameters incommensurate with lattice system." );
    }
}

// ********************************************************************************

char constraints( const double a, const double b, const double c, const Angle alpha, const Angle beta, const Angle gamma )
{
    bool angles_equal = triquality( alpha, beta, gamma );
    bool ab_equal = nearly_equal( a, b );
    bool alpha_is_90 = alpha.nearly_90();
// { TRICLINIC, MONOCLINIC, ORTHORHOMBIC, TETRAGONAL, HEXAGONAL, RHOMBOHEDRAL, CUBIC }
    if ( angles_equal )
    {
        bool lengths_equal = triquality( a, b, c );
        if ( alpha_is_90 )
        {
            if ( lengths_equal )
                return constraints( CrystalLattice::CUBIC );
            if ( ab_equal )
                return constraints( CrystalLattice::TETRAGONAL );
            else
                return constraints( CrystalLattice::ORTHORHOMBIC );
        }
        if ( lengths_equal )
            return constraints( CrystalLattice::RHOMBOHEDRAL );
        else
            return constraints( CrystalLattice::TRICLINIC );
    }
// { TRICLINIC, MONOCLINIC, HEXAGONAL }
    bool beta_is_90 = beta.nearly_90();
    if ( alpha_is_90 && beta_is_90 )
    {
        if ( ab_equal && nearly_equal( gamma, Angle::angle_120_degrees() ) )
            return constraints( CrystalLattice::HEXAGONAL );
        else
            return constraints( CrystalLattice::MONOCLINIC_C );
    }
// { TRICLINIC, MONOCLINIC }
    if ( gamma.nearly_90() )
    {
        if ( beta_is_90  )
            return constraints( CrystalLattice::MONOCLINIC_A );
        if ( alpha_is_90  )
            return constraints( CrystalLattice::MONOCLINIC_B );
    }
    return constraints( CrystalLattice::TRICLINIC );
}

// ********************************************************************************

CrystalLattice average( const CrystalLattice & lhs, const CrystalLattice & rhs, const double weight )
{
    if ( ! nearly_integer( weight ) )
        std::cout << "average( CrystalLattice, CrystalLattice, double ): warning: weight is not an integer, is that intended?" << std::endl;
    CrystalLattice result( average( lhs.a(), rhs.a(), weight ),
                           average( lhs.b(), rhs.b(), weight ),
                           average( lhs.c(), rhs.c(), weight ),
                           average( lhs.alpha(), rhs.alpha(), weight ),
                           average( lhs.beta() , rhs.beta() , weight ),
                           average( lhs.gamma(), rhs.gamma(), weight ) );
    if ( ! nearly_equal( lhs, result ) )
        std::cout << "average( CrystalLattice, CrystalLattice, double ): warning: lattices differ."<< std::endl;
    return result;
}

// ********************************************************************************

bool nearly_equal( const CrystalLattice & lhs, const CrystalLattice & rhs, double length_tolerance_percentage, const Angle angle_tolerance )
{
    if ( length_tolerance_percentage < 0.0 )
        throw std::runtime_error( "nearly_equal( CrystalLattice, CrystalLattice ): error: length tolerance is negative." );
    if ( length_tolerance_percentage > 100.0 )
        throw std::runtime_error( "nearly_equal( CrystalLattice, CrystalLattice ): error: length tolerance > 100%." );
    if ( angle_tolerance < Angle::from_degrees( 0.0 ) )
        throw std::runtime_error( "nearly_equal( CrystalLattice, CrystalLattice ): error: angle tolerance is negative." );
    if ( angle_tolerance > Angle::angle_90_degrees() )
        throw std::runtime_error( "nearly_equal( CrystalLattice, CrystalLattice ): error: angle tolerance > 90 degrees." );
    length_tolerance_percentage = length_tolerance_percentage / 100.0;
    if ( absolute_relative_difference( lhs.a(), rhs.a() ) > length_tolerance_percentage )
        return false;
    if ( absolute_relative_difference( lhs.b(), rhs.b() ) > length_tolerance_percentage )
        return false;
    if ( absolute_relative_difference( lhs.c(), rhs.c() ) > length_tolerance_percentage )
        return false;
    if ( ! nearly_equal( lhs.alpha(), rhs.alpha(), angle_tolerance ) )
        return false;
    if ( ! nearly_equal( lhs.beta(), rhs.beta(), angle_tolerance ) )
        return false;
    if ( ! nearly_equal( lhs.gamma(), rhs.gamma(), angle_tolerance ) )
        return false;
    return true;
}

// ********************************************************************************

