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

#include "CrystalStructure.h"
#include "3DCalculations.h"
#include "BasicMathsFunctions.h"
#include "ChemicalFormula.h"
#include "ConnectivityTable.h"
#include "FileName.h"
#include "Mapping.h"
#include "PhysicalConstants.h"
#include "PointGroup.h"
#include "RunningAverageAndESD.h"
#include "TextFileWriter.h"
#include "Utilities.h"

#include <iostream>
#include <stdexcept>

// ********************************************************************************

CrystalStructure::CrystalStructure(): space_group_symmetry_has_been_applied_(false)
{
}

// ********************************************************************************

Atom CrystalStructure::atom( const size_t i ) const
{
    if ( i >= atoms_.size() )
        throw std::runtime_error( "CrystalStructure::atom( size_t ): i >= atoms_.size()" );
    return atoms_[i];
}

// ********************************************************************************

ChemicalFormula CrystalStructure::chemical_formula() const
{
    ChemicalFormula result;
    for ( size_t i( 0 ); i != natoms(); ++i )
        result.add_element( atoms_[i].element() );
    return result;
}

// ********************************************************************************

void CrystalStructure::add_atom( const Atom & atom )
{
    atoms_.push_back( atom );
    suppressed_.push_back( false );
    basic_checks();
}

// ********************************************************************************

void CrystalStructure::add_atoms( const std::vector< Atom > & atoms )
{
    atoms_.reserve( atoms_.size() + atoms.size() );
    atoms_.insert( atoms_.end(), atoms.begin(), atoms.end() );
    for ( size_t i( 0 ); i != atoms.size(); ++i )
        suppressed_.push_back( false );
    basic_checks();
}

// ********************************************************************************

void CrystalStructure::remove_H_and_D()
{
    std::vector< Atom > new_atoms;
    for ( size_t i( 0 ); i != natoms(); ++i )
    {
        if ( ! atom( i ).element().is_H_or_D() )
            new_atoms.push_back( atom( i ) );
    }
    atoms_ = new_atoms;
}

// ********************************************************************************

void CrystalStructure::set_atom( const size_t i, const Atom & atom )
{
    atoms_[i] = atom;
    basic_checks();
}

// ********************************************************************************

void CrystalStructure::basic_checks() const
{

}

// ********************************************************************************

size_t CrystalStructure::find_label( const std::string & label ) const
{
    for ( size_t i( 0 ); i != atoms_.size(); ++i )
    {
        if ( atoms_[i].label() == label )
            return i;
    }
    return natoms();
}

// ********************************************************************************

size_t CrystalStructure::atom( const std::string & atom_label ) const
{
    for ( size_t i( 0 ); i != atoms_.size(); ++i )
    {
        if ( atoms_[ i ].label() == atom_label )
            return i;
    }
    throw std::runtime_error( "CrystalStructure::atom( const std::string & label ): label >" + atom_label + "< not found." );
}

// ********************************************************************************

void CrystalStructure::make_atom_labels_unique()
{
    for ( size_t i( 0 ); i != atoms_.size(); ++i )
        atoms_[ i ].set_label( atoms_[ i ].element().symbol() + size_t2string( i ) );
}

// ********************************************************************************

std::set< Element > CrystalStructure::elements() const
{
    std::set< Element > result;
    for ( std::vector< Atom >::const_iterator it( atoms_.begin() ); it != atoms_.end(); ++it )
        result.insert( it->element() );
    return result;
}

// ********************************************************************************

void CrystalStructure::list_all_bonds( std::vector< std::string > & labels_1, std::vector< std::string > & labels_2, std::vector< double > & bonds ) const
{
    
    for ( size_t i( 0 ); i != connectivity_table_.size(); ++i )
    {
        for ( size_t j( i + 1 ); j != connectivity_table_.size(); ++j )
        {
            if ( connectivity_table_.value( i, j ) != 0 )
            {
                labels_1.push_back( atom( i ).label() );
                labels_2.push_back( atom( j ).label() );
                double distance;
                Vector3D difference_vector; // @@ This is a bit wasteful...
                crystal_lattice_.shortest_distance( atom( i ).position(), atom( j ).position(), distance, difference_vector );
                bonds.push_back( distance );
            }
        }
    }
}

// ********************************************************************************

void CrystalStructure::list_all_angles( std::vector< std::string > & labels_1, std::vector< std::string > & labels_2, std::vector< std::string > & labels_3, std::vector< double > & angles ) const
{
    
}

// ********************************************************************************

// Only the asymmetric unit is kept, everything else is deleted.
void CrystalStructure::reduce_to_asymmetric_unit( const double tolerance )
{
    std::vector< Atom > new_atoms;
    std::vector< bool > is_duplicate( natoms(), false );
    for ( size_t i( 0 ); i != natoms(); ++i )
    {
        if ( is_duplicate[i] )
            continue;
        for ( size_t j( i+1 ); j != natoms(); ++j )
        {
            if ( is_duplicate[j] )
                continue;
            if ( atom(i).element() != atom(j).element() )
                continue;
            if ( shortest_distance2( atom(i).position(), atom(j).position() ) < square( tolerance ) )
            {
                if ( ! nearly_equal( atom(i).occupancy(), atom(j).occupancy() ) )
                    std::cout << "CrystalStructure::reduce_to_asymmetric_unit(): warning: duplicate atoms have different occupancies." << std::endl;
                is_duplicate[j] = true;
            }
        }
        new_atoms.push_back( atom(i) );
    }
    atoms_ = new_atoms;
    space_group_symmetry_has_been_applied_ = false;
}

// ********************************************************************************

void CrystalStructure::apply_space_group_symmetry()
{
    if ( space_group_symmetry_has_been_applied_ )
        std::cout << "CrystalStructure::apply_space_group_symmetry(): WARNING: space group has already been applied." << std::endl;
    std::vector< Atom > atoms;

    for ( size_t j( 1 ); j != space_group_.nsymmetry_operators(); ++j )
    {
        for ( size_t i( 0 ); i != natoms(); ++i )
        {
            Vector3D original_position = atom( i ).position();
            Vector3D new_position = space_group_.symmetry_operator( j ) * original_position;
            double distance = crystal_lattice_.shortest_distance( original_position, new_position );
            // Is it a special position?
            if ( distance > 0.1 )
            {
                Atom new_atom = atom( i );
                new_atom.set_position( new_position );
                if ( new_atom.ADPs_type() == Atom::ANISOTROPIC )
                {
                    new_atom.set_anisotropic_displacement_parameters( rotate_adps( new_atom.anisotropic_displacement_parameters(), space_group_.symmetry_operator( j ).rotation(), crystal_lattice_ ) );
                }
                new_atom.set_label( atom( i ).label() + "_" + size_t2string( j ) );
                atoms.push_back( new_atom );
            }
        }
    }

    add_atoms( atoms );
    space_group_symmetry_has_been_applied_ = true;
}

// ********************************************************************************

void CrystalStructure::move_atoms_to_form_molecules( const bool include_symmetry_operators )
{
    std::vector< bool > has_been_connected( natoms(), false );
    connectivity_table_ = ConnectivityTable( natoms() );
    for ( size_t i( 0 ); i != natoms(); ++i )
    {
        Atom iAtom = atom( i );
        for ( size_t j( i+1 ); j != natoms(); ++j )
        {
            Atom jAtom = atom( j );
            double distance2;
            if ( include_symmetry_operators )
                distance2 = shortest_distance2( iAtom.position(), jAtom.position() );
            else
                distance2 = crystal_lattice_.shortest_distance2( iAtom.position(), jAtom.position() );
            if ( are_bonded( iAtom.element(), jAtom.element(), distance2 ) )
            {
                // Add this one to the connectivity table.
                connectivity_table_.set_value( i, j, 1 );
                // Move atom j so that it really bonds to atom i
                double distance;
                Vector3D difference_vector;
                if ( include_symmetry_operators )
                    shortest_distance( iAtom.position(), jAtom.position(), distance, difference_vector );
                else
                    crystal_lattice_.shortest_distance( iAtom.position(), jAtom.position(), distance, difference_vector );
                if ( ! nearly_equal( jAtom.position(), iAtom.position() + difference_vector ) )
                {
                    // If we are here, we have to move atom j to connect it to atom i
                    if ( ! has_been_connected[j] )
                    {
                        jAtom.set_position( iAtom.position() + difference_vector );
                        set_atom( j, jAtom );
                    }
                    else if ( ! has_been_connected[i] )
                    {
                        iAtom.set_position( jAtom.position() - difference_vector );
                        set_atom( i, iAtom );
                    }
                    else
                        throw std::runtime_error( "CrystalStructure::move_atoms_to_form_molecules(): atoms i and j have both been moved but are not bonded." );
                }
                has_been_connected[i] = true;
                has_been_connected[j] = true;
            }
        }
    }
}

// ********************************************************************************

void CrystalStructure::perceive_molecules( const bool I_know_Zprime_is_one )
{
    // In crystals with a floating axis, mapping two crystals requires a shift over the c.o.m. of the molecule.
    // The following two commands are absolutely necessary to avoid a number of difficult complications.
    // 1. .cif files saved by Mercury probably have molecules on special positions expanded into full molecules.
    // So we cannot rely on the cif only containing the asymmetric unit, but we cannot rely on the cif containing
    // expanded molecules either.
    // 2. If the crystal structure is a polymer, expanding the asymmetric unit to build a molecule would never terminate.
    // (we could add some kind of loop counter, or a test to see if an atom we are about to add is in fact a translationally
    // equivalent atom of an atom that is already part of the molecule).
    // The following two commands guarantee that our list of atoms consists of exactly the atoms that fill one unit cell.
    if ( ! I_know_Zprime_is_one )
    {
        reduce_to_asymmetric_unit();
        apply_space_group_symmetry();
    }
    move_atoms_to_form_molecules( false );
    std::vector< std::vector< size_t > > molecules = split( connectivity_table_ );
    for ( size_t i( 0 ); i != molecules.size(); ++i )
    {
        MoleculeInCrystal molecule_in_crystal;
        for ( size_t j( 0 ); j != molecules[i].size(); ++j )
        {
            molecule_in_crystal.add_atom( atom( molecules[i][j] ) );
        }
        molecules_.push_back( molecule_in_crystal );
    }
}

// ********************************************************************************

// @@ This requires that you run the molecule perception method first
void CrystalStructure::remove_symmetry_related_molecules()
{
    
}

// ********************************************************************************

// This returns a copy, so would copy all atoms.
MoleculeInCrystal CrystalStructure::molecule_in_crystal( const size_t i ) const
{
    if ( i < molecules_.size() )
        return molecules_[i];
    throw std::runtime_error( "CrystalStructure::molecule_in_crystal(): i >= molecules_.size()." );
}

// ********************************************************************************

// point is in fractional coordinates.
// Tolerance is in Angstrom.
// On return, point contains the point moved to the exact special position.
PointGroup CrystalStructure::point_is_on_special_position( Vector3D & point, const double tolerance ) const
{
    std::vector< bool > used( space_group_.nsymmetry_operators(), false );
    bool a_change_was_made( true );
    while ( a_change_was_made )
    {
        a_change_was_made = false;
        // We start at index 1 because we skip the identity
        for ( size_t i( 1 ); i != space_group_.nsymmetry_operators(); ++i )
        {
            if ( used[i] )
                continue;
            Vector3D point_2 = space_group_.symmetry_operator( i ) * point;
            double distance;
            Vector3D difference_vector;
            crystal_lattice_.shortest_distance( point, point_2, distance, difference_vector );
            if ( distance < tolerance )
            {
                used[i] = true;
                a_change_was_made = true;
                // Check that the symmetry operator does not have an intrinsic translation--that would be weird
                if ( space_group_.symmetry_operator( i ).has_intrinsic_translation() )
                {
                    std::cout << "CrystalStructure::point_is_on_special_position( Vector3D, double ) : Warning: a symmetry operator with a non-zero intrinsic translation mapped an atom onto itself." << std::endl;
                }
                else
                {
                    // Average point and point_2
                    point += difference_vector / 2.0;
                }
            }
        }
    }
    std::vector< Matrix3D > point_group_symmetry_operators;
    point_group_symmetry_operators.push_back( Matrix3D() );
    for ( size_t i( 1 ); i != space_group_.nsymmetry_operators(); ++i )
    {
        Vector3D point_2 = space_group_.symmetry_operator( i ) * point;
        double distance;
        Vector3D difference_vector;
        crystal_lattice_.shortest_distance( point, point_2, distance, difference_vector );
        if ( distance < tolerance )
        {
            // Check that the symmetry operator does not have an intrinsic translation--that would be weird
            if ( space_group_.symmetry_operator( i ).has_intrinsic_translation() )
            {
                std::cout << "CrystalStructure::point_is_on_special_position( Vector3D, double ) : Warning: a symmetry operator with a non-zero intrinsic translation mapped an atom onto itself." << std::endl;
            }
            else
                point_group_symmetry_operators.push_back( space_group_.symmetry_operator( i ).rotation() );
        }
    }
    check_if_closed( point_group_symmetry_operators );
    PointGroup point_group( point_group_symmetry_operators );
    return point_group;
}

// ********************************************************************************

PointGroup CrystalStructure::point_is_on_special_position_const( Vector3D point, const double tolerance ) const
{
    return point_is_on_special_position( point, tolerance );
}

// ********************************************************************************

PointGroup CrystalStructure::molecule_is_on_special_position( const size_t i, const double tolerance ) const
{
    return point_is_on_special_position_const( molecule_in_crystal( i ).centre_of_mass(), tolerance );
}

// ********************************************************************************

//struct SpecialPositionsReport
//{
//      size_t number_of_atoms_in_unit_cell_; // Could make it a ChemicalFormula.
//    // We could enumerate the number of atoms found per point group, but two atoms that are on two special positions with the same point group are not necesarily the same atom.
//      std::vector< size_t > point_group_orders_;
//      std::vector< size_t > nmolecules_per_point_group_order_; // Same size as point_group_orders_
//      size_t nmolecules_on_special_positions_; // Sum of nmolecules_per_point_group_order_
//      size_t nsymmetry_operators_;
//    size_t nmolecules_;
//};

SpecialPositionsReport CrystalStructure::special_positions_report() const
{
    SpecialPositionsReport result;
 //   perceive_molecules();
    result.nsymmetry_operators_ = space_group_.nsymmetry_operators();
    //reduce_to_asymmetric_unit();
    //apply_space_group_symmetry();
    result.number_of_atoms_in_unit_cell_ = natoms();
    for ( size_t i( 0 ); i != nmolecules(); ++i )
    {
        PointGroup point_group = molecule_is_on_special_position( i );
        if ( point_group.nsymmetry_operators() != 1 )
        {
            bool found( false );
            for ( size_t j( 0 ); j != result.point_group_orders_.size(); ++j )
            {
                if ( result.point_group_orders_[j] == point_group.nsymmetry_operators() )
                {
                    ++result.nmolecules_per_point_group_order_[j];
                    found = true;
                }
            }
            if ( ! found )
            {
                result.point_group_orders_.push_back( point_group.nsymmetry_operators() );
                result.nmolecules_per_point_group_order_.push_back( 1 );
            }
        }
    }
    return result;
}

// ********************************************************************************

// @@ This is problematic because MoleculeInCrystal stores copies of all the atoms.
void CrystalStructure::move_molecule( const size_t i, const Vector3D shift )
{
    for ( size_t j( 0 ); j != molecules_[i].natoms(); ++j )
    {
        Atom new_atom( molecules_[i].atom( j ) );
        new_atom.set_position( new_atom.position() + shift );
        molecules_[i].set_atom( j, new_atom );
    }
}

// ********************************************************************************

void CrystalStructure::convert_to_P1()
{
    supercell( 1, 1, 1 );
}

// ********************************************************************************

void CrystalStructure::supercell( const size_t u, const size_t v, const size_t w )
{
    if ( ( u == 0 ) || ( v == 0 ) || ( w == 0 ) )
        throw std::runtime_error( "CrystalStructure::superstructure( u, v, w ): u, v, w cannot be 0." );
    if ( ! space_group_symmetry_has_been_applied() )
        apply_space_group_symmetry();
    CrystalStructure result;
    result.set_space_group( SpaceGroup() );
    result.set_name( name() );
    CrystalLattice new_crystal_lattice( crystal_lattice_.a() * u, crystal_lattice_.b() * v, crystal_lattice_.c() * w, crystal_lattice_.alpha(), crystal_lattice_.beta(), crystal_lattice_.gamma() );
    result.set_crystal_lattice( new_crystal_lattice );
    result.reserve_natoms( natoms() * u * v * w );
    for ( size_t i( 0 ); i != u; ++i )
    {
        for ( size_t j( 0 ); j != v; ++j )
        {
            for ( size_t k( 0 ); k != w; ++k )
            {
                for ( size_t l( 0 ); l != natoms(); ++l )
                {
                    Atom new_atom( atoms_[l] );
                    Vector3D new_position = atoms_[l].position(); // Fractional coordinates in the old unit cell
                    new_position = crystal_lattice_.fractional_to_orthogonal_matrix() * new_position; // Orthogonal coordinates (independent of unit cell)
                    new_position += i * crystal_lattice_.a_vector() + j * crystal_lattice_.b_vector() + k * crystal_lattice_.c_vector();
                    new_position = new_crystal_lattice.orthogonal_to_fractional( new_position ); // Fractional coordinates in the new unit cell
                    new_atom.set_position( new_position );
                    new_atom.set_label( atoms_[l].label() + "_" + size_t2string( i ) + "_" + size_t2string( j ) + "_" + size_t2string( k ) );
                    result.add_atom( new_atom );
                }
            }
        }
    }
    *this = result;
}

// ********************************************************************************

// Unit cell, atomic coordinates, ADPs and space group.
void CrystalStructure::transform( const Matrix3D & transformation_matrix )
{
    if ( ! nearly_equal( transformation_matrix.determinant(), 1.0 ) )
        std::cout << "Warning: CrystalStructure::transform(): the determinant of the transformation matrix is not 1." << std::endl;
    Matrix3D transformation_matrix_inverse_transpose( transformation_matrix );
    transformation_matrix_inverse_transpose.invert();
    transformation_matrix_inverse_transpose.transpose();
    CrystalLattice new_lattice( crystal_lattice_ );
    new_lattice.transform( transformation_matrix );
    for ( size_t i( 0 ); i != this->natoms(); ++i )
    {
        Atom new_atom( atoms_[i] );
        new_atom.set_position( transformation_matrix_inverse_transpose * new_atom.position() );
        if ( new_atom.ADPs_type() == Atom::ANISOTROPIC )
            new_atom.set_anisotropic_displacement_parameters( transform_adps( new_atom.anisotropic_displacement_parameters(), transformation_matrix, crystal_lattice_ ) );
        this->set_atom( i, new_atom );
    }
    space_group_.apply_similarity_transformation( SymmetryOperator( transformation_matrix_inverse_transpose, Vector3D() ) );
    crystal_lattice_ = new_lattice;
}

// ********************************************************************************

// Only atomic coordinates and ADPs.
// Takes output from find_match();
void CrystalStructure::transform( const SymmetryOperator & symmetry_operator, const std::vector< int > & integer_shifts )
{
        Matrix3D rotation = symmetry_operator.rotation();
        Vector3D shift = symmetry_operator.translation();
        shift.set_x( shift.x() + integer_shifts[0] );
        shift.set_y( shift.y() + integer_shifts[1] );
        shift.set_z( shift.z() + integer_shifts[2] );
        for ( size_t i( 0 ); i != this->natoms(); ++i )
        {
            Atom new_atom( this->atom( i ) );
            new_atom.set_position( ( rotation * new_atom.position() ) + shift );
            if ( new_atom.ADPs_type() == Atom::ANISOTROPIC )
                new_atom.set_anisotropic_displacement_parameters( rotate_adps( new_atom.anisotropic_displacement_parameters(), rotation, crystal_lattice_ ) );
            this->set_atom( i, new_atom );
        }
}

// ********************************************************************************

void CrystalStructure::position_all_atoms_within_unit_cell()
{
    for ( size_t i( 0 ); i != atoms_.size(); ++i )
        atoms_[ i ].set_position( adjust_for_translations( atoms_[ i ].position() ) );
}

// ********************************************************************************

Vector3D CrystalStructure::centre_of_mass() const
{
    if ( atoms_.empty() )
        throw std::runtime_error( "CrystalStructure::centre_of_mass(): there are no atoms, centre of mass is undefined." );
    bool ignore_hydrogens( false ); // Ignore deuteriums as well?
    bool weigh_by_atomic_weight( false );
    Vector3D result;
    double sum_of_weights( 0.0 );
    for ( size_t i( 0 ); i != atoms_.size(); ++i )
    {
        if ( ignore_hydrogens && ( atoms_[ i ].element() == Element( "H" ) ) )
            continue;
        double weight = ( weigh_by_atomic_weight ) ? atoms_[ i ].element().atomic_weight() : 1.0;
        result += weight * atoms_[ i ].position();
        sum_of_weights += weight;
    }
    result /= sum_of_weights;
    return result;
}

// ********************************************************************************

// @@ Undefined for a crystal structure, gives the wrong answer.
double CrystalStructure::dipole_moment() const
{
    // We must first get rid of the nett charge
    double nett_charge( 0.0 );
    for ( size_t i( 0 ); i != atoms_.size(); ++i )
        nett_charge += atoms_[ i ].charge();
    Vector3D centre_of_mass_negative_charge;
    Vector3D centre_of_mass_positive_charge;
    double sum_of_negative_charges( 0.0 );
    double sum_of_positive_charges( 0.0 );
    size_t natoms_zero_charge( 0 );
    Vector3D result;
    for ( size_t i( 0 ); i != atoms_.size(); ++i )
    {
        double charge = atoms_[ i ].charge() - ( nett_charge / atoms_.size() );
        if ( nearly_zero( charge ) )
            ++natoms_zero_charge;
        Vector3D position = atoms_[ i ].position();
//        position = adjust_for_translations( position ); // Move the atom to be within the unit cell
        position = crystal_lattice_.fractional_to_orthogonal( position );
        result += charge * position;
        if ( charge < 0.0 )
        {
            centre_of_mass_negative_charge += charge * position;
            sum_of_negative_charges += charge;
        }
        else
        {
            centre_of_mass_positive_charge += charge * position;
            sum_of_positive_charges += charge;
        }
    }
    centre_of_mass_negative_charge = -1.0*centre_of_mass_negative_charge;
    if ( natoms_zero_charge != 0 )
        std::cout <<  "CrystalStructure::dipole_moment(): Warning: " + size_t2string( natoms_zero_charge ) + " atoms have a charge equal 0.0." << std::endl;
    std::cout <<  "CrystalStructure::dipole_moment(): centre_of_mass_negative_charge: " << crystal_lattice_.orthogonal_to_fractional( centre_of_mass_negative_charge ) << std::endl;
    std::cout <<  "CrystalStructure::dipole_moment(): centre_of_mass_negative_charge: " << crystal_lattice_.orthogonal_to_fractional( centre_of_mass_negative_charge ) << std::endl;
    std::cout <<  "dipole moment: " << result << std::endl;
    return sum_of_positive_charges * ( centre_of_mass_negative_charge - centre_of_mass_positive_charge ).length();
}

// ********************************************************************************

double CrystalStructure::density() const
{
    if ( ! space_group_symmetry_has_been_applied() )
        std::cout << "CrystalStructure::density(): WARNING: space-group symmetry has not been applied, result will be nonsensical." << std::endl;
    ChemicalFormula chemical_formula;
    for ( size_t i( 0 ); i != atoms_.size(); ++i )
        chemical_formula.add_element( atoms_[ i ].element() );
    return ( chemical_formula.molecular_weight() / crystal_lattice_.volume() ) / ( Avogadros_constant / 1.0E24 );
}

// ********************************************************************************

void CrystalStructure::shortest_distance( const Vector3D & lhs, const Vector3D & rhs, double & shortest_distance, Vector3D & shortest_difference_vector )
{
    crystal_lattice_.shortest_distance( lhs, rhs, shortest_distance, shortest_difference_vector );
    // Loop over symmetry operators.
    for ( size_t k( 1 ); k != space_group_.nsymmetry_operators(); ++k )
    {
        // Fractional coordinates.
        Vector3D current_position = space_group_.symmetry_operator( k ) * rhs;
        // Adjust for translations and convert to Cartesian coordinates.
        double distance;
        Vector3D difference_vector; // Fractional coordinates.
        crystal_lattice_.shortest_distance( lhs, current_position, distance, difference_vector );
        if ( distance < shortest_distance )
        {
            shortest_distance = distance;
            shortest_difference_vector = difference_vector;
        }
    }
}

// ********************************************************************************

void CrystalStructure::second_shortest_distance( const Vector3D & lhs, const Vector3D & rhs, double & second_shortest_distance, Vector3D & second_shortest_difference_vector )
{
    double shortest_distance;
    crystal_lattice_.shortest_distance( lhs, rhs, shortest_distance, second_shortest_difference_vector );
    // Loop over symmetry operators.
    for ( size_t k( 1 ); k != space_group_.nsymmetry_operators(); ++k )
    {
        // Fractional coordinates.
        Vector3D current_position = space_group_.symmetry_operator( k ) * rhs;
        // Adjust for translations and convert to Cartesian coordinates.
        double distance;
        crystal_lattice_.shortest_distance( lhs, current_position, distance, second_shortest_difference_vector );
        if ( distance < shortest_distance )
            shortest_distance = distance;
    }
    second_shortest_distance = 1.0E12;
    // Loop over symmetry operators.
    for ( size_t k( 0 ); k != space_group_.nsymmetry_operators(); ++k )
    {
        // Fractional coordinates.
        Vector3D current_position = space_group_.symmetry_operator( k ) * rhs;
        // Adjust for translations and convert to Cartesian coordinates.
        double distance;
        Vector3D difference_vector; // Fractional coordinates.
        crystal_lattice_.shortest_distance( lhs, current_position, distance, difference_vector );
        if ( nearly_equal( distance, shortest_distance ) )
            continue;
        if ( distance < second_shortest_distance )
        {
            second_shortest_distance = distance;
            second_shortest_difference_vector = difference_vector;
        }
    }
}

// ********************************************************************************

// Finds shortest distance squared, in Angstrom squared, between two positions given in fractional coordinates.
// All space-group symmetry operators are taken into account; if this is undesired, use CrystalLattice::shortest_distance2().
double CrystalStructure::shortest_distance2( const Vector3D & lhs, const Vector3D & rhs )
{
    double shortest_distance2 = crystal_lattice_.shortest_distance2( lhs, rhs );
    // Loop over symmetry operators.
    for ( size_t k( 1 ); k != space_group_.nsymmetry_operators(); ++k )
    {
        // Fractional coordinates.
        Vector3D current_position = space_group_.symmetry_operator( k ) * rhs;
        // Adjust for translations and convert to Cartesian coordinates.
        double distance2 = crystal_lattice_.shortest_distance2( lhs, current_position );
        if ( distance2 < shortest_distance2 )
            shortest_distance2 = distance2;
    }
    return shortest_distance2;
}

// ********************************************************************************

// The current space group should be P1. u, v, w are the dimensions of the supercell with respect to
// the original unit cell, space_group is the space group of the original unit cell.
void CrystalStructure::collapse_supercell( const size_t u, const size_t v, const size_t w, const SpaceGroup & space_group )
{
    for ( size_t i( 0 ); i != atoms_.size(); ++i )
    {
        Vector3D position = atoms_[ i ].position();
        position.set_x( u * position.x() );
        position.set_y( v * position.y() );
        position.set_z( w * position.z() );
        atoms_[ i ].set_position( position );
    }
    position_all_atoms_within_unit_cell();
    CrystalLattice crystal_lattice( crystal_lattice_.a() / u,
                                    crystal_lattice_.b() / v,
                                    crystal_lattice_.c() / w,
                                    crystal_lattice_.alpha(),
                                    crystal_lattice_.beta(),
                                    crystal_lattice_.gamma() );
    crystal_lattice_ = crystal_lattice;
    // Now apply the symmetry operators (including unit-cell translations) to position each
    // atom as close as possible to 0,0,0 with all positive coordinates.
    for ( size_t i( 0 ); i != atoms_.size(); ++i )
    {
        Vector3D shortest_position = atoms_[ i ].position();
        double shortest_distance = shortest_position.length();
        for ( size_t j( 0 ); j != space_group.nsymmetry_operators(); ++j )
        {
            Vector3D new_position = space_group.symmetry_operator( j ) * atoms_[ i ].position();
            new_position = adjust_for_translations( new_position );
            if ( new_position.length() < shortest_distance )
            {
                shortest_position = new_position;
                shortest_distance = shortest_position.length();
            }
        }
        atoms_[ i ].set_position( shortest_position );
    }
}

// ********************************************************************************

// The current space group should be P1. lattice is the lattice of the original unit cell,
// from which the dimensions of the supercell are calculated.
// space_group is the space group of the original unit cell.
void CrystalStructure::collapse_supercell( const CrystalLattice & crystal_lattice, const SpaceGroup & space_group )
{
    size_t u = round_to_int( crystal_lattice_.a() / crystal_lattice.a() );
    size_t v = round_to_int( crystal_lattice_.b() / crystal_lattice.b() );
    size_t w = round_to_int( crystal_lattice_.c() / crystal_lattice.c() );
    collapse_supercell( u, v, w, space_group );
}

// ********************************************************************************

// The current space group should be P1. u, v, w are the dimensions of the supercell with respect to
// the original unit cell.
void CrystalStructure::collapse_supercell( const size_t u, const size_t v, const size_t w )
{
    for ( size_t i( 0 ); i != atoms_.size(); ++i )
    {
        Vector3D position = atoms_[ i ].position();
        position.set_x( u * position.x() );
        position.set_y( v * position.y() );
        position.set_z( w * position.z() );
        atoms_[ i ].set_position( position );
    }
    position_all_atoms_within_unit_cell();
    CrystalLattice crystal_lattice( crystal_lattice_.a() / u,
                                    crystal_lattice_.b() / v,
                                    crystal_lattice_.c() / w,
                                    crystal_lattice_.alpha(),
                                    crystal_lattice_.beta(),
                                    crystal_lattice_.gamma() );
    crystal_lattice_ = crystal_lattice;
    // Average the atomic coordinates.
    std::vector< Atom > new_atoms;
    size_t multiplicity = u * v * w;
    new_atoms.reserve( atoms_.size() / multiplicity );
    std::vector< bool > done( atoms_.size(), false );
    for ( size_t i( 0 ); i != atoms_.size(); ++i )
    {
        if ( done[ i ] )
            continue;
        RunningAverageAndESD< Vector3D > average_position;
        average_position.add_value( atoms_[ i ].position() );
        size_t natoms_for_average = 1;
        done[ i ] = true;
        for ( size_t j( i+1 ); j != atoms_.size(); ++j )
        {
            if ( done[ j ] )
                continue;
            double distance;
            Vector3D difference_vector;
            crystal_lattice_.shortest_distance( average_position.average(), atoms_[ j ].position(), distance, difference_vector );
            if ( distance < 0.3 )
            {
                if ( atoms_[ i ].element() != atoms_[ j ].element() )
                    std::cout << "CrystalStructure::collapse_supercell( ): Warning: the atoms to be averaged have different elements." << std::endl;
                ++natoms_for_average;
                average_position.add_value( average_position.average() + difference_vector );
                done[ j ] = true;
            }
        }
        if ( natoms_for_average != multiplicity )
            std::cout << "CrystalStructure::collapse_supercell( ): Warning: the number of averaged atoms (" + size_t2string(natoms_for_average) + ") is not equal to the multiplicity (" + size_t2string(multiplicity) + ")." << std::endl;
        new_atoms.push_back( Atom( atoms_[ i ].element(), average_position.average(), atoms_[ i ].label() ) );
    }
    atoms_ = new_atoms;
}

// ********************************************************************************

// The current space group should be P1. lattice is the lattice of the original unit cell,
// from which the dimensions of the supercell are calculated.
void CrystalStructure::collapse_supercell( const CrystalLattice & crystal_lattice )
{
    size_t u = round_to_int( crystal_lattice_.a() / crystal_lattice.a() );
    size_t v = round_to_int( crystal_lattice_.b() / crystal_lattice.b() );
    size_t w = round_to_int( crystal_lattice_.c() / crystal_lattice.c() );
    collapse_supercell( u, v, w );
}

// ********************************************************************************

// The current space group should be P1.
// Collapse supercell, assume order *in the unit cell* (not in the molecule) can be trusted
// (if there are n atoms in a unit cell, then atom n+1 corresponds to atom 1 in unit cell 1)
void CrystalStructure::collapse_supercell( const size_t u, const size_t v, const size_t w, const size_t natoms )
{
    for ( size_t i( 0 ); i != atoms_.size(); ++i )
    {
        Vector3D position = atoms_[ i ].position();
        position.set_x( u * position.x() );
        position.set_y( v * position.y() );
        position.set_z( w * position.z() );
        atoms_[ i ].set_position( position );
    }
    CrystalLattice crystal_lattice( crystal_lattice_.a() / u,
                                    crystal_lattice_.b() / v,
                                    crystal_lattice_.c() / w,
                                    crystal_lattice_.alpha(),
                                    crystal_lattice_.beta(),
                                    crystal_lattice_.gamma() );
    crystal_lattice_ = crystal_lattice;
    // Average the atomic coordinates
    size_t multiplicity = u * v * w;
    size_t natoms_per_unit_cell = atoms_.size() / multiplicity;
    std::vector< Atom > new_atoms;
    new_atoms.reserve( natoms_per_unit_cell );
    for ( size_t i( 0 ); i != natoms_per_unit_cell; ++i )
    {
        RunningAverageAndESD< Vector3D > average_position;
        average_position.add_value( atoms_[ i ].position() );
        for ( size_t j( 1 ); j != multiplicity; ++j )
        {
            size_t jatom = natoms_per_unit_cell * j + i;
            // Determine u, v and w for x, y and z.
            Vector3D jatom_position( atoms_[ jatom ].position() );
            int i_u = round_to_int( jatom_position.x() - atoms_[ i ].position().x() );
            int i_v = round_to_int( jatom_position.y() - atoms_[ i ].position().y() );
            int i_w = round_to_int( jatom_position.z() - atoms_[ i ].position().z() );
            jatom_position = Vector3D( jatom_position.x() - i_u, jatom_position.y() - i_v, jatom_position.z() - i_w );
            average_position.add_value( jatom_position );
            if ( atoms_[ i ].element() != atoms_[ jatom ].element() )
                std::cout << "CrystalStructure::collapse_supercell( ): Warning: the atoms to be averaged have different elements." << std::endl;
        }
        new_atoms.push_back( Atom( atoms_[ i ].element(), average_position.average(), atoms_[ i ].label() ) );
    }
    atoms_ = new_atoms;
}

// ********************************************************************************

// Collapse supercell, assume order *in the unit cell* (not in the molecule) can be trusted
// (if there are n atoms in a unit cell, then atom n+1 corresponds to atom 1 in unit cell 1)
void CrystalStructure::collapse_supercell( const size_t u,
                                           const size_t v,
                                           const size_t w,
                                           const int drift_correction,
                                           const Vector3D & target_centre,
                                           Matrix3D & transformation,
                                           Vector3D & actual_centre,
                                           std::vector< std::vector< Vector3D > > & positions )
{
    positions.clear();
    // The following loop is necessary but screws up the current crystal structure;
    // this method should essentially be const...
    // Correct for drift.
    actual_centre = Vector3D();
    for ( size_t i( 0 ); i != atoms_.size(); ++i )
        actual_centre += atoms_[ i ].position();
    actual_centre /= atoms_.size();
    if ( drift_correction != 0 )
    {
        for ( size_t i( 0 ); i != atoms_.size(); ++i )
        {
            Vector3D position = atoms_[ i ].position() - actual_centre + target_centre;
            atoms_[ i ].set_position( position );
        }
    }
    for ( size_t i( 0 ); i != atoms_.size(); ++i )
    {
        Vector3D position = atoms_[ i ].position();
        position.set_x( u * position.x() );
        position.set_y( v * position.y() );
        position.set_z( w * position.z() );
        atoms_[ i ].set_position( position );
    }
    CrystalLattice crystal_lattice( crystal_lattice_.a() / u,
                                    crystal_lattice_.b() / v,
                                    crystal_lattice_.c() / w,
                                    crystal_lattice_.alpha(),
                                    crystal_lattice_.beta(),
                                    crystal_lattice_.gamma() );
    crystal_lattice_ = crystal_lattice;
    size_t multiplicity = u * v * w * space_group_.nsymmetry_operators();
    size_t natoms_per_asymmetric_unit = atoms_.size() / multiplicity;
    positions.reserve( natoms_per_asymmetric_unit );
    size_t ndistances_gt_5( 0 );
    for ( size_t i( 0 ); i != natoms_per_asymmetric_unit; ++i )
    {
        Vector3D iatom_position( atoms_[ i ].position() );
        std::vector< Vector3D > positions_2;
        positions_2.push_back( iatom_position );
        for ( size_t j( 1 ); j != multiplicity; ++j )
        {
            size_t jatom = natoms_per_asymmetric_unit * j + i;
            if ( atoms_[ i ].element() != atoms_[ jatom ].element() )
                    std::cout << "CrystalStructure::collapse_supercell( ): Warning: the atoms to be averaged have different elements." << std::endl;
            double smallest_norm2 = 10000000.0;
            Vector3D smallest_norm2_position;
            for ( size_t k( 0 ); k != space_group_.nsymmetry_operators(); ++k ) // SpaceGroup::symmetry_operator( i ) returns a copy, so this is extremely wasteful.
            {
                Vector3D jatom_position = space_group_.symmetry_operator( k ) * atoms_[ jatom ].position();
                // Determine u, v and w for x, y and z.
                int i_u = round_to_int( jatom_position.x() - atoms_[ i ].position().x() );
                int i_v = round_to_int( jatom_position.y() - atoms_[ i ].position().y() );
                int i_w = round_to_int( jatom_position.z() - atoms_[ i ].position().z() );
                jatom_position = Vector3D( jatom_position.x() - i_u, jatom_position.y() - i_v, jatom_position.z() - i_w );
                Vector3D difference = jatom_position - iatom_position;
                difference = crystal_lattice_.fractional_to_orthogonal( difference );
                double norm2 = difference.norm2();
                if ( norm2 < smallest_norm2 )
                {
                    smallest_norm2 = norm2;
                    smallest_norm2_position = jatom_position;
                }
            }
            if ( smallest_norm2 > 25.0 )
                ++ndistances_gt_5;
            if ( smallest_norm2 == 10000000.0 )
                    std::cout << "Oops..." << std::endl;
            positions_2.push_back( smallest_norm2_position );
        }
        positions.push_back( positions_2 );
    }
    if ( ndistances_gt_5 > 0 )
        std::cout << "Number of distances > 5.0 A = " << size_t2string( ndistances_gt_5 ) << std::endl;
}

// ********************************************************************************

/*
160
WUBDOM
C      -1.031928   -0.865325    1.394524
C      -2.766437    0.404252    0.425293
C      -0.901004    1.591030    1.463277
C       0.838377    0.204561    2.425207
C       0.933200    2.644682    2.515255
O      -0.667888    3.912635    1.491873
C      -2.196361   -0.787397    0.723425
...
*/
void CrystalStructure::save_xyz( const FileName & file_name ) const
{
    TextFileWriter text_file_writer( file_name );
    text_file_writer.write_line( size_t2string( atoms_.size() ) );
    // Mercury gets confused when comment line empty (this is a bug in Mercury).
    if ( name_.empty() )
        text_file_writer.write_line( "Comment" );
    else
        text_file_writer.write_line( name_ );
    for ( size_t i( 0 ); i != atoms_.size(); ++i )
    {
        if ( suppressed_[i] )
            continue;
        Vector3D position = crystal_lattice_.fractional_to_orthogonal( atoms_[ i ].position() );
        text_file_writer.write_line( atoms_[ i ].element().symbol() + " " +
                                     double2string( position.x() ) + " " +
                                     double2string( position.y() ) + " " +
                                     double2string( position.z() ) );
    }
}

// ********************************************************************************

void CrystalStructure::save_cif( const FileName & file_name ) const
{
    TextFileWriter text_file_writer( file_name );
    text_file_writer.write_line( "data_" + name_ );
    if ( ! space_group_.name().empty() )
        text_file_writer.write_line( "_symmetry_space_group_name_H-M  '" + space_group_.name() + "'" );
//    text_file_writer.write_line( "_symmetry_Int_Tables_number     1" );
    text_file_writer.write_line( "_symmetry_cell_setting          " + LatticeSystem2string( crystal_lattice_.lattice_system() ) );
    text_file_writer.write_line( "_cell_length_a    " + double2string( crystal_lattice_.a(), 5 ) );
    text_file_writer.write_line( "_cell_length_b    " + double2string( crystal_lattice_.b(), 5 ) );
    text_file_writer.write_line( "_cell_length_c    " + double2string( crystal_lattice_.c(), 5 ) );
    text_file_writer.write_line( "_cell_angle_alpha " + double2string( crystal_lattice_.alpha().value_in_degrees(), 5 ) );
    text_file_writer.write_line( "_cell_angle_beta  " + double2string( crystal_lattice_.beta().value_in_degrees(), 5  ) );
    text_file_writer.write_line( "_cell_angle_gamma " + double2string( crystal_lattice_.gamma().value_in_degrees(), 5 ) );
    text_file_writer.write_line( "_cell_volume      " + double2string( crystal_lattice_.volume(), 5 ) );
    text_file_writer.write_line( "loop_" );
    text_file_writer.write_line( "_symmetry_equiv_pos_as_xyz" );
    for ( size_t i( 0 ); i != space_group_.nsymmetry_operators(); ++i )
        text_file_writer.write_line( space_group_.symmetry_operator( i ).to_string() );
    bool at_least_one_atom_has_anisotropic_ADPs( false );
    for ( size_t i( 0 ); i != atoms_.size(); ++i )
    {
        if ( suppressed_[i] )
            continue;
        if ( atoms_[i].ADPs_type() == Atom::ANISOTROPIC )
        {
            at_least_one_atom_has_anisotropic_ADPs = true;
            break;
        }
    }
    bool at_least_one_atom_has_isotropic_ADPs( at_least_one_atom_has_anisotropic_ADPs );
    if ( ! at_least_one_atom_has_anisotropic_ADPs )
    {
        for ( size_t i( 0 ); i != atoms_.size(); ++i )
        {
            if ( suppressed_[i] )
                continue;
            if ( atoms_[i].ADPs_type() == Atom::ISOTROPIC )
            {
                at_least_one_atom_has_isotropic_ADPs = true;
                break;
            }
        }
    }
    text_file_writer.write_line( "loop_" );
    text_file_writer.write_line( "_atom_site_label" ); // Needed for Materials Studio.
    text_file_writer.write_line( "_atom_site_type_symbol" );
    text_file_writer.write_line( "_atom_site_fract_x" );
    text_file_writer.write_line( "_atom_site_fract_y" );
    text_file_writer.write_line( "_atom_site_fract_z" );
    text_file_writer.write_line( "_atom_site_occupancy" );
    if ( at_least_one_atom_has_isotropic_ADPs )
        text_file_writer.write_line( "_atom_site_U_iso_or_equiv" );
    if ( at_least_one_atom_has_anisotropic_ADPs )
        text_file_writer.write_line( "_atom_site_adp_type" );
    size_t len( 2 );
    size_t current_size = 99;
    while ( atoms_.size() >= current_size )
    {
        ++len;
        current_size = 10 * current_size + 9;
    }
    for ( size_t i( 0 ); i != atoms_.size(); ++i )
    {
        if ( suppressed_[i] )
            continue;
        std::string label;
        if ( atoms_[i].label().empty() )
            label = atoms_[i].element().symbol() + size_t2string( i + 1, len, '0' );
        else
            label = atoms_[i].label();
        text_file_writer.write( label + " " +
                                     atoms_[i].element().symbol() + " " +
                                     double2string_pad_plus( atoms_[ i ].position().x(), 5, ' ' ) + " " +
                                     double2string_pad_plus( atoms_[ i ].position().y(), 5, ' ' ) + " " +
                                     double2string_pad_plus( atoms_[ i ].position().z(), 5, ' ' ) + " " +
                                     double2string_2( atoms_[ i ].occupancy(), 4 ) );
        if ( at_least_one_atom_has_isotropic_ADPs )
            text_file_writer.write( " " + double2string( atoms_[i].Uiso() ) );
        if ( at_least_one_atom_has_anisotropic_ADPs )
        {
            if ( atoms_[i].ADPs_type() == Atom::ANISOTROPIC )
                text_file_writer.write( " Uani" );
            else
                text_file_writer.write( " Uiso" );
        }
        text_file_writer.write_line();
    }
    if ( at_least_one_atom_has_anisotropic_ADPs )
    {
        text_file_writer.write_line( "loop_" );
        text_file_writer.write_line( "_atom_site_aniso_label" );
        text_file_writer.write_line( "_atom_site_aniso_U_11" );
        text_file_writer.write_line( "_atom_site_aniso_U_22" );
        text_file_writer.write_line( "_atom_site_aniso_U_33" );
        text_file_writer.write_line( "_atom_site_aniso_U_23" ); // U32 U13 U12 is the normal order in .cif files.
        text_file_writer.write_line( "_atom_site_aniso_U_13" );
        text_file_writer.write_line( "_atom_site_aniso_U_12" );
        for ( size_t i( 0 ); i != atoms_.size(); ++i )
        {
            if ( suppressed_[i] )
                continue;
            if ( atoms_[i].ADPs_type() == Atom::ANISOTROPIC )
            {
                std::string label;
                if ( atoms_[i].label().empty() )
                    label = atoms_[i].element().symbol() + size_t2string( i + 1, len, '0' );
                else
                    label = atoms_[i].label();
                SymmetricMatrix3D Ucif = atoms_[i].anisotropic_displacement_parameters().U_cif( crystal_lattice_ );
                text_file_writer.write_line( label + " " +
                                             double2string( Ucif.value( 0, 0 ) ) + " " +
                                             double2string( Ucif.value( 1, 1 ) ) + " " +
                                             double2string( Ucif.value( 2, 2 ) ) + " " +
                                             double2string( Ucif.value( 1, 2 ) ) + " " +
                                             double2string( Ucif.value( 0, 2 ) ) + " " +
                                             double2string( Ucif.value( 0, 1 ) ) );
            }
        }
    }
    text_file_writer.write_line( "#END" );
}

// ********************************************************************************

// | element | element | number of bonded atoms | number of bonded hydrogen/deuterium atoms | first bonded atom by atomic number, H == D, '0' if none | second bonded atom by atomic number
// | third bonded atom by atomic number | fourth bonded atom by atomic number | member of three-membered ring | member of four-membered ring
// | member of five-membered ring | member of six-membered ring | member of seven-membered ring | cyclic
// All these properties must be independent of the presence of 3D coordinates, they are topological attributes.
void CrystalStructure::calculate_topological_attributes()
{
    for ( size_t i( 0 ); i != natoms(); ++i )
    {
        Atom new_atom( atoms_[i] );
        std::string topological_attributes;
        topological_attributes = new_atom.element().symbol();
        topological_attributes = pad( topological_attributes, 2 );


        new_atom.set_topological_attributes( topological_attributes );
        this->set_atom( i, new_atom );
    }
}

// ********************************************************************************

double root_mean_square_Cartesian_displacement( const CrystalStructure & lhs, const CrystalStructure & rhs, const bool include_hydrogens )
{
    // Check if the number of atoms is the same.
    if ( lhs.natoms() != rhs.natoms() )
        throw std::runtime_error( "root_mean_square_Cartesian_displacement(): error: number of atoms is not the same." );
    // We assume that we are working only with the asymmetric unit here, so we have to check the space group,
    // whether or not space-group symmetry has been applied and the special positions.
    if ( lhs.space_group_symmetry_has_been_applied() != rhs.space_group_symmetry_has_been_applied() )
        throw std::runtime_error( "root_mean_square_Cartesian_displacement(): error: space-group symmetry should have been applied either to neither or to both." );
    if ( lhs.space_group_symmetry_has_been_applied() )
        std::cout << "root_mean_square_Cartesian_displacement(): warning: space-group symmetry has been applied to both." << std::endl;
    if ( ! same_symmetry_operators( lhs.space_group(), rhs.space_group() ) )
        std::cout << "root_mean_square_Cartesian_displacement(): warning: space-group symmetry operators are different." << std::endl;
    double result( 0.0 );
    // Atoms on special positions contribute fractionally.
    double nnon_H_atoms( 0.0 );
    for ( size_t i( 0 ); i != lhs.natoms(); ++i )
    {
        // If both H or D, skip this atom.
        if ( ( ! include_hydrogens ) && lhs.atom( i ).element().is_H_or_D() && rhs.atom( i ).element().is_H_or_D() )
            continue;
        // Check if elements the same.
        if ( lhs.atom( i ).element() != rhs.atom( i ).element() )
            throw std::runtime_error( "root_mean_square_Cartesian_displacement(): error: elements are not the same." );
        // Check that occupancies are the same
        if ( ! nearly_equal( lhs.atom( i ).occupancy(), rhs.atom( i ).occupancy() ) )
            throw std::runtime_error( "root_mean_square_Cartesian_displacement(): error: occupancies are not the same." );
        if ( ! nearly_equal( lhs.atom( i ).occupancy(), 1.0 ) )
            std::cout << "root_mean_square_Cartesian_displacement(): warning: occupancies equal, but not 1.0." << std::endl;
        double occupancy = lhs.atom( i ).occupancy();
        // No need to check special positions if space group has been applied.
        if ( ! lhs.space_group_symmetry_has_been_applied() )
        {
            Vector3D dummy = lhs.atom( i ).position();
            PointGroup point_group_lhs = lhs.point_is_on_special_position( dummy );
            dummy = rhs.atom( i ).position();
            PointGroup point_group_rhs = rhs.point_is_on_special_position( dummy );
            if ( ! same_symmetry_operators( point_group_lhs, point_group_rhs ) )
                throw std::runtime_error( "root_mean_square_Cartesian_displacement(): error: special positions are not the same." );
            if ( point_group_lhs.nsymmetry_operators() != 1 )
            {
                std::cout << "root_mean_square_Cartesian_displacement(): info: atom " + size_t2string( i ) + " on special position of order " + size_t2string( point_group_lhs.nsymmetry_operators() ) + "." << std::endl;
                occupancy /= static_cast<double>( point_group_lhs.nsymmetry_operators() );
            }
        }
        nnon_H_atoms += occupancy;
        double displacement;
        //  Cartesian displacement = (|G1 * r1 - G1 * r2| + |G2 * r1 - G2 * r2|) / 2.
        displacement = ( (lhs.crystal_lattice().fractional_to_orthogonal( lhs.atom( i ).position() ) - lhs.crystal_lattice().fractional_to_orthogonal( rhs.atom( i ).position() )).length() +
                         (rhs.crystal_lattice().fractional_to_orthogonal( lhs.atom( i ).position() ) - rhs.crystal_lattice().fractional_to_orthogonal( rhs.atom( i ).position() )).length() ) / 2.0;
        result += occupancy * square( displacement );
    }
    if ( nearly_zero( nnon_H_atoms ) )
        return 0.0;
    result /= nnon_H_atoms;
    result = sqrt( result );
    return result;
}

// ********************************************************************************

double RMSCD_with_matching( const CrystalStructure & lhs, const CrystalStructure & rhs, const size_t shift_steps, const bool include_hydrogens )
{
    // Some simple checks:
    size_t natoms = rhs.natoms();
    if ( natoms != lhs.natoms() )
        throw std::runtime_error( "RMSCD_with_matching(): numbers of atoms are not the same." );
    if ( natoms == 0 )
        return 0.0;
    CrystalLattice average_lattice = average( lhs.crystal_lattice(), rhs.crystal_lattice() );
    if ( ! nearly_equal( lhs.crystal_lattice(), rhs.crystal_lattice() ) )
        std::cout << "RMSCD_with_matching(): WARNING: unit cells differ." << std::endl;

    // Loop over all symmetry operators.
    // First find all floating axes; these are a problem if there is more than one residue in the asymmetric unit,
    // but we cannot detect that at the moment.

    std::vector< Vector3D > best_matches; // Fractional coordinates, from rhs.
    best_matches.reserve( natoms );
    std::vector< size_t > matching_indices; // Indices from rhs.
    matching_indices.reserve( natoms );
    // Find closest match
    std::vector< bool > done( natoms, false );
    // In principle, the two structures could have different space groups,
    // but for the moment they must have the same space group.
    if ( ! same_symmetry_operators( lhs.space_group(), rhs.space_group() ) )
        std::cout << "find_match(): WARNING: Space groups are different, this will give non-sensical results." << std::endl;
    SpaceGroup space_group = rhs.space_group();
    std::vector< Vector3D > shifts;
    if ( ( shift_steps == 0 ) || ( shift_steps == 1 ) )
        shifts.push_back( Vector3D( 0.0, 0.0, 0.0 ) );
    else
    {
        shifts.reserve( shift_steps*shift_steps*shift_steps );
        for ( size_t i1( 0 ); i1 != shift_steps; ++i1 )
        {
            for ( size_t i2( 0 ); i2 != shift_steps; ++i2 )
            {
                for ( size_t i3( 0 ); i3 != shift_steps; ++i3 )
                {
                    shifts.push_back( Vector3D( static_cast<double>(i1)/static_cast<double>(shift_steps), static_cast<double>(i2)/static_cast<double>(shift_steps), static_cast<double>(i3)/static_cast<double>(shift_steps) ) );
                }
            }
        }
    }
    for ( size_t i( 0 ); i != natoms; ++i )
    {
        double smallest_distance( 10000.0 );
        Vector3D best_match;
  //      size_t best_symmetry_operator;
  //      size_t best_shift;
        size_t matching_index = natoms + 1;
        // Loop over atoms.
        for ( size_t j( 0 ); j != natoms; ++j )
        {
            // Check if element the same.
            if ( lhs.atom( i ).element() != rhs.atom( j ).element() )
                continue;
            // Loop over symmetry operators.
            for ( size_t k( 0 ); k != space_group.nsymmetry_operators(); ++k )
            {
                // Loop over shifts.
                for ( size_t m( 0 ); m != shifts.size(); ++m )
                {
                    // Fractional coordinates.
                    Vector3D current_position = space_group.symmetry_operator( k ) * ( rhs.atom( j ).position() + shifts[m] );
                    // Adjust for translations and convert to Cartesian coordinates.
                    double shortest_distance;
                    Vector3D difference_vector; // Fractional coordinates.
                    average_lattice.shortest_distance( lhs.atom( i ).position(), current_position, shortest_distance, difference_vector );
                    if ( shortest_distance < smallest_distance )
                    {
                        smallest_distance = shortest_distance;
                        matching_index = j;
                        best_match = lhs.atom( i ).position() + difference_vector;
                //        best_symmetry_operator = k;
               //         best_shift = m;
                    }
                }
            }
        }
        std::cout << lhs.atom( i ).element().symbol() << " smallest distance = " << smallest_distance << std::endl;
        if ( done[ matching_index ] && ( ! lhs.atom( i ).element().is_H_or_D() ) )
        {
            std::cout << size_t2string( i ) << std::endl;
            throw std::runtime_error( "RMSCD_with_matching(): an atom has two matches." );
        }
        done[ matching_index ] = true;
        best_matches.push_back( best_match );
        matching_indices.push_back( matching_index );
    }
    // Save cif for this match.
    CrystalStructure reordered_crystal_structure;
    reordered_crystal_structure.set_crystal_lattice( rhs.crystal_lattice() );
    reordered_crystal_structure.set_space_group( rhs.space_group() );
    for ( size_t i( 0 ); i != natoms; ++i )
    {
        Atom new_atom( rhs.atom( matching_indices[i] ) );
        new_atom.set_position( best_matches[i] );
        reordered_crystal_structure.add_atom( new_atom );
    }
    reordered_crystal_structure.save_cif( FileName( "C:\\Users\\jacco\\Documents\\reordered.cif" ) );
    return root_mean_square_Cartesian_displacement( lhs, reordered_crystal_structure, include_hydrogens );
}

// ********************************************************************************

// Hydrogen / Deuterium is ignored
// To go from rhs to lhs, so rhs is changed and lhs is the target
SymmetryOperator find_match( const CrystalStructure & lhs, const CrystalStructure & rhs, const size_t shift_steps, std::vector< int > & integer_shifts, const bool add_inversion, const bool correct_floating_axes )
{
    // Some simple checks:
    size_t natoms = rhs.natoms();
    if ( natoms != lhs.natoms() )
        throw std::runtime_error( "find_match(): numbers of atoms are not the same." );
    CrystalLattice average_lattice = average( lhs.crystal_lattice(), rhs.crystal_lattice() );
    if ( ! nearly_equal( lhs.crystal_lattice(), rhs.crystal_lattice() ) )
        std::cout << "find_match(): WARNING: unit cells differ." << std::endl;
    if ( natoms == 0 )
        return SymmetryOperator();
    // Find closest match
    std::vector< bool > done( natoms, false );
    // In principle, the two structures could have different space groups,
    // but for the moment they must have the same space group.
    if ( ! same_symmetry_operators( lhs.space_group(), rhs.space_group() ) )
        std::cout << "find_match(): WARNING: Space groups are different, this will give non-sensical results." << std::endl;
    SpaceGroup space_group = rhs.space_group();
    // First find all floating axes; @@ these are a problem if there is more than one residue in the asymmetric unit,
    // but we cannot detect that at the moment.
    // @@ There is also a bug here for floating axes along a diagonal as found in cubic space groups.
    Vector3D floating_axes_correction;
    if ( correct_floating_axes )
    {
        Matrix3D sum( 0.0 );
        for ( size_t k( 0 ); k != space_group.nsymmetry_operators(); ++k )
            sum += space_group.symmetry_operator( k ).rotation();
        Vector3D com_lhs = lhs.centre_of_mass(); // Fractional coordinates
        Vector3D com_rhs = rhs.centre_of_mass();
        for ( size_t i( 0 ); i != 3; ++i )
        {
            if ( ! nearly_zero( sum.value( i, i ) ) )
            {
                std::cout << "Floating axis found " << i << std::endl;
                floating_axes_correction.set_value( i, com_lhs.value(i)-com_rhs.value(i));
            }
        }
    }
    if ( add_inversion && ( ! space_group.has_inversion_at_origin() ) )
        space_group.add_inversion_at_origin();
    // Add all combinations of shifts of 1/shift_steps along a, b and c.
    std::vector< Vector3D > shifts;
    if ( ( shift_steps == 0 ) || ( shift_steps == 1 ) )
        shifts.push_back( floating_axes_correction );
    else
    {
        shifts.reserve( shift_steps*shift_steps*shift_steps );
        for ( size_t i1( 0 ); i1 != shift_steps; ++i1 )
        {
            for ( size_t i2( 0 ); i2 != shift_steps; ++i2 )
            {
                for ( size_t i3( 0 ); i3 != shift_steps; ++i3 )
                {
                    shifts.push_back( floating_axes_correction + Vector3D( static_cast<double>(i1)/static_cast<double>(shift_steps), static_cast<double>(i2)/static_cast<double>(shift_steps), static_cast<double>(i3)/static_cast<double>(shift_steps) ) );
                }
            }
        }
    }
    std::vector< size_t > symmetry_operators_frequencies( space_group.nsymmetry_operators(), 0 );
    std::vector< size_t > shifts_frequencies( shifts.size(), 0 );
    std::vector< size_t > inversion_frequencies( 2, 0 );
    for ( size_t i( 0 ); i != natoms; ++i )
    {
        if ( lhs.atom( i ).element().is_H_or_D() )
            continue;
        double smallest_distance( 10000.0 );
        Vector3D best_match;
        size_t best_symmetry_operator( 0 ); // Stupid initialisation to silence compiler warnings
        size_t best_shift( 0 ); // Stupid initialisation to silence compiler warnings
        size_t matching_index = natoms + 1;
        // Loop over atoms.
        for ( size_t j( 0 ); j != natoms; ++j )
        {
            // Check if element the same.
            if ( lhs.atom( i ).element() != rhs.atom( j ).element() )
                continue;
            // Loop over symmetry operators.
            for ( size_t k( 0 ); k != space_group.nsymmetry_operators(); ++k )
            {
                // Loop over shifts.
                for ( size_t m( 0 ); m != shifts.size(); ++m )
                {
                    // Fractional coordinates.
                    Vector3D current_position = space_group.symmetry_operator( k ) * ( rhs.atom( j ).position() + shifts[m] );
                    // Adjust for translations and convert to Cartesian coordinates.
                    double shortest_distance;
                    Vector3D difference_vector; // Fractional coordinates.
                    average_lattice.shortest_distance( lhs.atom( i ).position(), current_position, shortest_distance, difference_vector );
                    if ( shortest_distance < smallest_distance )
                    {
                        smallest_distance = shortest_distance;
                        matching_index = j;
                        best_match = lhs.atom( i ).position() + difference_vector;
                        best_symmetry_operator = k;
                        best_shift = m;
                    }
                }
            }
        }
        std::cout << lhs.atom( i ).element().symbol() << " smallest distance = " << smallest_distance << std::endl;
        if ( done[ matching_index ] && ( ! lhs.atom( i ).element().is_H_or_D() ) )
        {
            std::cout << "find_match(): an atom has two matches." << std::endl;
          //  throw std::runtime_error( "find_match(): an atom has two matches." );
        }
        done[ matching_index ] = true;
        std::cout << "Symmetry operator = " << best_symmetry_operator << std::endl;
        std::cout << "Shift = " << best_shift << std::endl;
        ++symmetry_operators_frequencies[best_symmetry_operator];
        ++shifts_frequencies[best_shift];
    }
    // Check if there is one and only one symmetry operator that has the greatest frequency
    size_t most_common_space_group_index = symmetry_operators_frequencies.size();
    size_t most_common_space_group_frequency = 0;
    for ( size_t i( 0 ); i != symmetry_operators_frequencies.size(); ++i )
    {
        if ( symmetry_operators_frequencies[i] > most_common_space_group_frequency )
        {
            most_common_space_group_frequency = symmetry_operators_frequencies[i];
            most_common_space_group_index = i;
        }
    }
    // Check if any other symmetry operators are ever found
    for ( size_t i( 0 ); i != symmetry_operators_frequencies.size(); ++i )
    {
        if ( i == most_common_space_group_index )
            continue;
        if ( symmetry_operators_frequencies[i] != 0 )
            std::cout << "symmetry_operators_frequencies[" << i << "] = " << symmetry_operators_frequencies[i] << " " << space_group.symmetry_operator( i ).to_string() << std::endl;
    }
    // Check if there is one and only one shift that has the greatest frequency.
    size_t most_common_shift_index = shifts_frequencies.size();
    size_t most_common_shift_frequency = 0;
    for ( size_t i( 0 ); i != shifts_frequencies.size(); ++i )
    {
        if ( shifts_frequencies[i] > most_common_shift_frequency )
        {
            most_common_shift_frequency = shifts_frequencies[i];
            most_common_shift_index = i;
        }
    }
    // Check if any other shifts are ever found
    for ( size_t i( 0 ); i != shifts_frequencies.size(); ++i )
    {
        if ( i == most_common_shift_index )
            continue;
        if ( shifts_frequencies[i] != 0 )
            std::cout << "shifts_frequencies[" << i << "] = " << shifts_frequencies[i] << " " << shifts[i].to_string() << std::endl;
    }
    std::cout << "The most common symmetry operator, with " << most_common_space_group_frequency << " occurrences, is:" << std::endl;
    std::cout << space_group.symmetry_operator( most_common_space_group_index ).to_string() << std::endl;
    std::cout << "The most common shift, with " << most_common_shift_frequency << " occurrences, is:" << std::endl;
    shifts[ most_common_shift_index ].show();
    SymmetryOperator result( space_group.symmetry_operator( most_common_space_group_index ).rotation(), space_group.symmetry_operator( most_common_space_group_index ).rotation() * shifts[ most_common_shift_index ] + space_group.symmetry_operator( most_common_space_group_index ).translation() );
    // Apply found transformation to rhs
    Vector3D integer_shifts_2 = lhs.centre_of_mass() - ( result * rhs.centre_of_mass() ); // Fractional coordinates. lhs is the target.
    integer_shifts.clear();
    integer_shifts.push_back( round_to_int( integer_shifts_2.x() ) );
    integer_shifts.push_back( round_to_int( integer_shifts_2.y() ) );
    integer_shifts.push_back( round_to_int( integer_shifts_2.z() ) );
    std::cout << "Integer shifts are " << integer_shifts[0] << " " << integer_shifts[1] << " " << integer_shifts[2] << std::endl;
    std::cout << "Rotation = " << std::endl;
    space_group.symmetry_operator( most_common_space_group_index ).rotation().show();
    // @@ Can still be off by integer multiples!
    std::cout << "Translation = " << std::endl;
    result.translation().show();
    return result;
}

// ********************************************************************************

void CrystalStructure::match( const CrystalStructure & rhs, const size_t shift_steps, const bool allow_inversion, const bool correct_floating_axes )
{
    // Hydrogen / Deuterium is ignored
    // Some simple checks:
    size_t natoms = rhs.natoms();
    if ( natoms != this->natoms() )
        throw std::runtime_error( "CrystalStructure::match(): numbers of atoms are not the same." );        
    CrystalLattice average_lattice = average( crystal_lattice_, rhs.crystal_lattice() );
    if ( natoms == 0 )
        return;
    // Find closest match
    std::vector< bool > done( natoms, false );
    // In principle, the two structures could have different space groups,
    // but for the moment they must have the same space group.
    if ( ! same_symmetry_operators( space_group_, rhs.space_group() ) )
        std::cout << "CrystalStructure::match(): WARNING: Space groups are different, this will give non-sensical results." << std::endl;
    SpaceGroup space_group = space_group_;
    // First find all floating axes; @@ these are a problem if there is more than one residue in the asymmetric unit,
    // but we cannot detect that at the moment.
    // @@ There is also a bug here for floating axes along a diagonal as found in cubic space groups.
    Vector3D floating_axes_correction;
    if ( correct_floating_axes )
    {
        Matrix3D sum( 0.0 );
        for ( size_t k( 0 ); k != space_group.nsymmetry_operators(); ++k )
            sum += space_group.symmetry_operator( k ).rotation();
        Vector3D com_lhs = centre_of_mass(); // Fractional coordinates
        Vector3D com_rhs = rhs.centre_of_mass();
        for ( size_t i( 0 ); i != 3; ++i )
        {
            if ( ! nearly_zero( sum.value( i, i ) ) )
            {
                std::cout << "Floating axis found " << i << std::endl;
                floating_axes_correction.set_value( i, com_lhs.value(i)-com_rhs.value(i) );
            }
        }
    }
    if ( allow_inversion && ( ! space_group.has_inversion_at_origin() ) )
        space_group.add_inversion_at_origin();
    // Add all combinations of shifts of 1/shift_steps along a, b and c.
    std::vector< Vector3D > shifts;
    if ( ( shift_steps == 0 ) || ( shift_steps == 1 ) )
        shifts.push_back( floating_axes_correction );
    else
    {
        shifts.reserve( shift_steps*shift_steps*shift_steps );
        for ( size_t i1( 0 ); i1 != shift_steps; ++i1 )
        {
            for ( size_t i2( 0 ); i2 != shift_steps; ++i2 )
            {
                for ( size_t i3( 0 ); i3 != shift_steps; ++i3 )
                {
                    shifts.push_back( floating_axes_correction + Vector3D( static_cast<double>(i1)/static_cast<double>(shift_steps), static_cast<double>(i2)/static_cast<double>(shift_steps), static_cast<double>(i3)/static_cast<double>(shift_steps) ) );
                }
            }
        }
    }
    std::vector< size_t > symmetry_operators_frequencies( space_group.nsymmetry_operators(), 0 );
    std::vector< size_t > shifts_frequencies( shifts.size(), 0 );
    for ( size_t i( 0 ); i != natoms; ++i )
    {
        if ( atom( i ).element().is_H_or_D() )
            continue;
        double smallest_distance( 10000.0 );
        Vector3D best_match;
        size_t best_symmetry_operator( 0 ); // Stupid initialisation to silence compiler warnings
        size_t best_shift( 0 ); // Stupid initialisation to silence compiler warnings
        size_t matching_index = natoms + 1;
        // Loop over atoms.
        for ( size_t j( 0 ); j != natoms; ++j )
        {
            // Check if element the same.
            if ( atom( i ).element() != rhs.atom( j ).element() )
                continue;
            // Loop over symmetry operators.
            for ( size_t k( 0 ); k != space_group.nsymmetry_operators(); ++k )
            {
                // Loop over shifts.
                for ( size_t m( 0 ); m != shifts.size(); ++m )
                {
                    // Fractional coordinates.
                    Vector3D current_position = space_group.symmetry_operator( k ) * ( atom( i ).position() + shifts[m] );
                    // Adjust for translations and convert to Cartesian coordinates.
                    double shortest_distance;
                    Vector3D difference_vector; // Fractional coordinates.
                    average_lattice.shortest_distance( rhs.atom( j ).position(), current_position, shortest_distance, difference_vector );
                    if ( shortest_distance < smallest_distance )
                    {
                        smallest_distance = shortest_distance;
                        matching_index = j;
                        best_symmetry_operator = k;
                        best_shift = m;
                    }
                }
            }
        }
        std::cout << atom( i ).element().symbol() << " smallest distance = " << smallest_distance << std::endl;
        if ( done[ matching_index ] && ( ! atom( i ).element().is_H_or_D() ) )
        {
            std::cout << "CrystalStructure::match(): an atom has two matches." << std::endl;
          //  throw std::runtime_error( "find_match(): an atom has two matches." );
        }
        done[ matching_index ] = true;
        std::cout << "Symmetry operator = " << best_symmetry_operator << std::endl;
        std::cout << "Shift = " << best_shift << std::endl;
        ++symmetry_operators_frequencies[best_symmetry_operator];
        ++shifts_frequencies[best_shift];
    }
    // Check if there is one and only one symmetry operator that has the greatest frequency
    size_t most_common_space_group_index = symmetry_operators_frequencies.size();
    size_t most_common_space_group_frequency = 0;
    for ( size_t i( 0 ); i != symmetry_operators_frequencies.size(); ++i )
    {
        if ( symmetry_operators_frequencies[i] > most_common_space_group_frequency )
        {
            most_common_space_group_frequency = symmetry_operators_frequencies[i];
            most_common_space_group_index = i;
        }
    }
    // Check if any other symmetry operators are ever found.
    for ( size_t i( 0 ); i != symmetry_operators_frequencies.size(); ++i )
    {
        if ( i == most_common_space_group_index )
            continue;
        if ( symmetry_operators_frequencies[i] != 0 )
            std::cout << "symmetry_operators_frequencies[" << i << "] = " << symmetry_operators_frequencies[i] << " " << space_group.symmetry_operator( i ).to_string() << std::endl;
    }
    // Check if there is one and only one shift that has the greatest frequency.
    size_t most_common_shift_index = shifts_frequencies.size();
    size_t most_common_shift_frequency = 0;
    for ( size_t i( 0 ); i != shifts_frequencies.size(); ++i )
    {
        if ( shifts_frequencies[i] > most_common_shift_frequency )
        {
            most_common_shift_frequency = shifts_frequencies[i];
            most_common_shift_index = i;
        }
    }
    // Check if any other shifts are ever found.
    for ( size_t i( 0 ); i != shifts_frequencies.size(); ++i )
    {
        if ( i == most_common_shift_index )
            continue;
        if ( shifts_frequencies[i] != 0 )
            std::cout << "shifts_frequencies[" << i << "] = " << shifts_frequencies[i] << " " << shifts[i].to_string() << std::endl;
    }
    std::cout << "The most common symmetry operator, with " << most_common_space_group_frequency << " occurrences, is:" << std::endl;
    std::cout << space_group.symmetry_operator( most_common_space_group_index ).to_string() << std::endl;
    std::cout << "The most common shift, with " << most_common_shift_frequency << " occurrences, is:" << std::endl;
    shifts[ most_common_shift_index ].show();
    // Now apply the transformation to *this.
    rhs.centre_of_mass().show();
    centre_of_mass().show();
  //  com_diff.show();
    for ( size_t i( 0 ); i != natoms; ++i )
    {
        Atom new_atom( atom( i ) );
        new_atom.set_position( space_group.symmetry_operator( most_common_space_group_index ) * ( atom( i ).position() + shifts[most_common_shift_index] ) );
        if ( new_atom.ADPs_type() == Atom::ANISOTROPIC )
            new_atom.set_anisotropic_displacement_parameters( rotate_adps( new_atom.anisotropic_displacement_parameters(), space_group.symmetry_operator( most_common_space_group_index ).rotation(), crystal_lattice_ ) );
        set_atom( i, new_atom );
    }
    // We now only have the integer translations left.
    Vector3D com_diff = rhs.centre_of_mass() - centre_of_mass();
    rhs.centre_of_mass().show();
    centre_of_mass().show();
    com_diff.show();
    com_diff.set_x( round_to_int( com_diff.x() ) );
    com_diff.set_y( round_to_int( com_diff.y() ) );
    com_diff.set_z( round_to_int( com_diff.z() ) );
    for ( size_t i( 0 ); i != natoms; ++i )
    {
        Atom new_atom( atom( i ) );
        new_atom.set_position( atom( i ).position() + com_diff );
        set_atom( i, new_atom );
    }
}

// ********************************************************************************

// Hydrogen / Deuterium is ignored
void map( const CrystalStructure & to_be_changed, const CrystalStructure & target, const size_t shift_steps, Mapping & mapping, SymmetryOperator & symmetry_operator, std::vector< Vector3D > & translations, const bool allow_inversion, const bool correct_floating_axes )
{
    bool debug_output( false );
    // In principle, the two structures could have different space groups,
    // but for the moment they must have the same space group.
    if ( ! same_symmetry_operators( to_be_changed.space_group(), target.space_group() ) )
        std::cout << "map( CrystalStructure, CrystalStructure ): WARNING: Space groups are different, this will give non-sensical results." << std::endl;
    SpaceGroup space_group = to_be_changed.space_group();
    // Some simple checks:
    size_t natoms = target.natoms();
    if ( natoms != to_be_changed.natoms() )
        throw std::runtime_error( "map( CrystalStructure, CrystalStructure ): numbers of atoms are not the same." );
    CrystalLattice average_lattice = average( to_be_changed.crystal_lattice(), target.crystal_lattice() );
    // @@ The following is wrong because it fails if an H in one corresponds to a D in the other.
    if ( to_be_changed.chemical_formula().to_string() != target.chemical_formula().to_string() )
        throw std::runtime_error( "map( CrystalStructure, CrystalStructure ): chemical formulae are not the same." );
    if ( natoms == 0 )
        return;
    std::vector< size_t > mapping_temp;
 //   mapping.clear();
    translations.clear();
    // Find closest match
    std::vector< bool > done( natoms, false );
    // First find all floating axes; @@ these are a problem if there is more than one residue in the asymmetric unit,
    // but we cannot detect that at the moment.
    // @@ There is also a bug here for floating axes along a diagonal as found in cubic space groups.
    Vector3D floating_axes_correction;
    if ( correct_floating_axes )
    {
        Matrix3D sum( 0.0 );
        for ( size_t k( 0 ); k != space_group.nsymmetry_operators(); ++k )
            sum += space_group.symmetry_operator( k ).rotation();
        Vector3D com_lhs = to_be_changed.centre_of_mass(); // Fractional coordinates
        Vector3D com_rhs = target.centre_of_mass();
        for ( size_t i( 0 ); i != 3; ++i )
        {
            if ( ! nearly_zero( sum.value( i, i ) ) )
            {
                if ( debug_output )
                {
                    std::cout << "Floating axis found " << i << std::endl;
                    std::cout << "com_lhs.value(i) = " << com_lhs.value(i) << std::endl;
                    std::cout << "com_rhs.value(i) = " << com_rhs.value(i) << std::endl;
                }
                floating_axes_correction.set_value( i, com_rhs.value(i)-com_lhs.value(i) );
            }
        }
    }
    if ( allow_inversion && ( ! space_group.has_inversion_at_origin() ) )
        space_group.add_inversion_at_origin();
    // Add all combinations of shifts of 1/shift_steps along a, b and c.
    // If the space group contains centring vectors, do not add shifts that duplicate the centring vectors
    // I checked that this should not upset any floating origin corrections
    std::vector< Vector3D > centring_vectors = space_group.centring_vectors();
    std::vector< Vector3D > shifts;
    if ( ( shift_steps == 0 ) || ( shift_steps == 1 ) )
        shifts.push_back( floating_axes_correction );
    else
    {
        shifts.reserve( shift_steps*shift_steps*shift_steps );
        for ( size_t i1( 0 ); i1 != shift_steps; ++i1 )
        {
            for ( size_t i2( 0 ); i2 != shift_steps; ++i2 )
            {
                for ( size_t i3( 0 ); i3 != shift_steps; ++i3 )
                {
                    Vector3D shift( static_cast<double>(i1)/static_cast<double>(shift_steps), static_cast<double>(i2)/static_cast<double>(shift_steps), static_cast<double>(i3)/static_cast<double>(shift_steps) );
                    bool centring_vector_found( false );
                    for ( size_t k( 0 ); k != centring_vectors.size(); ++k )
                    {
                        if ( nearly_equal( shift, centring_vectors[k] ) )
                        {
                            centring_vector_found = true;
                            break;
                        }
                        
                    }
                    if ( ! centring_vector_found )
                        shifts.push_back( floating_axes_correction + shift );
                }
            }
        }
    }
    std::vector< size_t > symmetry_operators_frequencies( space_group.nsymmetry_operators(), 0 );
    std::vector< size_t > shifts_frequencies( shifts.size(), 0 );
    for ( size_t i( 0 ); i != natoms; ++i )
    {
        if ( to_be_changed.atom( i ).element().is_H_or_D() )
            continue;
        double smallest_distance( 10000.0 );
        Vector3D best_match;
        size_t best_symmetry_operator( 0 ); // Stupid initialisation to silence compiler warnings
        size_t best_shift( 0 ); // Stupid initialisation to silence compiler warnings
        size_t matching_index = natoms + 1;
        // Loop over atoms.
        for ( size_t j( 0 ); j != natoms; ++j )
        {
            // Check if element the same.
            if ( to_be_changed.atom( i ).element() != target.atom( j ).element() )
                continue;
            // Loop over symmetry operators.
            for ( size_t k( 0 ); k != space_group.nsymmetry_operators(); ++k )
            {
                // Loop over shifts.
                for ( size_t m( 0 ); m != shifts.size(); ++m )
                {
                    // Fractional coordinates.
                    Vector3D current_position = space_group.symmetry_operator( k ) * ( to_be_changed.atom( i ).position() + shifts[m] );
                    // Adjust for translations and convert to Cartesian coordinates.
                    double shortest_distance;
                    Vector3D difference_vector; // Fractional coordinates.
                    average_lattice.shortest_distance( target.atom( j ).position(), current_position, shortest_distance, difference_vector );
                    if ( shortest_distance < smallest_distance )
                    {
                        smallest_distance = shortest_distance;
                        matching_index = j;
                        best_symmetry_operator = k;
                        best_shift = m;
                    }
                }
            }
        }
        if ( debug_output )
            std::cout << to_be_changed.atom( i ).element().symbol() << " smallest distance = " << smallest_distance << std::endl;
        if ( matching_index == natoms + 1 ) // This is impossible at the moment because we check via the chemical formula that the same elements are found in both structures
            throw std::runtime_error( "map( CrystalStructure, CrystalStructure ): an atom could not be matched." );
        if ( done[ matching_index ] )
        {
            std::cout << "map( CrystalStructure, CrystalStructure ): an atom has two matches." << std::endl;
          //  throw std::runtime_error( "map( CrystalStructure, CrystalStructure ): an atom has two matches." );
        }
        done[ matching_index ] = true;
        ++symmetry_operators_frequencies[best_symmetry_operator];
        ++shifts_frequencies[best_shift];
        mapping_temp.push_back( matching_index );
    }
    // Check if there is one and only one symmetry operator that has the greatest frequency
    size_t most_common_symmetry_operator_index = symmetry_operators_frequencies.size();
    size_t most_common_symmetry_operator_frequency = 0;
    for ( size_t i( 0 ); i != symmetry_operators_frequencies.size(); ++i )
    {
        if ( symmetry_operators_frequencies[i] > most_common_symmetry_operator_frequency )
        {
            most_common_symmetry_operator_frequency = symmetry_operators_frequencies[i];
            most_common_symmetry_operator_index = i;
        }
    }
    symmetry_operator = space_group.symmetry_operator( most_common_symmetry_operator_index );
    // Check if any other symmetry operators are ever found.
    for ( size_t i( 0 ); i != symmetry_operators_frequencies.size(); ++i )
    {
        if ( i == most_common_symmetry_operator_index )
            continue;
        if ( symmetry_operators_frequencies[i] != 0 )
            std::cout << "symmetry_operators_frequencies[" << i << "] = " << symmetry_operators_frequencies[i] << " " << space_group.symmetry_operator( i ).to_string() << std::endl;
    }
    // Check if there is one and only one shift that has the greatest frequency.
    size_t most_common_shift_index = shifts_frequencies.size();
    size_t most_common_shift_frequency = 0;
    for ( size_t i( 0 ); i != shifts_frequencies.size(); ++i )
    {
        if ( shifts_frequencies[i] > most_common_shift_frequency )
        {
            most_common_shift_frequency = shifts_frequencies[i];
            most_common_shift_index = i;
        }
    }
    // Check if any other shifts are ever found.
    for ( size_t i( 0 ); i != shifts_frequencies.size(); ++i )
    {
        if ( i == most_common_shift_index )
            continue;
        if ( shifts_frequencies[i] != 0 )
            std::cout << "shifts_frequencies[" << i << "] = " << shifts_frequencies[i] << " " << shifts[i].to_string() << std::endl;
    }
    if ( debug_output )
    {
        std::cout << "The most common symmetry operator, with " << most_common_symmetry_operator_frequency << " occurrences, is:" << std::endl;
        std::cout << symmetry_operator.to_string() << std::endl;
        std::cout << "The most common shift, with " << most_common_shift_frequency << " occurrences, is:" << std::endl;
        shifts[ most_common_shift_index ].show();
    }
    mapping = Mapping( mapping_temp );
    for ( size_t i( 0 ); i != natoms; ++i )
    {
        // Fractional coordinates.
        Vector3D current_position = symmetry_operator * ( to_be_changed.atom( i ).position() + shifts[most_common_shift_index] );
        // Adjust for translations and convert to Cartesian coordinates.
        double shortest_distance;
        Vector3D difference_vector; // Fractional coordinates.
        average_lattice.shortest_distance( current_position, target.atom( mapping[ i ] ).position() , shortest_distance, difference_vector );
        // difference_vector is now the shortest distance (vector) with all integer translations factored out.
        // We also know the actual difference vector. The difference between the two is the integer translations.
        Vector3D integer_translations = target.atom( mapping[ i ] ).position() - difference_vector - current_position;
        translations.push_back( shifts[most_common_shift_index] + inverse( symmetry_operator.rotation() ) * integer_translations );
    }
    mapping.invert();
}

// ********************************************************************************

void CrystalStructure::apply_map( const Mapping & mapping, const SymmetryOperator & symmetry_operator, const std::vector< Vector3D > & translations )
{
    std::vector< Atom > new_atoms;
    new_atoms.reserve( natoms() );
    for ( size_t i( 0 ); i != natoms(); ++i )
    {
        Atom new_atom( atom( mapping[ i ] ) );
        new_atom.set_position( symmetry_operator * ( new_atom.position() + translations[ mapping[ i ] ] ) );
        if ( new_atom.ADPs_type() == Atom::ANISOTROPIC )
            new_atom.set_anisotropic_displacement_parameters( rotate_adps( new_atom.anisotropic_displacement_parameters(), symmetry_operator.rotation(), crystal_lattice_ ) );
        new_atoms.push_back( new_atom );
    }
    atoms_ = new_atoms;
}

// ********************************************************************************

