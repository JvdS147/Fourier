#ifndef CRYSTALSTRUCTURE_H
#define CRYSTALSTRUCTURE_H

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

class ChemicalFormula;
class FileName;
class Mapping;

#include "Atom.h"
#include "ConnectivityTable.h"
#include "CrystalLattice.h"
#include "MoleculeInCrystal.h"
#include "SpaceGroup.h"

#include <set>
#include <vector>

struct SpecialPositionsReport
{
    size_t number_of_atoms_in_unit_cell_; // Could make it a ChemicalFormula.
    // We could enumerate the number of atoms found per point group, but two atoms that are on two special positions with the same point group are not necesarily the same atom.
    std::vector< size_t > point_group_orders_;
    std::vector< size_t > nmolecules_per_point_group_order_; // Same size as point_group_orders_
    size_t nmolecules_on_special_positions_; // Sum of nmolecules_per_point_group_order_
    size_t nsymmetry_operators_;
    size_t nmolecules_;
};

//enum DriftCorrection { NONE, USE_FIRST_FRAME, USE_VECTOR };

/*
  A crystal structure

  Class is schizofrenic regarding what is stored: asymmetric unit or all atoms in the unit cell
  (or anything in between or even more atoms)

*/
class CrystalStructure
{
public:

    // Default constructor
    CrystalStructure();

    std::string name() const { return name_; }
    void set_name( const std::string & name ) { name_ = name; }

    size_t natoms() const { return atoms_.size(); }

    Atom atom( const size_t i ) const;

    // Returns natoms() when label not found.
    size_t find_label( const std::string & label ) const;

    // Label must be an exact match, i.e. "C11" does not match "C1" and "C1_0" does not match "C1".
    // Zero-based, throws if no match found.
    size_t atom( const std::string & atom_label ) const;

    // This is an expensive method...
    std::vector< Atom > atoms() const { return atoms_; }
    
    ChemicalFormula chemical_formula() const;

    void reserve_natoms( const size_t value ) { atoms_.reserve( value ); suppressed_.reserve( value ); }
    void add_atom( const Atom & atom );
    void add_atoms( const std::vector< Atom > & atoms );

    void remove_H_and_D();

    // Replaces an existing atom, enables making changes to atoms in the crystal
    // The alternative would have been to make atom(size_t) return a reference.
    // To keep all other attributes of the atom use something like:
    // Atom new_atom = crystal_structure.atom( i );
    // new_atom.set_X( X_new ); // Change property X to X_new
    // crystal_structure.set_atom( i, new_atom );
    void set_atom( const size_t i, const Atom & atom );

    // The following is a bit of a hack to avoid copying crystal structures all the time.
    // Individual atoms can be switched on and off. Switched off atoms are not saved.
    bool suppressed( const size_t i ) const { return suppressed_[i]; }
    void set_suppressed( const size_t i, const bool value ) { suppressed_[i] = value; }

    // Checks if there are overlapping atoms
    // Checks if any atom labels are duplicate
    void basic_checks() const;
    
    void make_atom_labels_unique();

    std::set< Element > elements() const;

    SpaceGroup space_group() const { return space_group_; }

    void set_space_group( const SpaceGroup & space_group ) { space_group_ = space_group; }

    CrystalLattice crystal_lattice() const { return crystal_lattice_; }

    // This should somehow be cross-checked with the space group, which suggests that the two should be combined into a class.
    void set_crystal_lattice( const CrystalLattice & crystal_lattice ) { crystal_lattice_ = crystal_lattice; }

    void list_all_bonds( std::vector< std::string > & labels_1, std::vector< std::string > & labels_2, std::vector< double > & bonds ) const;
    void list_all_angles( std::vector< std::string > & labels_1, std::vector< std::string > & labels_2, std::vector< std::string > & labels_3, std::vector< double > & angles ) const;

    // @@ The only correct way to interpret a CrystalStructure correctly is to use the series of commands:
    // reduce_to_asymmetric_unit();
    // apply_space_group_symmetry();
    // perceive_molecules();
    // remove_symmetry_related_molecules();
    //
    // And some intermediate results may actually be very useful

    // Only the asymmetric unit is kept, everything else is deleted.
    void reduce_to_asymmetric_unit( const double tolerance = 0.001 );

    bool space_group_symmetry_has_been_applied() const { return space_group_symmetry_has_been_applied_; }

    // For each atom, adds all symmetry-related atoms.
    // If the symmetry-related atom is less than 0.1 A away from the original atom (taking
    // periodicity into account), it is assumed to be on a special position and the symmetry-related
    // copy is discarded.
    // @@ How is this different from convert_to_P1()? (Space group is not reset...)
    // @@ How is this different from supercell( 1, 1, 1 ) ?
    // @@ If a molecule on a special position is present and it has been expanded and saved to cif by Mercury,
    // this method gives the wrong answer because.
    // @@ Should this be made such that it can only be applied once? (I.e. is ignored when called again after first call.)
    void apply_space_group_symmetry();

    // Moves atoms, currently only using integer translations, but symmetry operators
    // should also be included, so that the atoms form connected molecules.
    // This is important when dealing with floating axes.
    // If a molecule is on a special position, this function does NOT generate the second half.
    void move_atoms_to_form_molecules( const bool include_symmetry_operators );

    void perceive_molecules( const bool I_know_Zprime_is_one = false );

    // @@ This requires that you run the molecule preception method first
    void remove_symmetry_related_molecules();

    size_t nmolecules() const { return molecules_.size(); }

    // This returns a copy, so would copy all atoms.
    MoleculeInCrystal molecule_in_crystal( const size_t i ) const;

    void set_molecule_in_crystal( const size_t i, const MoleculeInCrystal & molecule_in_crystal ) { molecules_[i] = molecule_in_crystal; }

    // point is in fractional coordinates.
    // Tolerance is in Angstrom.
    // On return, point contains the point moved to the exact special position.
    PointGroup point_is_on_special_position( Vector3D & point, const double tolerance = 0.01 ) const;

    PointGroup point_is_on_special_position_const( Vector3D point, const double tolerance = 0.01 ) const;
    
    // @@ Can't be const because it calls perceive_molecules()
    PointGroup molecule_is_on_special_position( const size_t i, const double tolerance = 0.01 ) const;

    // Must have called perceive_molecules() first.
    SpecialPositionsReport special_positions_report() const;
    
    // @@ This is problematic because MoleculeInCrystal stores copies of all the atoms.
    void move_molecule( const size_t i, const Vector3D shift );

    void convert_to_P1();

    // @@ Atom labels are currently not updated. But uniqueness of atom labels is not enforced anyway (probably should be)
    void supercell( const size_t u, const size_t v, const size_t w );

    // Unit cell, atomic coordinates, ADPs and space group.
    void transform( const Matrix3D & transformation_matrix );

    // Only atomic coordinates and ADPs.
    // Takes output from find_match();
    void transform( const SymmetryOperator & symmetry_operator, const std::vector< int > & integer_shifts );

    // Atoms with coordinates like 0.999999: keep at 0.999999 or move to 0.0 or move to -0.000001?
    void position_all_atoms_within_unit_cell();

    // Fractional coordinates.
    Vector3D centre_of_mass() const;

    // Unit: eA
    // @@ Undefined for a crystal structure, gives the wrong answer.
    double dipole_moment() const;

    // g/cm3
    // Must call CrystalStructure::apply_space_group_symmetry() first.
    double density() const;

    // Finds shortest distance, in Angstrom, between two positions given in fractional coordinates.
    // All space-group symmetry operators are taken into account; if this is undesired, use CrystalLattice::shortest_distance().
    // Returns the shortest distance (in Angstrom) and the shortest difference vector (defined as rhs - lhs, in fractional coordinates).
    void shortest_distance( const Vector3D & lhs, const Vector3D & rhs, double & distance, Vector3D & difference_vector );

    void second_shortest_distance( const Vector3D & lhs, const Vector3D & rhs, double & second_shortest_distance, Vector3D & second_shortest_difference_vector );

    // Finds shortest distance squared, in Angstrom squared, between two positions given in fractional coordinates.
    // All space-group symmetry operators are taken into account; if this is undesired, use CrystalLattice::shortest_distance2().
    double shortest_distance2( const Vector3D & lhs, const Vector3D & rhs );

    // The current space group should be P1. u, v, w are the dimensions of the supercell with respect to
    // the original unit cell, space_group is the space group of the original unit cell.
    void collapse_supercell( const size_t u, const size_t v, const size_t w, const SpaceGroup & space_group );

    // The current space group should be P1. lattice is the lattice of the original unit cell,
    // from which the dimensions of the supercell are calculated.
    // space_group is the space group of the original unit cell.
    void collapse_supercell( const CrystalLattice & crystal_lattice, const SpaceGroup & space_group );

    // The current space group should be P1. u, v, w are the dimensions of the supercell with respect to
    // the original unit cell.
    void collapse_supercell( const size_t u, const size_t v, const size_t w );

    // The current space group should be P1. lattice is the lattice of the original unit cell,
    // from which the dimensions of the supercell are calculated.
    void collapse_supercell( const CrystalLattice & crystal_lattice );

// @@ natoms is currently not used but is used to disambiguate the overload...

    // The current space group should be P1. u, v, w are the dimensions of the supercell with respect to
    // the original unit cell.
    // Collapse supercell, assume order *in the unit cell* (not in the molecule) can be trusted
    // (if there are n atoms in a unit cell, then atom n+1 corresponds to atom 1 in unit cell 1)
    void collapse_supercell( const size_t u, const size_t v, const size_t w, const size_t natoms );

    // This function is mainly for use in the AnalyseTrajectory class. You probably do not need to call it directly but you should
    // use an AnalyseTrajectory object instead.
    // u, v, w are the dimensions of the supercell with respect to the original unit cell.
    // Collapse supercell, assume order *in the unit cell* (not in the molecule) can be trusted
    // (if there are n atoms in a unit cell, then atom n+1 corresponds to atom 1 in unit cell 1).
    // The method assumes that you have not repositioned atoms to lie within the unit cell.
    // The method assumes that we are dealing with a solid where atomic coordinates are fairly constant and symmetry-related copies can be
    // easily identified from simple geometric considerations.
    // Space-group symmetry is also used to group *all* symmetry copies of each atom--so make sure that the space group is correct,
    // i.e. that no phase transition has taken place. Of course, if the space group is P1 there is no problem.
    // Setting the space group screws up the structure! This method should essentially be const.
    // positions contains fractional coordinates with respect to one unit cell (so after collapsing the supercell).
    // You'll get horrible results when applying the space-group symmetry if you repositioned the molecules so that all molecules were comfortably within the unit cell.
    // transformation does not work
    void collapse_supercell( const size_t u,
                             const size_t v,
                             const size_t w,
                             const int drift_correction,
                             const Vector3D & target_centre,
                             Matrix3D & transformation,
                             Vector3D & actual_centre,
                             std::vector< std::vector< Vector3D > > & positions );

    void save_xyz( const FileName & file_name ) const;
    
    void save_cif( const FileName & file_name ) const;

    void calculate_topological_attributes();

    // Change *this to correspond to the other crystal structure (supplied as argument).
    void match( const CrystalStructure & rhs, const size_t shift_steps, const bool allow_inversion, const bool correct_floating_axes );

    // The matching is becoming a bit of a mess.
    // Applies result of map() to *this.
    void apply_map( const Mapping & mapping, const SymmetryOperator & symmetry_operator, const std::vector< Vector3D > & translations );

private:
    SpaceGroup space_group_;
    CrystalLattice crystal_lattice_;
    std::vector< Atom > atoms_;
    std::vector< MoleculeInCrystal > molecules_;
    std::vector< bool > suppressed_; //
    std::string name_;
    bool space_group_symmetry_has_been_applied_;
    ConnectivityTable connectivity_table_;
};

// Atoms on special positions contribute fractionally if space group has not been applied.
// If there are atoms on special positions and the space group has not been applied then corresponding atoms must be on the same special position.
// Atoms with partial occupancies contribute fractionally if the occupancies are the same in both structures.
// This function in principle operates on the asymmetric unit, so space group should not have been applied.
// If space group has been applied, it should have been applied to both.
double root_mean_square_Cartesian_displacement( const CrystalStructure & lhs, const CrystalStructure & rhs, const bool include_hydrogens );

// This really should not be used any more, use find_match() followed by root_mean_square_Cartesian_displacement();
double RMSCD_with_matching( const CrystalStructure & lhs, const CrystalStructure & rhs, const size_t shift_steps, const bool include_hydrogens );

// Hydrogen / Deuterium is ignored
// What is returned is not necessarily a symmetry operator, it is a combination of a rotation matrix and a translation vector that may or may not correspond to a symmetry operator.
// Note that symmetry operators are canonicalised and the translational part must always be [0,1>. integer_shifts carries the remainder and is a vector of dimension three in fractional coordinates.
// To go from rhs to lhs, so rhs is changed and lhs is the target
// @@ Should also incorporate transformation of the unit cell
// @@ Should probably be a class
// @@ currently the result is returned in two variables, perhaps better to simply return the matching crystal structure.
SymmetryOperator find_match( const CrystalStructure & lhs, const CrystalStructure & rhs, const size_t shift_steps, std::vector< int > & integer_shifts, const bool add_inversion, const bool correct_floating_axes );

// Hydrogen and deuterium are ignored
// Maybe this should be a class
void map( const CrystalStructure & to_be_changed, const CrystalStructure & target, const size_t shift_steps, Mapping & mapping, SymmetryOperator & symmetry_operator, std::vector< Vector3D > & translations, const bool allow_inversion, const bool correct_floating_axes );

#endif // CRYSTALSTRUCTURE_H

