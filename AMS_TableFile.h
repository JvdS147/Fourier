#ifndef AMS_TABLEFILE_H
#define AMS_TABLEFILE_H

#include "CrystalLattice.h"
#include "FileName.h"
#include "Fraction.h"

#include <string>
#include <vector>

struct AMS_TableFileEntry
{

    double energy() const { return energy_; }
    double density() const { return density_; }
    double volume() const { return volume_; }
    double natoms() const { return natoms_; }
    std::string space_group() const { return space_group_; }
    CrystalLattice crystal_lattice() const { return crystal_lattice_; }

    bool operator< ( const AMS_TableFileEntry & rhs ) const { return ( this->energy() < rhs.energy() ); }

    double energy_;
    double density_;
    double volume_;
    double natoms_;
    std::string space_group_;
    CrystalLattice crystal_lattice_;

};

/*
  A table file

  Contents can be added (to enable merging two files) but not changed.
  
*/
class AMS_TableFile
{
public:

    AMS_TableFile();

    explicit AMS_TableFile( const FileName & file_name, const size_t natoms_per_molecule );

    void read_file( const FileName & file_name );

    size_t size() const { return sorted_map_.size(); }

    void add_entry( const AMS_TableFileEntry & entry, const FileName & file_name );

    AMS_TableFileEntry entry( const size_t i ) const { return entries_[sorted_map_[i]]; }
    
    size_t natoms_per_molecule() const { return natoms_per_molecule_; }
    
    void set_natoms_per_molecule( const size_t natoms_per_molecule ) { natoms_per_molecule_ = natoms_per_molecule; }

    FileName file_name( const size_t i ) const { return file_names_[sorted_map_[i]]; }

    // For consistency with the AMS_RESULTS_File interface
    double lowest_energy() const { return entries_[sorted_map_[0]].energy(); }
    
    // Energy in kcal/mol. Returns the number of predicted structures in the bottom "energy" kcal/mol
    size_t npolymorphs_in_bottom( const double energy ) const;
    
    void debug() const;
    
    // ############## The user must create the structures directory
    void save( const FileName & file_name, const std::string & structures_directory ) const;

private:
    std::vector< AMS_TableFileEntry > entries_;
    std::vector< FileName > file_names_;
    size_t natoms_per_molecule_;
    std::vector< Fraction > Z_primes_;
    // We don't actually sort the list, but create a sorted map
    std::vector< size_t > sorted_map_;
    
    void sort_by_energy();

};

AMS_TableFile merge( const AMS_TableFile & lhs, const AMS_TableFile & rhs );

#endif // AMS_TABLEFILE_H

