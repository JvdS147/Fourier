#include "AMS_TableFile.h"
#include "CopyTextFile.h"
#include "FileList.h"
#include "FileName.h"
#include "Sort.h"
#include "TextFileReader_2.h"
#include "TextFileWriter.h"
#include "Utilities.h"

#include <iostream>
#include <stdexcept>

// ********************************************************************************

AMS_TableFile::AMS_TableFile():
natoms_per_molecule_(0)
{
}

// ********************************************************************************

AMS_TableFile::AMS_TableFile( const FileName & file_name, const size_t natoms_per_molecule ):
natoms_per_molecule_(natoms_per_molecule)
{
    read_file( file_name );
}

// ********************************************************************************

// rank |       energy        |      density       |       volume       | atoms_per_asym_unit |  space_group  |         a          |         b          |         c          |       alpha        |        beta        |       gamma        |  ID
//      |   [kcal/mol/atom]   |      [g/cm3]       |        [A3]        |                     |               |        [A]         |        [A]         |        [A]         |     [degrees]      |     [degrees]      |     [degrees]      |
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//  1   | -147.27727963570999 | 1.5432619149167712 | 2250.1702278782445 |         101         |     P_2_1     | 10.097744266882842 | 21.480522266328531 | 10.606498849034837 |         90         | 102.01877453727563 |         90         | 32171
//  2   | -147.27644561860046 | 1.5517080301831858 | 2237.9223070426024 |         101         |     P_2_1     | 10.145541065570267 | 21.519383596782387 | 10.492298936172645 |         90         | 102.32750692394498 |         90         | 32172
//  3   |  -147.27572122395   | 1.5495514850386483 | 2241.0368731165295 |         101         |     P_2_1     | 10.236961973858019 | 21.524559306512003 | 10.43701142963701  |         90         | 77.02491219320639  |         90         | 32173
//  0             1                     2                   3                      4                  5                  6                    7                   8                     9                    10                  11             12
void AMS_TableFile::read_file( const FileName & file_name )
{
    TextFileReader_2 text_file_reader( file_name );
    if ( text_file_reader.size() < 4 )
        throw std::runtime_error( "AMS_TableFile::read_file(): file contains less than four lines." );
    Splitter splitter( "|" );
    bool Zprime_1_found( false );
    std::vector< std::string > words = splitter.split( text_file_reader.line( 3 ) );
    for ( size_t i( 0 ); i != words.size(); ++i )
        words[i] = remove( words[i], ' ' );
    double lowest_energy = string2double( words[1] );
    size_t natoms_0 = string2integer( words[4] );
    size_t natoms_1( 0 );
    std::vector< Fraction > natoms_fractions;
    bool Z_prime_consistent( true );
    entries_.reserve( text_file_reader.size() - 3 );
    for ( size_t iLine( 3 ); iLine != text_file_reader.size(); ++iLine )
    {
        words = splitter.split( text_file_reader.line( iLine ) );
        for ( size_t i( 0 ); i != words.size(); ++i )
            words[i] = remove( words[i], ' ' );
        AMS_TableFileEntry table_file_entry;
        table_file_entry.energy_ = string2double( words[1] );
        table_file_entry.density_ = string2double( words[2] );
        table_file_entry.volume_ = string2double( words[3] );
        table_file_entry.natoms_ = string2integer( words[4] );
        table_file_entry.space_group_ = words[5];
        table_file_entry.crystal_lattice_ = CrystalLattice( string2double( words[6] ),
                                                            string2double( words[7] ),
                                                            string2double( words[8] ),
                                                            Angle::from_degrees( string2double( words[9] ) ),
                                                            Angle::from_degrees( string2double( words[10] ) ) ,
                                                            Angle::from_degrees( string2double( words[11] ) ) );
        entries_.push_back( table_file_entry );
        // The asymmetric unit can be e.g. 50.5 (for e.g. a hemihydrate with the H2O on a two-fold axis).
        double natoms_double = string2double( words[4] );
        Fraction natoms_fraction = double2fraction( natoms_double, Fraction( 1, 24 ) );
        natoms_fractions.push_back( natoms_fraction );
        if ( natoms_fraction.is_integer() )
        {
            size_t natoms = natoms_fraction.integer_part();
            if ( natoms != natoms_0 )
            {
                if ( natoms_1 == 0 )
                {
                    natoms_1 = natoms;
                    if ( natoms_0 == 2 * natoms_1 )
                        std::swap( natoms_0, natoms_1 );
                    else if ( natoms_1 != 2 * natoms_0 )
                    {
                        std::cout << "Warning: Number of atoms in asymmetric unit does not seem to be Z'=1 or Z'=2." << std::endl;
                        Z_prime_consistent = false;
                    }
                }
                else if ( natoms != natoms_1 )
                {
                    std::cout << "Warning: Number of atoms in asymmetric unit does not seem to be Z'=1 or Z'=2." << std::endl;
                    Z_prime_consistent = false;
                }
            }
            if ( natoms_per_molecule_ != 0 )
            {
                if ( natoms == natoms_per_molecule_ )
                    Zprime_1_found = true;
                else if ( natoms != ( 2 * natoms_per_molecule_ ) )
                {
                    std::cout << "Warning: Number of atoms in asymmetric unit incompatible with specified number of atoms per molecule." << std::endl;
                    Z_prime_consistent = false;
                }

            }
        }
    }
    if ( natoms_per_molecule_ != 0 )
    {
        if ( ! Zprime_1_found )
        {
            std::cout << "Warning: no Z'=1 structures found." << std::endl;
            Z_prime_consistent = false;
        }
    }
    for ( size_t i( 0 ); i != entries_.size(); ++i )
        file_names_.push_back( FileName( file_name.directory() + "structures", "structure_" + size_t2string( i+1, 6, '0' ), "cif" ) );
    sort_by_energy();
}

// ********************************************************************************

void AMS_TableFile::add_entry( const AMS_TableFileEntry & entry, const FileName & file_name )
{
    entries_.push_back( entry );
    file_names_.push_back( file_name );
    sort_by_energy();
}

// ********************************************************************************

size_t AMS_TableFile::npolymorphs_in_bottom( const double energy ) const
{
    size_t i( 0 );
    while ( ( (i) < size() ) && ( ( ( entry(i).energy() - lowest_energy() ) * natoms_per_molecule() ) < energy ) )
        ++i;
    return i;
}

// ********************************************************************************

void AMS_TableFile::save( const FileName & file_name, const std::string & structures_directory ) const
{
    std::string structures_directory_2 = append_backslash( structures_directory );
    TextFileWriter text_file_writer( file_name );
    text_file_writer.write_line( " rank |       energy        |      density       |       volume       | atoms_per_asym_unit |  space_group  |         a          |         b          |         c          |       alpha        |        beta        |       gamma        |  ID" );
    text_file_writer.write_line( "      |   [kcal/mol/atom]   |      [g/cm3]       |        [A3]        |                     |               |        [A]         |        [A]         |        [A]         |     [degrees]      |     [degrees]      |     [degrees]      |" );
    text_file_writer.write_line( "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" );
//  1   | -147.27727963570999 | 1.5432619149167712 | 2250.1702278782445 |         101         |     P_2_1     | 10.097744266882842 | 21.480522266328531 | 10.606498849034837 |         90         | 102.01877453727563 |         90         | 32171

    FileName file_list_name( structures_directory_2, "FileList", "txt" );
    TextFileWriter file_list_writer( file_list_name );
    for ( size_t i( 0 ); i != this->size(); ++i )
    {
        text_file_writer.write_line( centre( size_t2string( i+1 ), 5 ) + " | " +
                                     double2string( entry(i).energy(), 14, 19 ) + " | " +
                                     double2string( entry(i).density(), 16, 18 ) + " | " +
                                     double2string( entry(i).volume(), 13, 18 ) + " | " +
                                     centre( size_t2string( entry(i).natoms() ), 19 ) + " | " +
                                     centre( entry(i).space_group(), 13 ) + " | " +
                                     double2string( entry(i).crystal_lattice().a(), 15, 18 ) + " | " +
                                     double2string( entry(i).crystal_lattice().b(), 15, 18 ) + " | " +
                                     double2string( entry(i).crystal_lattice().c(), 15, 18 ) + " | " +
                                     double2string( entry(i).crystal_lattice().alpha().value_in_degrees(), 14, 18 ) + " | " +
                                     double2string( entry(i).crystal_lattice().beta().value_in_degrees(), 14, 18 ) + " | " +
                                     double2string( entry(i).crystal_lattice().gamma().value_in_degrees(), 14, 18 ) + " | " +
                                     centre( size_t2string( i+1 ), 5 ) );
        copy_text_file( this->file_name( i ), FileName( structures_directory_2, "structure_" + size_t2string( i+1, 6, '0' ), "cif" ) );
        file_list_writer.write_line( "structure_" + size_t2string( i+1, 6, '0' ) + ".cif" );
    }
    file_list_writer.~TextFileWriter();
    // Insert the correct data name into the cif files
    FileList file_list( file_list_name );
    file_list.set_prepend_file_name_with_basedirectory( true );
    for ( size_t i( 0 ); i != file_list.size(); ++i )
    {
        TextFileReader_2 input_file( file_list.value( i ) );
        TextFileWriter output_file( file_list.value( i ) );
        for ( size_t iLine( 0 ); iLine != input_file.size(); ++iLine )
        {
            if ( input_file.line( iLine ).substr( 0, 5 ) == "data_" )
                output_file.write_line( "data_" + size_t2string( i + 1, 6, '0' ) );
            else
                output_file.write_line( input_file.line( iLine ) );
        }
    }
    
}

// ********************************************************************************

void AMS_TableFile::sort_by_energy()
{
// We don't actually sort the lists, but create a sorted map
    sorted_map_ = sort( entries_ );
}

// ********************************************************************************

AMS_TableFile merge( const AMS_TableFile & lhs, const AMS_TableFile & rhs )
{
    if ( lhs.natoms_per_molecule() != rhs.natoms_per_molecule() )
        throw std::runtime_error( "merge(): ERROR number of atoms differ." );
    AMS_TableFile result( lhs );
    for ( size_t i( 0 ); i != rhs.size(); ++i )
        result.add_entry( rhs.entry( i ), rhs.file_name( i ) );
    return result;
}

// ********************************************************************************

