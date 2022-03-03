#include "RunTests.h"
#include "TextFileReader.h"
#include "TextFileReader_2.h"
#include "TextFileWriter.h"
#include "VoidsFinder.h"
#include "FileName.h"
#include "FileList.h"
#include "ReadCif.h"
#include "Sort.h"
#include "CrystalStructure.h"
#include "EndGame.h"
#include "Utilities.h"

#define MACRO_ONE_FILELISTNAME_AS_ARGUMENT \
        if ( argc != 2 ) \
            throw std::runtime_error( "Please give the name of a FileList.txt file." ); \
        FileName file_list_file_name( argv[ 1 ] ); \
        FileList file_list( file_list_file_name ); \
        if ( file_list.empty() ) \
            throw std::runtime_error( std::string( "No files in file list " ) + file_list_file_name.full_name() );

int main( int argc, char** argv )
{

    try // Run tests.
    {
        run_tests();
    }
    catch ( std::exception & e )
    {
        std::cout << "An exception was thrown" << std::endl;
        std::cout << e.what() << std::endl;
    }
    try // Find voids for FileList.txt.
    {
        MACRO_ONE_FILELISTNAME_AS_ARGUMENT
        std::cout << "WARNING: the molecular volume is estimated assuming that the smallest molecular volume corresponds to Z'=1." << std::endl;
        std::cout << "WARNING: if the smallest molecular volume corresponds to Z'>1 or Z'<1 then the results will be wrong." << std::endl;
        TextFileWriter text_file_writer( FileName( file_list_file_name.directory(), "Voids", "txt" ) );
        std::vector< double > total_voids_volumes_per_symmetry_operator;
        std::vector< std::string > identifiers;
        std::vector< double > total_void_volumes;
        std::vector< double > molecular_volumes;
        std::vector< double > unit_cell_volumes;
        text_file_writer.write_line( "Identifier | total void volume | ( unit-cell volume - total void volume) / number of symmetry operators" );
        if ( file_list.empty() )
            return 0;
        size_t nfiles = file_list.size();
        total_voids_volumes_per_symmetry_operator.reserve( nfiles );
        identifiers.reserve( nfiles );
        total_void_volumes.reserve( nfiles );
        molecular_volumes.reserve( nfiles );
        unit_cell_volumes.reserve( nfiles );
        double smallest_molecular_volume( 0.0 );
        for ( size_t i( 0 ); i != nfiles; ++i )
        {
            identifiers.push_back( FileName( "", file_list.value( i ).file_name(), file_list.value( i ).extension() ).full_name() );
            CrystalStructure crystal_structure;
            std::cout << "Now reading cif... " + file_list.value( i ).full_name() << std::endl;
            read_cif( file_list.value( i ), crystal_structure );
            unit_cell_volumes.push_back( crystal_structure.crystal_lattice().volume() );
            crystal_structure.apply_space_group_symmetry();
            double total_void_volume = find_voids( crystal_structure );
            total_void_volumes.push_back( total_void_volume );
            total_voids_volumes_per_symmetry_operator.push_back( total_void_volume / crystal_structure.space_group().nsymmetry_operators() );
            double molecular_volume = ( crystal_structure.crystal_lattice().volume() - total_void_volume ) / crystal_structure.space_group().nsymmetry_operators();
            molecular_volumes.push_back( molecular_volume );
            if ( ( i == 0 ) || ( molecular_volume < smallest_molecular_volume ) )
                smallest_molecular_volume = molecular_volume;
            text_file_writer.write_line( FileName( "", file_list.value( i ).file_name(), file_list.value( i ).extension() ).full_name() + " " +
                                         double2string( total_void_volume ) + " " +
                                         double2string( molecular_volume ) );
        }
        std::vector< double > voids_volumes_per_Z;
        for ( size_t i( 0 ); i != nfiles; ++i )
        {
            // round_to_int( molecular_volumes[i] / smallest_molecular_volume ) = Z'
            voids_volumes_per_Z.push_back( total_voids_volumes_per_symmetry_operator[i] / round_to_int( molecular_volumes[i] / smallest_molecular_volume ) );
        }
        text_file_writer.write_line( "##### customer specific #####" );
        text_file_writer.write_line( "Rank/Form Void volume Void fraction" );
        text_file_writer.write_line( "              [A3/Z]              " );
        for ( size_t i( 0 ); i != nfiles; ++i )
        {
            if ( voids_volumes_per_Z[i] < 0.000001 )
                continue;
            text_file_writer.write_line( FileName( "", file_list.value( i ).file_name(), "" ).full_name() + " " +
                                         double2string_2( voids_volumes_per_Z[i], 0 ) + " " +
                                         double2string_2( 100.0 * ( total_void_volumes[i]/unit_cell_volumes[i] ), 2 ) + "%" );
        }
        Mapping sorted_map = sort( voids_volumes_per_Z );
        size_t iStart;
        for ( iStart = 0; iStart != nfiles; ++iStart )
        {
            if ( voids_volumes_per_Z[ sorted_map[iStart] ] > 20.0 )
                break;
        }

        if ( iStart == nfiles )
        {
            text_file_writer.write_line( "There are no voids greater than 20 A3/Z." );
        }
        else
        {
            text_file_writer.write_line( "##### sorted #####" );
            for ( size_t i( iStart ); i != nfiles; ++i )
                text_file_writer.write_line( identifiers[ sorted_map[i] ] + " " + double2string( voids_volumes_per_Z[ sorted_map[i] ] ) );
            text_file_writer.write_line();
            if ( (nfiles - iStart) == 1 )
            {
                text_file_writer.write( "Rank " );
                text_file_writer.write( size_t2string( sorted_map[iStart] + 1 ) );
                text_file_writer.write( " contains voids amounting to " );
                text_file_writer.write( double2string_2( voids_volumes_per_Z[sorted_map[iStart]], 0 ) );
                text_file_writer.write( " A3/Z." );
            }
            else
            {
//        Ranks 12, 22 5, 17, 1, 9 and 10 contain voids amounting to 20, 21, 21, 24, 28, 40 and 45 �3/Z, respectively.
                text_file_writer.write( "Ranks " );
                for ( size_t i( iStart ); i != nfiles; ++i )
                {
                    if ( i == nfiles - 1 )
                        text_file_writer.write( " and "  );
                    else if ( i != iStart )
                        text_file_writer.write( ", "  );
                    text_file_writer.write( size_t2string( sorted_map[i] + 1 ) );
                }
                text_file_writer.write( " contain voids amounting to " );
                for ( size_t i( iStart ); i != nfiles; ++i )
                {
                    if ( i == nfiles - 1 )
                        text_file_writer.write( " and "  );
                    else if ( i != iStart )
                        text_file_writer.write( ", "  );
                    text_file_writer.write( double2string_2( voids_volumes_per_Z[ sorted_map[i] ], 0 ) );
                }
                text_file_writer.write( " A3/Z, respectively." );
            }
            text_file_writer.write( " Of interest are voids that are greater than about 20 A3/Z: 21.5 A3/Z suffices to store a water molecule (at least in terms of volume), a chloride ion is about 25 A3/Z." );
            text_file_writer.write_line( " Voids between 15 and 20 A3/Z are quite common, but voids over 25 A3/Z are rare." );
        }
    MACRO_END_GAME
}
