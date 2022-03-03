#include "RunTests.h"
#include "TextFileReader.h"
#include "TextFileReader_2.h"
#include "TextFileWriter.h"
#include "FileName.h"
#include "FileList.h"
#include "ReadCif.h"
#include "ReadCifOrCell.h"
#include "CrystalStructure.h"
#include "EndGame.h"
#include "Angle.h"
#include "CrystalLattice.h"
#include "Matrix3D.h"
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
    try // Find unit-cell transformation.
    {
        if ( argc != 3 )
            throw std::runtime_error( "Please give the name of two .cif files." );
        FileName input_file_name( argv[ 1 ] );
        CrystalStructure crystal_structure;
        read_cif_or_cell( input_file_name, crystal_structure );
        if ( false ) // C-centred to primitive
        {
            crystal_structure.transform( Matrix3D(  1.0,  0.0,  0.0,
                                                    0.0,  0.5, -0.5,
                                                    0.0,  0.5,  0.5 ) );
        }
        CrystalLattice old_crystal_lattice = crystal_structure.crystal_lattice();
        FileName file_name_2( argv[ 2 ] );
        CrystalStructure crystal_structure_2;
        read_cif_or_cell( file_name_2, crystal_structure_2 );
        if ( false ) // C-centred to primitive
        {
            crystal_structure_2.transform( Matrix3D(  0.5,  0.5,  0.0,
                                                     -0.5,  0.5,  0.0,
                                                      0.0,  0.0,  1.0 ) );
        }
        CrystalLattice target_crystal_lattice = crystal_structure_2.crystal_lattice();
        double length_tolerance_percent( 10.0 );
        Angle angle_tolerance = Angle::from_degrees( 10.0 );
        int limit = 5;
        for ( int i1( -limit ); i1 != limit+1; ++i1 )
        {
            std::cout << " ." << std::endl;
        for ( int i2( -limit ); i2 != limit+1; ++i2 )
        {
            std::cout << " .." << std::endl;
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
                    // Make a copy
                    CrystalLattice new_lattice( old_crystal_lattice );
                    // Transform
                    Matrix3D transformation_matrix( i1, i2, i3, j1, j2, j3, k1, k2, k3 );
                    if ( ! nearly_equal( transformation_matrix.determinant(), 1.0 ) )
                        continue;
                    new_lattice.transform( transformation_matrix );
                    if ( ! nearly_equal( new_lattice, target_crystal_lattice, length_tolerance_percent, angle_tolerance ) )
                        continue;
                    {
                        transformation_matrix.show();
                        new_lattice.print();
                        std::cout << ((target_crystal_lattice.a_vector()+target_crystal_lattice.b_vector()+target_crystal_lattice.c_vector()) - (new_lattice.a_vector()+new_lattice.b_vector()+new_lattice.c_vector())).length() << std::endl;
                        std::cout << std::endl;
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
    MACRO_END_GAME
}


