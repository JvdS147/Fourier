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
#include "SpaceGroup.h"
#include "Matrix3D.h"
#include "3DCalculations.h"
#include "Utilities.h"

#define MACRO_LIST_OF_FILES_AS_ARGUMENT \
        if ( argc == 1 ) \
            throw std::runtime_error( "Please give the names of one or more files as argument." ); \
        FileList file_list; \
        std::vector< FileName > files; \
        for ( int i( 1 ); i != argc; ++i ) \
            files.push_back( FileName( argv[ i ] ) ); \
        file_list = FileList( files );


int main( int argc, char** argv )
{
    try // Centred to primitive, list of files
    {
        MACRO_LIST_OF_FILES_AS_ARGUMENT
        for ( size_t i( 0 ); i != file_list.size(); ++i )
        {
            CrystalStructure crystal_structure;
            std::cout << "Now reading cif... " + file_list.value( i ).full_name() << std::endl;
            read_cif( file_list.value( i ), crystal_structure );

            // Extract centring vectors from space-group symmetry operators.
            std::string centring = crystal_structure.space_group().centring();
            Matrix3D centred2primitive;
            if ( centring == "A" )
                centred2primitive = A_centred_to_primitive();
            else if ( centring == "B" )
                centred2primitive = B_centred_to_primitive();
            else if ( centring == "C" )
                centred2primitive = C_centred_to_primitive();
            else if ( centring == "I" )
                centred2primitive = I_centred_to_primitive();
            else if ( centring == "F" )
                centred2primitive = F_centred_to_primitive();
            else if ( centring == "R" )
                centred2primitive = R_centred_to_primitive();
            else
                continue;
            // Apply transformation to primitive.
            crystal_structure.transform( centred2primitive );

            SpaceGroup space_group = crystal_structure.space_group();
            space_group.set_name( "" );
            space_group.remove_duplicate_symmetry_operators();
            crystal_structure.set_space_group( space_group );

            // Search for transformations that make the unit-cell angles closer to 90 degrees.
            CrystalLattice old_crystal_lattice = crystal_structure.crystal_lattice();
            Angle smallest_deviation = Angle::angle_180_degrees();
            int limit = 3;
            std::vector< Matrix3D > best_matrices;
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
                        CrystalLattice new_lattice( old_crystal_lattice );
                        // Transform.
                        Matrix3D transformation_matrix( i1, i2, i3, j1, j2, j3, k1, k2, k3 );
                        if ( ! nearly_equal( transformation_matrix.determinant(), 1.0 ) )
                            continue;
                        new_lattice.transform( transformation_matrix );
                        Angle greatest_deviation =                         absolute( new_lattice.alpha() - Angle::angle_90_degrees() );
                        greatest_deviation = std::max( greatest_deviation, absolute( new_lattice.beta()  - Angle::angle_90_degrees() ) );
                        greatest_deviation = std::max( greatest_deviation, absolute( new_lattice.gamma() - Angle::angle_90_degrees() ) );
                        // If this is a significant improvement, not just an equivalent solution that is numerically insignificantly better,
                        // wipe all we have stored and start from scratch.
                        if ( greatest_deviation < ( smallest_deviation - Angle::from_degrees( 0.1 ) ) )
                        {
                            best_matrices.clear();
                            smallest_deviation = greatest_deviation;
                        }
                        if ( nearly_equal( greatest_deviation, smallest_deviation, Angle::from_degrees( 0.001 ) ) )
                            best_matrices.push_back( transformation_matrix );
                    }
                    }
                    }
                }
                }
                }
            }
            }
            }
            // Select the best transformation. As close as possible to identity and all angles greater 90.
            Matrix3D identity_matrix;
            Matrix3D best_transformation_matrix;
            bool best_matrix_found( false );
            double smallest_sum_of_absolute_elements( 1000000.0 );
            for ( size_t j( 0 ); j != best_matrices.size(); ++j )
            {
                CrystalLattice new_lattice( old_crystal_lattice );
                // Transform.
                new_lattice.transform( best_matrices[j] );
                if ( ( ! nearly_equal( new_lattice.alpha(), Angle::angle_90_degrees(), Angle::from_degrees( 0.001 ) ) ) && ( new_lattice.alpha() < Angle::angle_90_degrees() ) )
                    continue;
                if ( ( ! nearly_equal( new_lattice.beta() , Angle::angle_90_degrees(), Angle::from_degrees( 0.001 ) ) ) && ( new_lattice.beta()  < Angle::angle_90_degrees() ) )
                    continue;
                if ( ( ! nearly_equal( new_lattice.gamma(), Angle::angle_90_degrees(), Angle::from_degrees( 0.001 ) ) ) && ( new_lattice.gamma() < Angle::angle_90_degrees() ) )
                    continue;
                // When we are here, alpha, beta and gamma are >= 90.
                best_matrix_found = true;
                double sum_of_absolute_elements = (best_matrices[j]-identity_matrix).sum_of_absolute_elements();
                if ( sum_of_absolute_elements < smallest_sum_of_absolute_elements )
                {
                    smallest_sum_of_absolute_elements = sum_of_absolute_elements;
                    best_transformation_matrix = best_matrices[j];
                }
            }
            if ( ! best_matrix_found )
            {
                // None of the transformations gave all angles >= 90.
                for ( size_t j( 0 ); j != best_matrices.size(); ++j )
                {
                    CrystalLattice new_lattice( old_crystal_lattice );
                    // Transform.
                    new_lattice.transform( best_matrices[j] );
                    double sum_of_absolute_elements = (best_matrices[j]-identity_matrix).sum_of_absolute_elements();
                    if ( sum_of_absolute_elements < smallest_sum_of_absolute_elements )
                    {
                        smallest_sum_of_absolute_elements = sum_of_absolute_elements;
                        best_transformation_matrix = best_matrices[j];
                    }
                }
            }
            // Apply best transformation.
            crystal_structure.transform( best_transformation_matrix );
            // Write out crystal structure with new symmetry operators and in P1 or P-1.

            crystal_structure.space_group().show();

            crystal_structure.save_cif( replace_extension( append_to_file_name( file_list.value( i ), "_reduced" ), "cif" ) );
            bool has_inversion_at_origin = crystal_structure.space_group().has_inversion_at_origin();
            crystal_structure.convert_to_P1();
            crystal_structure.save_cif( replace_extension( append_to_file_name( file_list.value( i ), "_reduced_P1" ), "cif" ) );
            if ( has_inversion_at_origin )
            {
                SpaceGroup space_group;
                space_group.add_inversion_at_origin();
                space_group.set_name( "P-1" );
                crystal_structure.set_space_group( space_group );
                crystal_structure.reduce_to_asymmetric_unit( 0.01 );
                crystal_structure.save_cif( replace_extension( append_to_file_name( file_list.value( i ), "_reduced_P-1" ), "cif" ) );
            }
            // Print inverse.
            std::cout << "Best transformation = " << std::endl;
            std::cout << best_transformation_matrix << std::endl;
            Matrix3D combined_transformation_matrix = best_transformation_matrix * centred2primitive;
            std::cout << "Combined transformation = " << std::endl;
            std::cout << combined_transformation_matrix << std::endl;
            std::cout << "Inverse = " << std::endl;
            std::cout << inverse( combined_transformation_matrix ) << std::endl;
        }
    MACRO_END_GAME
}
