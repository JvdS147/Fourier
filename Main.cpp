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

// Comment regarding the licence:
// These files are deliberately licensed under the most lenient licence I could find.
// The following is allowed:
// Copy the source code for free, change the source code, compile the source code, sell the compiled code for money. As long as C. J. van de Streek is acknowledged and is not used to endorse the new code.

#include "3DCalculations.h"
#include "AddClass.h"
#include "AnalyseRings.h"
#include "AnalyseTrajectory.h"
#include "Angle.h"
#include "AnisotropicDisplacementParameters.h"
//#include "BFDH.h"
#include "BondDetector.h"
#include "ChebyshevBackground.h"
#include "CheckFoundItem.h"
#include "ChemicalFormula.h"
#include "CollectionOfPoints.h"
#include "Constraints.h"
#include "CorrelationMatrix.h"
#include "CrystalStructure.h"
#include "CyclicInteger.h"
#include "DoubleWithESD.h"
#include "DrunkardsWalk.h"
#include "Eigenvalue.h"
#include "EndGame.h"
#include "FileList.h"
#include "FileName.h"
#include "Finish_inp.h"
#include "Fraction.h"
#include "GeneratePowderCIF.h"
#include "Histogram.h"
#include "InpWriter.h"
#include "LabelsAndShieldings.h"
#include "MathsFunctions.h"
#include "ModelBuilding.h"
#include "Plane.h"
#include "PowderMatchTable.h"
#include "PowderPattern.h"
#include "PowderPatternCalculator.h"
#include "RandomNumberGenerator.h"
#include "ReadCell.h"
#include "ReadCif.h"
#include "ReadCifOrCell.h"
#include "ReadXSD.h"
#include "ReadXYZ.h"
#include "Refcode.h"
#include "ReflectionList.h"
#include "RunningAverageAndESD.h"
#include "RunTests.h"
#include "SimilarityAnalysis.h"
#include "SkipBo.h"
#include "Sort.h"
#include "Sudoku.h"
#include "SudokuSolver.h"
#include "SymmetryOperator.h"
#include "TextFileReader.h"
#include "TextFileReader_2.h"
#include "TextFileWriter.h"
#include "TLS.h"
#include "TLSWriter.h"
#include "TOPAS.h"
#include "Utilities.h"
#include "VoidsFinder.h"
#include "WriteCASTEPFile.h"

#include <vector>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <fstream>

//template<typename T>
//std::vector< T > convert_file_list_to_vector( const FileList & file_list )
//{
//    std::vector< T > result;
//    result.reserve( file_list.size() );
//    for ( size_t i( 0 ); i != file_list.size(); ++i )
//    {
//        result.push_back( T( file_list.value( i ) ) );
//    }
//    return result;
//}

void compare( const std::string & name_structure_1, const std::string & name_structure_2, TextFileWriter & text_file_writer )
{
    CrystalStructure crystal_structure_1;
    read_cif( FileName( "C:\\Data\\ActaCryst_powder\\LSERIN\\", name_structure_1, "cif" ), crystal_structure_1 );
    CrystalStructure crystal_structure_2;
    read_cif( FileName( "C:\\Data\\ActaCryst_powder\\LSERIN\\", name_structure_2, "cif" ), crystal_structure_2 );
    double result = root_mean_square_Cartesian_displacement( crystal_structure_1, crystal_structure_2, false );
    text_file_writer.write_line( name_structure_1 + " " + name_structure_2 + " " + double2string( result ) );
}

struct C6Record
{
    size_t element_1;
    size_t element_2;
    double CN_1;
    double CN_2;
    double C6;
};

struct SimulatedPowderPatternCrystalStructure
{
    SimulatedPowderPatternCrystalStructure( CrystalStructure & crystal_structure ):
    crystal_structure_(crystal_structure),
    total_signal_normalisation_(10000.0),
    include_PO_(false),
    PO_direction_(MillerIndices(0,0,1)),
    PO_extent_(1.0),
    FWHM_(0.1),
    include_background_(true),
    background_total_signal_normalisation_(10000.0)
    {}


    // The PO does NOT apply to the amorphous contribution
    CrystalStructure crystal_structure_;
    double total_signal_normalisation_;
    bool   include_PO_;
    MillerIndices PO_direction_;
    double PO_extent_;
    double FWHM_;
    bool   include_background_;
    double background_total_signal_normalisation_;
};

#define MACRO_ONE_FILELISTNAME_AS_ARGUMENT \
        if ( argc != 2 ) \
            throw std::runtime_error( "Please give the name of a FileList.txt file." ); \
        FileName file_list_file_name( argv[ 1 ] ); \
        FileList file_list( file_list_file_name ); \
        if ( file_list.empty() ) \
            throw std::runtime_error( std::string( "No files in file list " ) + file_list_file_name.full_name() );

#define MACRO_ONE_CIFFILENAME_AS_ARGUMENT \
        if ( argc != 2 ) \
            throw std::runtime_error( "Please give the name of a .cif or .cell file." ); \
        FileName input_file_name( argv[ 1 ] ); \
        CrystalStructure crystal_structure; \
        read_cif_or_cell( input_file_name, crystal_structure );

#define MACRO_ONE_XYEFILENAME_AS_ARGUMENT \
        if ( argc != 2 ) \
            throw std::runtime_error( "Please give the name of a .xye file." ); \
        FileName input_file_name( argv[ 1 ] ); \
        PowderPattern powder_pattern; \
        powder_pattern.read_xye( input_file_name );

#ifndef BACKGROUND_TOTAL_SIGNAL_NORMALISATION
#define BACKGROUND_TOTAL_SIGNAL_NORMALISATION 5000.0
#endif

#ifndef NORMALISE_HIGHEST_PEAK
#define NORMALISE_HIGHEST_PEAK 10000.0
#endif

#ifndef ZERO_POINT_ERROR
#define ZERO_POINT_ERROR 0.02
#endif

#ifndef PREFERRED_ORIENTATION
#define PREFERRED_ORIENTATION 0.9
#endif

#ifndef FULL_WIDTH_HALF_MAX
#define FULL_WIDTH_HALF_MAX 0.25
#endif

// Test if .h files compile stand alone
// FractionalCoordinate / OrthonormalCoordinate
// Smart pointers. Because I currently have no smart pointers, it is almost impossible to hold objects by pointer so everything is always copied.
// Add class "CheckedItemReadFromFile" or "CheckedItemAssignedValue" Store reference (&) to an existing variable
// + bool, upon destruction writes error message if variable has not been assigned a value.
// Maybe a "format checker and beautifier", e.g. an atom label must be "Aa11", a refcode must be "AAAAAA11"
// A "Distance" and "Length" class, which for < and = uses the norm2 and calculates the norm on demand (so cached).
// void check_if_quotes_correct( const std::string & input );
// Need a "add_two_strings_with_quotes()"
// Alternatively, a QuotedString class.
// Add Marcus' cell deformation measure.
// Number of molecules for hemihydrate must be 1 API plus 1/2 H2O
// Could now add CrystalStructure::CrystalStructure( const FileName & )
// FileList::save()
// A "DoubleChecked" class: check that adding two numbers changes the result (i.e. check if one is too small to be added to the other),
//    a bool to keep track whether the value has been initialised. Function to convert everything smaller than e.g. 1E-6 to 0.0.
//    add specialised function for adding a std::vector of numbers. Is there anything special that can be done by alternating positive and negative numbers? Dimension (as in the unit, kg, meter, second) as std::string
//    internal bool for "very big" (infinity), although "not initialised" may fulfil that purpose
// If we want British English, Math should be Maths.
// Topological equivalence of torsion angles.
// Two kinds of crystal structure: asymmetric unit only and all molecules in one unit cell.
// Add appending + checks for existence to TextFileWriter
// Remove trailing whitespace from text file
// The "finish .inp after DASH 6" could reshuffle some keywords (move 2theta step to just after 2theta begin and 2theta end, for example) and beautify some others.
// Replace Rietveld by Loopstra-Rietveld
// Maybe a file, or class, CrystallographicCalculator. The space group and the lattice are always needed when e.g. calculating if two reflections are equivalent (currently member functions of PowderPatternCalculator).
// Clean up Hill System in ChemicalFormula.
// Could implement a "VeryBigInteger" class by storing a std::vector< size_t > which stores e.g. the number 8242 as [ 8, 2, 4, 2 ]. MAXINT = 2,147,483,647
// std::string centring_vectors_to_string( const std::vector< Vector3D > & centring_vectors )
// After transformation, check that the structure has not changed by calculating the XRPD patterns
// Simulate a single crystal diffraction pattern
// We currently have several different ways to check if two unit cells are the same. These could all be moved to the CrystalLattice class. Some tolerances are absolute, some are a percentage.
//    For lengths, a percentage makes sense.

void add_all_torsions_in_ring( const std::vector< std::string > & ring_labels, std::vector< std::vector< std::string > > & torsions )
{
    for ( size_t i( 0 ); i != ring_labels.size(); ++i )
    {
        std::vector< std::string > torsion;
        CyclicInteger ci( 0, ring_labels.size() - 1, i );
        torsion.push_back( ring_labels[ci.next_value()] );
        torsion.push_back( ring_labels[ci.next_value()] );
        torsion.push_back( ring_labels[ci.next_value()] );
        torsion.push_back( ring_labels[ci.next_value()] );
        torsions.push_back( torsion );
    }
}

struct DisorderGroup
{
    // Poor data structure design, should use std::pair to force the two to have equal size
    // and to make it easier to identify which atoms corresponds to which.
    std::vector< std::string > major_occupancy_labels_;
    std::vector< std::string > minor_occupancy_labels_;
};

// Helper function for converting .gzmat to .zmatrix.
std::string extract_variable_value( const std::string & line )
{
    Splitter splitter( "=" );
    splitter.set_merge_delimiters( false );
    std::vector< std::string > words = splitter.split( line );
    if ( words.size() != 2 )
        throw std::runtime_error( "extract_variable_value(): ERROR: no variable / value pair found." );
    return strip( words[1] );
}

std::string clean_string( const std::string & input, const size_t power )
{
    std::string result( input );
    bool negative = ( input[0] == '-' );
    if ( negative )
        result[0] = '0';
    result = "00000000" + result;
    size_t iPos = result.find( "(" );
    if ( iPos == std::string::npos )
        iPos = result.length();
    result.insert( iPos-power, "." );
    result = remove_from_start( result, '0' );
    if ( result[0] == '.' )
        result = "0" + result;
    if ( negative )
        result = "-" + result;
    return pad_plus( result, 12 );
}

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

    try // Average two unit cells.
    {
        if ( argc != 3 )
            throw std::runtime_error( "Please give the name of two .cif files." );
        FileName file_name_1( argv[ 1 ] );
        CrystalStructure crystal_structure_1;
        read_cif_or_cell( file_name_1, crystal_structure_1 );
        FileName file_name_2( argv[ 2 ] );
        CrystalStructure crystal_structure_2;
        read_cif_or_cell( file_name_2, crystal_structure_2 );
        CrystalLattice average_crystal_lattice = average( crystal_structure_1.crystal_lattice(), crystal_structure_2.crystal_lattice() );
        crystal_structure_1.set_crystal_lattice( average_crystal_lattice );
        crystal_structure_2.set_crystal_lattice( average_crystal_lattice );
        crystal_structure_1.save_cif( replace_extension( append_to_file_name( file_name_1, "_avguc" ) , "cif" ) );
        crystal_structure_2.save_cif( replace_extension( append_to_file_name( file_name_2, "_avguc" ) , "cif" ) );
    MACRO_END_GAME

    // Transformation of the crystal structure (unit cell + atomic coordinates including ADPs + space group)
    // followed by a transformation of the atomic coordinates including ADPs.
    try
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        Vector3D com = crystal_structure.centre_of_mass();
        std::cout << "Centre of mass = " << std::endl;
        com.show();
        Matrix3D tranformation_matrix(  4.0,  0.0,  1.0,
                                        0.0, -1.0,  0.0,
                                        1.0,  0.0,  0.0 );
        crystal_structure.transform( tranformation_matrix );
        std::cout << "Inverse transformation matrix:" << std::endl;
        std::cout << inverse( tranformation_matrix ) << std::endl;
        // R-centred to primitive
//        crystal_structure.transform( Matrix3D(  2.0/3.0,  1.0/3.0,  1.0/3.0,
//                                                1.0/3.0,  2.0/3.0,  2.0/3.0,
//                                                0.0,  0.0,  1.0 ) );
        if ( (false) )
        {
            SymmetryOperator symmetry_operator( "x-1/4,y,z-3/4" );
            SpaceGroup space_group = crystal_structure.space_group();
            space_group.set_name( "" );
            space_group.apply_similarity_transformation( symmetry_operator );
         //   space_group.add_inversion_at_origin();
            crystal_structure.set_space_group( space_group );
        }
        Matrix3D rotation(  1.0,  0.0,  0.0,
                            0.0,  1.0,  0.0,
                            0.0,  0.0,  1.0 );
        Vector3D shift( 0.0, 0.0, 0.0 );
//        shift -= com;
//        shift += com_2;
        for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
        {
            Atom new_atom( crystal_structure.atom( i ) );
//            new_atom.set_position( new_atom.position() - shift );
//            new_atom.set_position( rotation * new_atom.position() );
//            new_atom.set_position( new_atom.position() + shift );
            new_atom.set_position( ( rotation * crystal_structure.atom( i ).position() ) + shift );
            if ( new_atom.ADPs_type() == Atom::ANISOTROPIC )
                new_atom.set_anisotropic_displacement_parameters( rotate_adps( new_atom.anisotropic_displacement_parameters(), rotation, crystal_structure.crystal_lattice() ) );
            crystal_structure.set_atom( i, new_atom );
        }
        // In Mercury, if the space-group name and the set of symmetry operators do not match up,
        // the space-group name takes precedence, so we have to erase it to ensure that the
        // symmetry operators are used instead.
        if ( (true) )
        {
            SpaceGroup space_group = crystal_structure.space_group();
            space_group.set_name( "" );
            space_group.remove_duplicate_symmetry_operators();
            crystal_structure.set_space_group( space_group );
        }
        crystal_structure.save_cif( replace_extension( append_to_file_name( input_file_name, "_transformed" ) , "cif" ) );
    MACRO_END_GAME

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

    try // Find unit-cell angles close to 90 degrees.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        CrystalLattice old_crystal_lattice = crystal_structure.crystal_lattice();
        // If the unit-cell metric was monoclinic with alpha = gamma = 90, then we want to keep it that way.
        bool is_monoclinic( true );
        Angle monoclinic_tolerance = Angle::from_degrees( 0.001 );
        if ( ! nearly_equal( old_crystal_lattice.alpha(), Angle::from_degrees( 90.0 ), monoclinic_tolerance ) )
            is_monoclinic = false;
        if ( ! nearly_equal( old_crystal_lattice.gamma(), Angle::from_degrees( 90.0 ), monoclinic_tolerance ) )
            is_monoclinic = false;
        Matrix3D identity_matrix;
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
                    // Make a copy
                    CrystalLattice new_lattice( old_crystal_lattice );
                    // Transform
                    Matrix3D transformation_matrix( i1, i2, i3, j1, j2, j3, k1, k2, k3 );
                    if ( ! nearly_equal( transformation_matrix.determinant(), 1.0 ) )
                        continue;
                    new_lattice.transform( transformation_matrix );
                    // If the unit-cell metric was monoclinic with alpha = gamma = 90, then we want to keep it that way.
                    if ( is_monoclinic )
                    {
                        if ( ! nearly_equal( new_lattice.alpha(), Angle::from_degrees( 90.0 ), monoclinic_tolerance ) )
                            continue;
                        if ( ! nearly_equal( new_lattice.gamma(), Angle::from_degrees( 90.0 ), monoclinic_tolerance ) )
                            continue;
                    }
                    Angle limit = Angle::from_degrees( 80.0 );
                    if ( ( new_lattice.alpha() < limit ) || ( new_lattice.alpha() > ( Angle::from_degrees( 180.0 ) - limit ) ) )
                        continue;
                    if ( ( new_lattice.beta()  < limit ) || ( new_lattice.beta()  > ( Angle::from_degrees( 180.0 ) - limit ) ) )
                        continue;
                    if ( ( new_lattice.gamma() < limit ) || ( new_lattice.gamma() > ( Angle::from_degrees( 180.0 ) - limit ) ) )
                        continue;
                    {
                        transformation_matrix.show();
                        new_lattice.print();
                        std::cout << "  " << (transformation_matrix-identity_matrix).sum_of_absolute_elements() << std::endl;
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

    try // RMSCD with / without matching.
    {
        if ( argc != 3 )
            throw std::runtime_error( "Please give the names of two .cif or .cell files." );
        FileName file_name_1( argv[ 1 ] );

  //      FileName file_name_1( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Computational\\fixed_cell\\0028\\ADERIL\\ADERIL.cell" );
        CrystalStructure crystal_structure_1;
        read_cif_or_cell( file_name_1, crystal_structure_1 );
        if ( to_upper( file_name_1.extension() ) == "CELL" )
            crystal_structure_1.save_cif( replace_extension( file_name_1, "cif" ) );
        else
            crystal_structure_1.reduce_to_asymmetric_unit();
        FileName file_name_2( argv[ 2 ] );
//        FileName file_name_2( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Computational\\fixed_cell\\0028\\ADERIL01\\ADERIL01.cell" );
        CrystalStructure crystal_structure_2;
        read_cif_or_cell( file_name_2, crystal_structure_2 );
        if ( to_upper( file_name_2.extension() ) == "CELL" )
            crystal_structure_2.save_cif( replace_extension( file_name_2, "cif" ) );
        else
            crystal_structure_1.reduce_to_asymmetric_unit();
//        crystal_structure_1.set_space_group( SpaceGroup() );
//        crystal_structure_2.set_space_group( SpaceGroup() );
        double result = RMSCD_with_matching( crystal_structure_1, crystal_structure_2, 2, false );
//        double result = root_mean_square_Cartesian_displacement( crystal_structure_1, crystal_structure_2, false );
        std::cout << result << std::endl;
    MACRO_END_GAME

    try // RMSCD of a structure and its energy-minimised counterpart. Without matching.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        std::string file_name_str = input_file_name.file_name();
        FileName file_name_2;
        if ( ( file_name_str.length() > 8 ) && ( ( file_name_str.substr( file_name_str.length() - 8 ) == "_mi_ucfr" ) || ( file_name_str.substr( file_name_str.length() - 8 ) == "_mi_ucfx" ) ) )
            file_name_2 = FileName( input_file_name.directory(), file_name_str.substr( 0, file_name_str.length() - 8 ), input_file_name.extension() );
        else if ( append_to_file_name( input_file_name, "_mi_ucfr" ).exists() )
            file_name_2 = append_to_file_name( input_file_name, "_mi_ucfr" );
        else if ( append_to_file_name( input_file_name, "_mi_ucfx" ).exists() )
            file_name_2 = append_to_file_name( input_file_name, "_mi_ucfx" );
        else
            throw std::runtime_error( "Please give the name of a .cif file such that it can be matched with its energy minimised counterpart." );
        CrystalStructure crystal_structure_2;
        read_cif( file_name_2, crystal_structure_2 );
        double result = root_mean_square_Cartesian_displacement( crystal_structure, crystal_structure_2, false );
        std::cout << "Structure 1 = " << input_file_name.full_name() << std::endl;
        std::cout << "Structure 2 = " << file_name_2.full_name() << std::endl;
        std::cout << result << std::endl;
    MACRO_END_GAME

    // Simulate an experimental powder diffraction pattern
    try
    {

     //   MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        MACRO_ONE_FILELISTNAME_AS_ARGUMENT
     //   std::vector< SimulatedPowderPatternCrystalStructure > crystal_structures;
//    total_signal_normalisation_(10000.0),
//    include_PO_(false),
//    PO_direction_(MillerIndices(0,0,1)),
//    PO_extent_(1.0),
//    FWHM_(0.1),
//    include_background_(true),
//    background_total_signal_normalisation_(10000.0)
      //  SimulatedPowderPatternCrystalStructure sim_XRPD_crystal_structure( crystal_structure );
        double wavelength( 1.54056 );
        // ######################## CHANGE THIS ##################################
        double zero_point_error = ZERO_POINT_ERROR; // 0.06;
        Angle two_theta_start( 1.0, Angle::DEGREES );
        Angle two_theta_end(  35.0, Angle::DEGREES );
        Angle two_theta_step( 0.015, Angle::DEGREES );
        //std::cout << NORMALISE_HIGHEST_PEAK << " " << BACKGROUND_TOTAL_SIGNAL_NORMALISATION << std::endl;
        for ( size_t i( 0 ); i != file_list.size(); ++i )
        {
            CrystalStructure crystal_structure;
            std::cout << "Now reading cif... " + file_list.value( i ).full_name() << std::endl;
            read_cif( file_list.value( i ), crystal_structure );
            SimulatedPowderPatternCrystalStructure sim_XRPD_crystal_structure( crystal_structure );
            // ######################## CHANGE THIS ##################################
            sim_XRPD_crystal_structure.background_total_signal_normalisation_ = BACKGROUND_TOTAL_SIGNAL_NORMALISATION;
            sim_XRPD_crystal_structure.include_PO_ = true;
          //  sim_XRPD_crystal_structure.include_PO_ = false;
            if ( sim_XRPD_crystal_structure.include_PO_ )
            {
                // ######################## CHANGE THIS ##################################
                sim_XRPD_crystal_structure.PO_extent_= PREFERRED_ORIENTATION;
                switch ( crystal_structure.crystal_lattice().lattice_system() )
                {
                    case CrystalLattice::TRICLINIC    :
                    case CrystalLattice::ORTHORHOMBIC : {
                                                            CrystalLattice crystal_lattice = crystal_structure.crystal_lattice();
                                                            if ( crystal_lattice.a() < crystal_lattice.b() )
                                                            {
                                                                if ( crystal_lattice.a() < crystal_lattice.c() )
                                                                    sim_XRPD_crystal_structure.PO_direction_ = MillerIndices( 1, 0, 0 );
                                                            }
                                                            else // b < a
                                                            {
                                                                if ( crystal_lattice.b() < crystal_lattice.c() )
                                                                    sim_XRPD_crystal_structure.PO_direction_ = MillerIndices( 0, 1, 0 );
                                                            }
                                                            break;
                                                        }
                    case CrystalLattice::MONOCLINIC   : {
                                                            sim_XRPD_crystal_structure.PO_direction_ = MillerIndices( 0, 1, 0 );
                                                            break;
                                                        }
                    case CrystalLattice::TRIGONAL     :
                    case CrystalLattice::TETRAGONAL   :
                    case CrystalLattice::HEXAGONAL    : break; // Nothing to do
                    case CrystalLattice::RHOMBOHEDRAL :
                    case CrystalLattice::CUBIC        : sim_XRPD_crystal_structure.include_PO_ = false;
                }
            }
            // ######################## CHANGE THIS ##################################
            sim_XRPD_crystal_structure.FWHM_ = FULL_WIDTH_HALF_MAX; // 0.2; // ?
      //      crystal_structures.push_back( sim_XRPD_crystal_structure );
            PowderPattern result( two_theta_start, two_theta_end, two_theta_step );
    //        for ( size_t i( 0 ); i != crystal_structures.size(); ++i )
            {
                sim_XRPD_crystal_structure.crystal_structure_.apply_space_group_symmetry();
                std::cout << "Now calculating powder pattern... " << std::endl;
                PowderPatternCalculator powder_pattern_calculator( sim_XRPD_crystal_structure.crystal_structure_ );
                powder_pattern_calculator.set_wavelength( wavelength );
                powder_pattern_calculator.set_two_theta_start( two_theta_start );
                powder_pattern_calculator.set_two_theta_end( two_theta_end );
                powder_pattern_calculator.set_two_theta_step( two_theta_step );
                powder_pattern_calculator.set_FWHM( sim_XRPD_crystal_structure.FWHM_ );
                if ( sim_XRPD_crystal_structure.include_PO_ )
                    powder_pattern_calculator.set_preferred_orientation( sim_XRPD_crystal_structure.PO_direction_, sim_XRPD_crystal_structure.PO_extent_ );
                PowderPattern powder_pattern;
                powder_pattern_calculator.calculate( powder_pattern );
                powder_pattern.normalise_total_signal( sim_XRPD_crystal_structure.total_signal_normalisation_ );
                result += powder_pattern;
                if ( sim_XRPD_crystal_structure.include_background_ )
                {
                    PowderPatternCalculator background_powder_pattern_calculator( sim_XRPD_crystal_structure.crystal_structure_ );
                    background_powder_pattern_calculator.set_wavelength( wavelength );
                    background_powder_pattern_calculator.set_two_theta_start( two_theta_start );
                    background_powder_pattern_calculator.set_two_theta_end( two_theta_end );
                    background_powder_pattern_calculator.set_two_theta_step( two_theta_step );
                    background_powder_pattern_calculator.set_FWHM( 5.0 );
                    // We never include PO for the amorphous background
                    PowderPattern powder_pattern;
                    background_powder_pattern_calculator.calculate( powder_pattern );
                    powder_pattern.normalise_total_signal( sim_XRPD_crystal_structure.background_total_signal_normalisation_ );
                    result += powder_pattern;
                }
            }

          //  PowderPattern background_1;
          //  background_1.read_xye( FileName( "W:\\GC\\Form_I_BKGR.xye" ) );
          //  background_1.normalise( 1000.0 );
          //  result += background_1;

            // ######################## CHANGE THIS ##################################
            result.normalise_highest_peak( NORMALISE_HIGHEST_PEAK ); // 300.0
            result.add_constant_background( 20.0 );
            // Next line is partially superfluous: PowderPattern::add_Poisson_noise() converts counts to integers and the result of the Poisson distribution is a size_t.
            // So the final powder pattern does not change, but it changes the values stored in PowderPattern::noise_
            result.make_counts_integer();
            result.add_Poisson_noise();
            // They are *estimated* standard deviations, so they should be calculated *after* the noise has been introduced.
            result.recalculate_estimated_standard_deviations();
            result.correct_zero_point_error( Angle::from_degrees( -zero_point_error ) );

            PowderPattern background = calculate_Brueckner_background( result,
                                                                       50, // niterations
                                                                       round_to_int( 50.0 * ( Angle::from_degrees( 0.015 ) / result.average_two_theta_step() ) ), // window
                                                                       true, // apply_smoothing
                                                                       5 ); // smoothing_window
       //     result -= background;

            result.save_xye( replace_extension( file_list.value( i ), "xye" ), true );
        }
    MACRO_END_GAME

    try // Print chemical formula for molecule from .cif or .xyz file.
    {
        if ( argc == 1 )
        {
            char a;
            std::cout << "Usage:" << std::endl;
            std::cout << std::endl;
            std::cout << "ChemicalFormula.exe <CrystalStructure.cif>" << std::endl;
            std::cout << "ChemicalFormula.exe <Molecule.xyz>" << std::endl;
            std::cout << std::endl;
            std::cout << "Output: a chemical formula." << std::endl;
            std::cin >> a;
            return 0;
        }
        FileName file_name( argv[ 1 ] );
        std::string extension = to_lower( file_name.extension() );
        ChemicalFormula chemical_formula;
        size_t tot_atoms( 0 );
        if ( extension == "cif" )
        {
            std::cout << "Now reading cif... " + file_name.full_name() << std::endl;
            CrystalStructure crystal_structure;
            read_cif( file_name, crystal_structure );
            tot_atoms = crystal_structure.natoms();
            for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
                chemical_formula.add_element( crystal_structure.atom( i ).element() );
        }
        else if ( extension == "xyz" )
        {
            std::cout << "Now reading xyz... " + file_name.full_name() << std::endl;
            std::vector< Atom > atoms;
            read_xyz( file_name, atoms );
            tot_atoms = atoms.size();
            for ( size_t i( 0 ); i != atoms.size(); ++i )
                chemical_formula.add_element( atoms[ i ].element() );
		}
        std::cout << chemical_formula.to_string( true, true ) << std::endl;
        std::cout << chemical_formula.to_string( false, false ) << std::endl;
        std::cout << "Total number of atoms = " << tot_atoms << std::endl;
        std::cout << double2string( chemical_formula.solid_state_volume() ) << " A3" << std::endl;
        std::cout << double2string( chemical_formula.molecular_weight() ) << " g/mol" << std::endl;
    MACRO_END_GAME

    try // TRIZIN04. The ADPs never worked...
    {

        double a = 6.884;
        double b = 9.569;
        double c = 7.093;
        Angle beta = Angle::from_degrees( 126.61 );

        CrystalLattice crystal_lattice_reduced( a,
                                                b,
                                                c,
                                                Angle::angle_90_degrees(),
                                                beta,
                                                Angle::angle_90_degrees() );
   //     std::cout << "crystal_lattice_reduced.volume()" << crystal_lattice_reduced.volume() << std::endl;
        CrystalLattice crystal_lattice_non_reduced( sqrt( 9.0*square( a ) + 4.0*square( c ) + 12.0*a*c*beta.cosine() ),
                                                    b,
                                                    c,
                                                    Angle::angle_90_degrees(),
                                                    arccotangent( beta.cotangent() + 2.0*c/( 3.0*a*beta.sine() ) ),
                                                    Angle::angle_90_degrees() );
   //     std::cout << "crystal_lattice_non_reduced.volume()" << crystal_lattice_non_reduced.volume() << std::endl;

  //      crystal_lattice_non_reduced.print();

        std::vector< SymmetryOperator > symmetry_operators;
        symmetry_operators.push_back( SymmetryOperator( "x,y,z" ) );
        symmetry_operators.push_back( SymmetryOperator( "1/2+x,1/2+y,z" ) );
        symmetry_operators.push_back( SymmetryOperator( "1/2-x,1/2+y,1/2-z" ) );
        symmetry_operators.push_back( SymmetryOperator( "-x,y,1/2-z" ) );
        symmetry_operators.push_back( SymmetryOperator( "-x,-y,-z" ) );
        symmetry_operators.push_back( SymmetryOperator( "1/2-x,1/2-y,-z" ) );
        symmetry_operators.push_back( SymmetryOperator( "1/2+x,1/2-y,1/2+z" ) );
        symmetry_operators.push_back( SymmetryOperator( "x,-y,1/2+z" ) );
        SpaceGroup space_group_C2c( symmetry_operators, "C2/c");

        Angle phi( 11.1, Angle::DEGREES );
        Angle NCN( 125.2, Angle::DEGREES );
        Angle CNC = Angle::from_degrees( 240.0 ) - NCN;
        Angle hNCN = NCN / 2.0;
        Angle hCNC = CNC / 2.0;
        double CN( 1.338 );
        double CH( 0.96 );
    //    double x = 0.0;
        double y = -0.004; // This is in Angstrom (so Cartesian)
    //    double z = c/4.0;
        Vector3D O( 0.0, y, c/4.0 );
        Angle sixty( 60.0, Angle::DEGREES );
        double tan_60 = sixty.tangent();

        double C1_x = hCNC.sine() * CN * phi.cosine();
        double C1_y = ( hCNC.sine() / tan_60 ) * CN;
        double C1_z = phi.sine() * hCNC.sine() * CN;
        Vector3D C1( C1_x, C1_y, C1_z );
        double N1_x = hNCN.sine() * CN * phi.cosine();
        double N1_y = -( hNCN.sine() / tan_60 ) * CN;
        double N1_z = phi.sine() * hNCN.sine() * CN;
        Vector3D N1( N1_x, N1_y, N1_z );
        double C2_x = 0.0;
        double C2_y = -( hNCN.sine() / tan_60 ) * CN - hNCN.cosine() * CN;
        double C2_z = 0.0;
        Vector3D C2( C2_x, C2_y, C2_z );
        double N2_x = 0.0;
        double N2_y = ( hCNC.sine() / tan_60 ) * CN + hCNC.cosine() * CN;
        double N2_z = 0.0;
        Vector3D N2( N2_x, N2_y, N2_z );
        Vector3D H1 = C1 + C1 * ( CH / C1.length() );
        double H2_x = 0.0;
        double H2_y = C2_y - CH;
        double H2_z = 0.0;
        Vector3D H2( H2_x, H2_y, H2_z );
        C1 += O;
        N1 += O;
        C2 += O;
        N2 += O;
        H1 += O;
        H2 += O;

        double from_nm = 10.0;
        double from_rad = 1.0; // radians2degrees;

        double T11 = from_nm * from_nm * 0.00024;
        double T12 = from_nm * from_nm * 0.0;
        double T13 = from_nm * from_nm * -0.00026;
        double T23 = T13; // Why ?
        double T22 = from_nm * from_nm * 0.00035;
        double T33 = from_nm * from_nm * -0.0005;

        SymmetricMatrix3D T( T11, T22, T33, T12, T13, T23 );

        double L11 = from_rad * 0.05;
        double L12 = from_rad * 0.0;
        double L13 = from_rad * -0.002;
        double L23 = L13; // Why ?
        double L22 = from_rad * 0.06;
        double L33 = from_rad * 0.004;

        SymmetricMatrix3D L( L11, L22, L33, L12, L13, L23 );

    // The following is not possible: S33 = -S11 - S22, so if S11 is 0.0, S33 and S22 must be of opposite sign.
    // I think that they chose as constraint either that S11 is 0.0 or that S22 = S33
    // Because S11 does not have an ESD and because the ESDs of S22 and S33 are different, apparently they chose S11 = 0.0.

        double S11 = from_nm * from_rad * 0.0;
        double S12 = from_nm * from_rad * 0.0;
        double S13 = from_nm * from_rad * 0.0;
        double S21 = S12;
        double S22 = from_nm * from_rad * 0.0002;
        double S23 = from_nm * from_rad * 0.0;
        double S31 = from_nm * from_rad * 0.0;
        double S32 = S23;
        double S33 = from_nm * from_rad * 0.0002;

        Matrix3D S( S11, S12, S13,
                    S21, S22, S23,
                    S31, S32, S33 );

        TLS tls( T, L, S );

        std::cout << "Centre of reaction: " << tls.centre_of_reaction() << std::endl;

        CrystalLattice crystal_lattice( a,
                                        b,
                                        c,
                                        Angle::angle_90_degrees(),
                                        beta,
                                        Angle::angle_90_degrees() );

        Matrix3D orthogonal_to_fractional = crystal_lattice.for_CASTEP();
        orthogonal_to_fractional.transpose();
        orthogonal_to_fractional.invert();

        // If phi = 0.0 then the normal to the plane of the molecule is parallel to z
        double TLS_C1_x = hCNC.sine() * CN;
        double TLS_C1_y = ( hCNC.sine() / tan_60 ) * CN;
        double TLS_C1_z = 0.0;
        Vector3D TLS_C1( TLS_C1_x, TLS_C1_y, TLS_C1_z );
        double TLS_N1_x = hNCN.sine() * CN;
        double TLS_N1_y = -( hNCN.sine() / tan_60 ) * CN;
        double TLS_N1_z = 0.0;
        Vector3D TLS_N1( TLS_N1_x, TLS_N1_y, TLS_N1_z );
        double TLS_C2_x = 0.0;
        double TLS_C2_y = -( hNCN.sine() / tan_60 ) * CN - hNCN.cosine() * CN;
        double TLS_C2_z = 0.0;
        Vector3D TLS_C2( TLS_C2_x, TLS_C2_y, TLS_C2_z );
        double TLS_N2_x = 0.0;
        double TLS_N2_y = ( hCNC.sine() / tan_60 ) * CN + hCNC.cosine() * CN;
        double TLS_N2_z = 0.0;
        Vector3D TLS_N2( TLS_N2_x, TLS_N2_y, TLS_N2_z );
        Vector3D TLS_H1 = TLS_C1 + TLS_C1 * ( CH / TLS_C1.length() );
        double TLS_H2_x = 0.0;
        double TLS_H2_y = TLS_C2_y - CH;
        double TLS_H2_z = 0.0;
        Vector3D TLS_H2( TLS_H2_x, TLS_H2_y, TLS_H2_z );

        std::vector< double > eigenvalues;
        std::vector< NormalisedVector3D > eigenvectors;

        AnisotropicDisplacementParameters Ucart_C1 = tls.U( crystal_lattice_reduced.fractional_to_orthogonal( orthogonal_to_fractional * TLS_C1 ) );

        std::cout << "C1" << std::endl;
        Ucart_C1.U_cif( crystal_lattice_reduced ).show();
        Ucart_C1.show();
        calculate_eigenvalues( Ucart_C1.U_cart(), eigenvalues, eigenvectors );
        std::cout << "Eigenvalues : " + double2string( eigenvalues[ 0 ] ) + ", " + double2string( eigenvalues[ 1 ] ) + ", " + double2string( eigenvalues[ 2 ] ) << std::endl;
        std::cout << "Eigenvectors:" << std::endl;
        for ( size_t i( 0 ); i != eigenvectors.size(); ++i )
            eigenvectors[i].show();
        std::cout << std::endl;

        AnisotropicDisplacementParameters Ucart_N1 = tls.U( crystal_lattice_reduced.fractional_to_orthogonal( orthogonal_to_fractional * TLS_N1 ) );

        std::cout << "N1" << std::endl;
        Ucart_N1.U_cif( crystal_lattice_reduced ).show();
        Ucart_N1.show();
        calculate_eigenvalues( Ucart_N1.U_cart(), eigenvalues, eigenvectors );
        std::cout << "Eigenvalues : " + double2string( eigenvalues[ 0 ] ) + ", " + double2string( eigenvalues[ 1 ] ) + ", " + double2string( eigenvalues[ 2 ] ) << std::endl;
        std::cout << "Eigenvectors:" << std::endl;
        for ( size_t i( 0 ); i != eigenvectors.size(); ++i )
            eigenvectors[i].show();
        std::cout << std::endl;

        AnisotropicDisplacementParameters Ucart_C2 = tls.U( crystal_lattice_reduced.fractional_to_orthogonal( orthogonal_to_fractional * TLS_C2 ) );

        std::cout << "C2" << std::endl;
        Ucart_C2.U_cif( crystal_lattice_reduced ).show();
        Ucart_C2.show();
        calculate_eigenvalues( Ucart_C2.U_cart(), eigenvalues, eigenvectors );
        std::cout << "Eigenvalues : " + double2string( eigenvalues[ 0 ] ) + ", " + double2string( eigenvalues[ 1 ] ) + ", " + double2string( eigenvalues[ 2 ] ) << std::endl;
        std::cout << "Eigenvectors:" << std::endl;
        for ( size_t i( 0 ); i != eigenvectors.size(); ++i )
            eigenvectors[i].show();
        std::cout << std::endl;

        AnisotropicDisplacementParameters Ucart_N2 = tls.U( crystal_lattice_reduced.fractional_to_orthogonal( orthogonal_to_fractional * TLS_N2 ) );

        std::cout << "N2" << std::endl;
        Ucart_N2.U_cif( crystal_lattice_reduced ).show();
        Ucart_N2.show();
        calculate_eigenvalues( Ucart_N2.U_cart(), eigenvalues, eigenvectors );
        std::cout << "Eigenvalues : " + double2string( eigenvalues[ 0 ] ) + ", " + double2string( eigenvalues[ 1 ] ) + ", " + double2string( eigenvalues[ 2 ] ) << std::endl;
        std::cout << "Eigenvectors:" << std::endl;
        for ( size_t i( 0 ); i != eigenvectors.size(); ++i )
            eigenvectors[i].show();
        std::cout << std::endl;

        AnisotropicDisplacementParameters Ucart_H1 = tls.U( crystal_lattice_reduced.fractional_to_orthogonal( orthogonal_to_fractional * TLS_H1 ) );

        std::cout << "H1" << std::endl;
        Ucart_H1.U_cif( crystal_lattice_reduced ).show();
        Ucart_H1.show();
        calculate_eigenvalues( Ucart_H1.U_cart(), eigenvalues, eigenvectors );
        std::cout << "Eigenvalues : " + double2string( eigenvalues[ 0 ] ) + ", " + double2string( eigenvalues[ 1 ] ) + ", " + double2string( eigenvalues[ 2 ] ) << std::endl;
        std::cout << "Eigenvectors:" << std::endl;
        for ( size_t i( 0 ); i != eigenvectors.size(); ++i )
            eigenvectors[i].show();
        std::cout << std::endl;

        AnisotropicDisplacementParameters Ucart_H2 = tls.U( crystal_lattice_reduced.fractional_to_orthogonal( orthogonal_to_fractional * TLS_H2 ) );

        std::cout << "H2" << std::endl;
        Ucart_H2.U_cif( crystal_lattice_reduced ).show();
        Ucart_H2.show();
        calculate_eigenvalues( Ucart_H2.U_cart(), eigenvalues, eigenvectors );
        std::cout << "Eigenvalues : " + double2string( eigenvalues[ 0 ] ) + ", " + double2string( eigenvalues[ 1 ] ) + ", " + double2string( eigenvalues[ 2 ] ) << std::endl;
        std::cout << "Eigenvectors:" << std::endl;
        for ( size_t i( 0 ); i != eigenvectors.size(); ++i )
            eigenvectors[i].show();
        std::cout << std::endl;

        CrystalStructure crystal_structure;
        crystal_structure.set_crystal_lattice( crystal_lattice );
        crystal_structure.set_space_group( space_group_C2c );
        crystal_structure.add_atom( Atom( Element( "C" ), orthogonal_to_fractional * C1, "C1", Ucart_C1 ) );
        crystal_structure.add_atom( Atom( Element( "N" ), orthogonal_to_fractional * N1, "N1", Ucart_N1 ) );
        crystal_structure.add_atom( Atom( Element( "C" ), orthogonal_to_fractional * C2, "C2", Ucart_C2 ) );
        crystal_structure.add_atom( Atom( Element( "N" ), orthogonal_to_fractional * N2, "N2", Ucart_N2 ) );
        crystal_structure.add_atom( Atom( Element( "H" ), orthogonal_to_fractional * H1, "H1", Ucart_H1 ) );
        crystal_structure.add_atom( Atom( Element( "H" ), orthogonal_to_fractional * H2, "H2", Ucart_H2 ) );
        crystal_structure.save_cif( FileName( "C:\\Users\\jacco\\Documents\\TRIZIN04.cif" ) );
    MACRO_END_GAME

    // Simulate an experimental powder diffraction pattern.
    try
    {
        double wavelength( 1.54056 );
        // ######################## CHANGE THIS ##################################
        double zero_point_error = ZERO_POINT_ERROR; // 0.06;
        Angle two_theta_start( 1.0, Angle::DEGREES );
        Angle two_theta_end(  35.0, Angle::DEGREES );
        Angle two_theta_step( 0.015, Angle::DEGREES );

       CrystalStructure crystal_structure;
            SimulatedPowderPatternCrystalStructure sim_XRPD_crystal_structure( crystal_structure );
            // ######################## CHANGE THIS ##################################
            sim_XRPD_crystal_structure.background_total_signal_normalisation_ = BACKGROUND_TOTAL_SIGNAL_NORMALISATION;
            sim_XRPD_crystal_structure.include_PO_ = true;
            if ( sim_XRPD_crystal_structure.include_PO_ )
            {
                // ######################## CHANGE THIS ##################################
                sim_XRPD_crystal_structure.PO_extent_= PREFERRED_ORIENTATION;
                switch ( crystal_structure.crystal_lattice().lattice_system() )
                {
                    case CrystalLattice::TRICLINIC    :
                    case CrystalLattice::ORTHORHOMBIC : {
                                                            CrystalLattice crystal_lattice = crystal_structure.crystal_lattice();
                                                            if ( crystal_lattice.a() < crystal_lattice.b() )
                                                            {
                                                                if ( crystal_lattice.a() < crystal_lattice.c() )
                                                                    sim_XRPD_crystal_structure.PO_direction_ = MillerIndices( 1, 0, 0 );
                                                            }
                                                            else // b < a
                                                            {
                                                                if ( crystal_lattice.b() < crystal_lattice.c() )
                                                                    sim_XRPD_crystal_structure.PO_direction_ = MillerIndices( 0, 1, 0 );
                                                            }
                                                            break;
                                                        }
                    case CrystalLattice::MONOCLINIC   : {
                                                            sim_XRPD_crystal_structure.PO_direction_ = MillerIndices( 0, 1, 0 );
                                                            break;
                                                        }
                    case CrystalLattice::TRIGONAL     :
                    case CrystalLattice::TETRAGONAL   :
                    case CrystalLattice::HEXAGONAL    : break; // Nothing to do
                    case CrystalLattice::RHOMBOHEDRAL :
                    case CrystalLattice::CUBIC        : sim_XRPD_crystal_structure.include_PO_ = false;
                }
            }
            // ######################## CHANGE THIS ##################################
            sim_XRPD_crystal_structure.FWHM_ = FULL_WIDTH_HALF_MAX; // 0.2; // ?
            PowderPattern result( two_theta_start, two_theta_end, two_theta_step );

            {
                CrystalStructure crystal_structure;
                std::cout << "Now reading cif... " + FileName( "C:\\Users\\jacco\\Documents\\Bicalutamide_FormI_2.cif" ).full_name() << std::endl;
                read_cif( FileName( "C:\\Users\\jacco\\Documents\\Bicalutamide_FormI_2.cif" ), crystal_structure );
                crystal_structure.apply_space_group_symmetry();
                std::cout << "Now calculating powder pattern... " << std::endl;
                PowderPatternCalculator powder_pattern_calculator( crystal_structure );
                powder_pattern_calculator.set_wavelength( wavelength );
                powder_pattern_calculator.set_two_theta_start( two_theta_start );
                powder_pattern_calculator.set_two_theta_end( two_theta_end );
                powder_pattern_calculator.set_two_theta_step( two_theta_step );
                powder_pattern_calculator.set_FWHM( sim_XRPD_crystal_structure.FWHM_ );
                if ( sim_XRPD_crystal_structure.include_PO_ )
                    powder_pattern_calculator.set_preferred_orientation( sim_XRPD_crystal_structure.PO_direction_, sim_XRPD_crystal_structure.PO_extent_ );
                PowderPattern powder_pattern;
                powder_pattern_calculator.calculate( powder_pattern );
                powder_pattern.normalise_total_signal( sim_XRPD_crystal_structure.total_signal_normalisation_ * 3.0 );
                result += powder_pattern;
                if ( sim_XRPD_crystal_structure.include_background_ )
                {
                    PowderPatternCalculator background_powder_pattern_calculator( crystal_structure );
                    background_powder_pattern_calculator.set_wavelength( wavelength );
                    background_powder_pattern_calculator.set_two_theta_start( two_theta_start );
                    background_powder_pattern_calculator.set_two_theta_end( two_theta_end );
                    background_powder_pattern_calculator.set_two_theta_step( two_theta_step );
                    background_powder_pattern_calculator.set_FWHM( 5.0 );
                    // We never include PO for the amorphous background
                    PowderPattern powder_pattern;
                    background_powder_pattern_calculator.calculate( powder_pattern );
                    powder_pattern.normalise_total_signal( sim_XRPD_crystal_structure.background_total_signal_normalisation_ * 3.0 );
                    result += powder_pattern;
                }
            }

            {
                CrystalStructure crystal_structure;
                std::cout << "Now reading cif... " + FileName( "C:\\Users\\jacco\\Documents\\Bicalutamide_FormII_2.cif" ).full_name() << std::endl;
                read_cif( FileName( "C:\\Users\\jacco\\Documents\\Bicalutamide_FormII_2.cif" ), crystal_structure );
                crystal_structure.apply_space_group_symmetry();
                std::cout << "Now calculating powder pattern... " << std::endl;
                PowderPatternCalculator powder_pattern_calculator(crystal_structure );
                powder_pattern_calculator.set_wavelength( wavelength );
                powder_pattern_calculator.set_two_theta_start( two_theta_start );
                powder_pattern_calculator.set_two_theta_end( two_theta_end );
                powder_pattern_calculator.set_two_theta_step( two_theta_step );
                powder_pattern_calculator.set_FWHM( sim_XRPD_crystal_structure.FWHM_ );
                if ( sim_XRPD_crystal_structure.include_PO_ )
                    powder_pattern_calculator.set_preferred_orientation( sim_XRPD_crystal_structure.PO_direction_, sim_XRPD_crystal_structure.PO_extent_ );
                PowderPattern powder_pattern;
                powder_pattern_calculator.calculate( powder_pattern );
                powder_pattern.normalise_total_signal( sim_XRPD_crystal_structure.total_signal_normalisation_ );
                result += powder_pattern;
                if ( sim_XRPD_crystal_structure.include_background_ )
                {
                    PowderPatternCalculator background_powder_pattern_calculator( crystal_structure );
                    background_powder_pattern_calculator.set_wavelength( wavelength );
                    background_powder_pattern_calculator.set_two_theta_start( two_theta_start );
                    background_powder_pattern_calculator.set_two_theta_end( two_theta_end );
                    background_powder_pattern_calculator.set_two_theta_step( two_theta_step );
                    background_powder_pattern_calculator.set_FWHM( 5.0 );
                    // We never include PO for the amorphous background
                    PowderPattern powder_pattern;
                    background_powder_pattern_calculator.calculate( powder_pattern );
                    powder_pattern.normalise_total_signal( sim_XRPD_crystal_structure.background_total_signal_normalisation_ );
                    result += powder_pattern;
                }
            }

            // ######################## CHANGE THIS ##################################
            result.normalise_highest_peak( 10000.0 ); // 300.0
            result.add_constant_background( 20.0 );
            result.recalculate_estimated_standard_deviations();
            result.add_Poisson_noise();
            result.correct_zero_point_error( Angle::from_degrees( -zero_point_error ) );

            result.save_xye( FileName( "C:\\Users\\jacco\\Documents\\Bicalutamide_mixture.xye" ), true );
    MACRO_END_GAME

    try // Write .inp from .cif + two _restraints.txt files + .xye file.
    {
        if ( argc != 3 )
            throw std::runtime_error( "Please give the name of a .cif file and a .xye file that need to be converted to a .inp file." );
        inp_writer( FileName( argv[ 1 ] ), FileName( argv[ 2 ] ) );
    MACRO_END_GAME

    try // Loop over all space groups to test constraints.
    {
        TextFileReader_2 input_file( FileName( "C:\\Users\\jacco\\Documents\\Data\\CCDC\\IT.cif" ) );
        if ( input_file.size() != 8795 )
            throw std::runtime_error( "read_cif(): symmetry line must have same number of items as specified in loop." );
        size_t iLine( 0 );
        std::vector< std::string > words;
        for ( size_t i( 0 ); i != 230; ++i )
        {
            std::string space_group_name = input_file.line( iLine ).substr( 5 );
            space_group_name = space_group_name.substr( 0, space_group_name.find_first_of( "(" ) );
            ++iLine;
            words = split( input_file.line( iLine ) );
            double a = string2double( words[1] );
            ++iLine;
            words = split( input_file.line( iLine ) );
            double b = string2double( words[1] );
            ++iLine;
            words = split( input_file.line( iLine ) );
            double c = string2double( words[1] );
            ++iLine;
            words = split( input_file.line( iLine ) );
            Angle alpha = Angle::from_degrees( string2double( words[1] ) );
            ++iLine;
            words = split( input_file.line( iLine ) );
            Angle beta = Angle::from_degrees( string2double( words[1] ) );
            ++iLine;
            words = split( input_file.line( iLine ) );
            Angle gamma = Angle::from_degrees( string2double( words[1] ) );
            CrystalLattice crystal_lattice( a, b, c, alpha, beta, gamma );
            iLine += 3;
            std::vector< SymmetryOperator > symmetry_operators;
            bool finished( false );
            do // Read the symmetry operators
            {
                words = split( input_file.line( iLine ) );
                if ( ( words[0][0] == '_' ) || ( words[0] == "loop_" ) )
                {
                    finished = true;
                }
                else
                {
                    if ( words.size() != 1 )
                        throw std::runtime_error( "read_cif(): symmetry line must have same number of items as specified in loop." );
                    SymmetryOperator symmetry_operator( words[0] );
                    symmetry_operators.push_back( symmetry_operator );
                    ++iLine;
                }
            } while ( ! finished );
          //  if ( i == 197 )
            {
                SpaceGroup space_group( symmetry_operators, space_group_name );
                CrystalStructure crystal_structure;
                crystal_structure.set_crystal_lattice( crystal_lattice );
                crystal_structure.set_space_group( space_group );
                std::cout << space_group_name << std::endl;
                size_t step_size( 2 );
//                for ( size_t iX( 0 ); iX != step_size; ++iX )
//                {
//                    for ( size_t iY( 0 ); iY != step_size; ++iY )
//                    {
//                        for ( size_t iZ( 0 ); iZ != step_size; ++iZ )
//                        {

                            size_t iX = 1;
                            size_t iY = 0;
                            size_t iZ = 0;

                            Vector3D point( static_cast< double >( iX ) / static_cast< double >( step_size ),
                                            static_cast< double >( iY ) / static_cast< double >( step_size ),
                                            static_cast< double >( iZ ) / static_cast< double >( step_size ) );
                            std::string constraints = write_constraints( crystal_structure, point );
                            std::cout << constraints << std::endl;
//                        }
//                    }
//                }
            }
            iLine += 10;
        }
    MACRO_END_GAME

    try // Add class.
    {
        add_class( "TLS" );
    MACRO_END_GAME

    // Add methyl hydrogen atoms.
    try
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        if ( (true) )
        {
        // The methyl group is attached to atom_1
        std::string atom_1_label( "C4" );
        std::string atom_2_label( "N1" );
        Vector3D atom_1 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_1_label ) ).position() );
        Vector3D atom_2 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_2_label ) ).position() );
        std::vector< Vector3D > methyl_hydrogen_atoms = add_methyl_group( atom_1, atom_2, Angle::from_degrees( 60.0 ) );
        crystal_structure.add_atom( Atom( Element( "H" ), crystal_structure.crystal_lattice().orthogonal_to_fractional( methyl_hydrogen_atoms[0] ), "H" + size_t2string( crystal_structure.natoms() ) ) );
        crystal_structure.add_atom( Atom( Element( "H" ), crystal_structure.crystal_lattice().orthogonal_to_fractional( methyl_hydrogen_atoms[1] ), "H" + size_t2string( crystal_structure.natoms() ) ) );
        crystal_structure.add_atom( Atom( Element( "H" ), crystal_structure.crystal_lattice().orthogonal_to_fractional( methyl_hydrogen_atoms[2] ), "H" + size_t2string( crystal_structure.natoms() ) ) );
        }
        if ( (true) )
        {
        // The methyl group is attached to atom_1
        std::string atom_1_label( "C3" );
        std::string atom_2_label( "N1" );
        Vector3D atom_1 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_1_label ) ).position() );
        Vector3D atom_2 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_2_label ) ).position() );
        std::vector< Vector3D > methyl_hydrogen_atoms = add_methyl_group( atom_1, atom_2, Angle::from_degrees( 0.0 ) );
        crystal_structure.add_atom( Atom( Element( "H" ), crystal_structure.crystal_lattice().orthogonal_to_fractional( methyl_hydrogen_atoms[0] ), "H" + size_t2string( crystal_structure.natoms() ) ) );
        crystal_structure.add_atom( Atom( Element( "H" ), crystal_structure.crystal_lattice().orthogonal_to_fractional( methyl_hydrogen_atoms[1] ), "H" + size_t2string( crystal_structure.natoms() ) ) );
        crystal_structure.add_atom( Atom( Element( "H" ), crystal_structure.crystal_lattice().orthogonal_to_fractional( methyl_hydrogen_atoms[2] ), "H" + size_t2string( crystal_structure.natoms() ) ) );
        }
        if ( (false) )
        {
        // The methyl group is attached to atom_1
        std::string atom_1_label( "C1" );
        std::string atom_2_label( "N2" );
        Vector3D atom_1 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_1_label ) ).position() );
        Vector3D atom_2 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_2_label ) ).position() );
        std::vector< Vector3D > methyl_hydrogen_atoms = add_methyl_group( atom_1, atom_2, Angle::from_degrees( 0.0 ) );
        crystal_structure.add_atom( Atom( Element( "H" ), crystal_structure.crystal_lattice().orthogonal_to_fractional( methyl_hydrogen_atoms[0] ), "H" + size_t2string( crystal_structure.natoms() ) ) );
        crystal_structure.add_atom( Atom( Element( "H" ), crystal_structure.crystal_lattice().orthogonal_to_fractional( methyl_hydrogen_atoms[1] ), "H" + size_t2string( crystal_structure.natoms() ) ) );
        crystal_structure.add_atom( Atom( Element( "H" ), crystal_structure.crystal_lattice().orthogonal_to_fractional( methyl_hydrogen_atoms[2] ), "H" + size_t2string( crystal_structure.natoms() ) ) );
        }
        if ( (true) )
        {
        // The methyl group is attached to atom_1
        std::string atom_1_label( "C1" );
        std::string atom_2_label( "N2" );
        std::string atom_3_label( "C2" );
        Vector3D atom_1 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_1_label ) ).position() );
        Vector3D atom_2 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_2_label ) ).position() );
        Vector3D atom_3 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_3_label ) ).position() );
        std::vector< Vector3D > methyl_hydrogen_atoms = add_methyl_group( atom_1, atom_2, atom_3 );
        crystal_structure.add_atom( Atom( Element( "H" ), crystal_structure.crystal_lattice().orthogonal_to_fractional( methyl_hydrogen_atoms[0] ), "H" + size_t2string( crystal_structure.natoms() ) ) );
        crystal_structure.add_atom( Atom( Element( "H" ), crystal_structure.crystal_lattice().orthogonal_to_fractional( methyl_hydrogen_atoms[1] ), "H" + size_t2string( crystal_structure.natoms() ) ) );
        crystal_structure.add_atom( Atom( Element( "H" ), crystal_structure.crystal_lattice().orthogonal_to_fractional( methyl_hydrogen_atoms[2] ), "H" + size_t2string( crystal_structure.natoms() ) ) );
        }
        if ( (true) )
        {
        // The methyl group is attached to atom_1
        std::string atom_1_label( "C2" );
        std::string atom_2_label( "N2" );
        std::string atom_3_label( "C1" );
        Vector3D atom_1 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_1_label ) ).position() );
        Vector3D atom_2 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_2_label ) ).position() );
        Vector3D atom_3 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_3_label ) ).position() );
        std::vector< Vector3D > methyl_hydrogen_atoms = add_methyl_group( atom_1, atom_2, atom_3 );
        crystal_structure.add_atom( Atom( Element( "H" ), crystal_structure.crystal_lattice().orthogonal_to_fractional( methyl_hydrogen_atoms[0] ), "H" + size_t2string( crystal_structure.natoms() ) ) );
        crystal_structure.add_atom( Atom( Element( "H" ), crystal_structure.crystal_lattice().orthogonal_to_fractional( methyl_hydrogen_atoms[1] ), "H" + size_t2string( crystal_structure.natoms() ) ) );
        crystal_structure.add_atom( Atom( Element( "H" ), crystal_structure.crystal_lattice().orthogonal_to_fractional( methyl_hydrogen_atoms[2] ), "H" + size_t2string( crystal_structure.natoms() ) ) );
        }
        crystal_structure.save_cif( append_to_file_name( input_file_name, "_Me_inserted" ) );
    MACRO_END_GAME

    try // Varenicline.
    {
        TextFileReader_2 input_file( FileName( "C:\\Users\\jacco\\Documents\\FormBC.txt" ) );
        TextFileWriter text_file_writer( FileName( "C:\\Users\\jacco\\Documents\\FormBC_output.txt" ) );
        std::vector< std::string > words;
        size_t iSection( 0 );
        size_t iLine( 0 );
        while ( iLine != input_file.size() )
        {
            words = split( input_file.line( iLine ) );
            if ( words[0] != "TABLE" )
            {
                ++iLine;
                continue;
            }
            ++iSection;
            ++iLine;
            ++iLine;
            if ( iSection == 1 )
            {
                text_file_writer.write_line( "data_FormB" );
                text_file_writer.write_line( "_cell_length_a  7.0753(5)" );
                text_file_writer.write_line( "_cell_length_b  7.7846(5)" );
                text_file_writer.write_line( "_cell_length_c  29.870(2)" );
                text_file_writer.write_line( "_cell_angle_alpha  90" );
                text_file_writer.write_line( "_cell_angle_beta   90" );
                text_file_writer.write_line( "_cell_angle_gamma  90" );
                text_file_writer.write_line( "_symmetry_cell_setting orthorhombic" );
                text_file_writer.write_line( "_symmetry_space_group_name_H-M \"P212121\"" );
                text_file_writer.write_line( "loop_" );
                text_file_writer.write_line( "_space_group_symop_id" );
                text_file_writer.write_line( "_space_group_symop_operation_xyz" );
                text_file_writer.write_line( "1 x,y,z" );
                text_file_writer.write_line( "2 -x+1/2,-y,z+1/2" );
                text_file_writer.write_line( "3 -x,y+1/2,-z+1/2" );
                text_file_writer.write_line( "4 x+1/2,-y+1/2,-z" );
            }
            if ( iSection == 6 )
            {
                text_file_writer.write_line( "data_FormC" );
                text_file_writer.write_line( "_cell_length_a  7.5120" );
                text_file_writer.write_line( "_cell_length_b  29.854" );
                text_file_writer.write_line( "_cell_length_c  7.671" );
                text_file_writer.write_line( "_cell_angle_alpha  90" );
                text_file_writer.write_line( "_cell_angle_beta   90.40" );
                text_file_writer.write_line( "_cell_angle_gamma  90" );
                text_file_writer.write_line( "_symmetry_cell_setting monoclinic" );
                text_file_writer.write_line( "_symmetry_space_group_name_H-M \"P21\"" );
                text_file_writer.write_line( "loop_" );
                text_file_writer.write_line( "_space_group_symop_id" );
                text_file_writer.write_line( "_space_group_symop_operation_xyz" );
                text_file_writer.write_line( "1 x,y,z" );
                text_file_writer.write_line( "2 -x,y+1/2,-z" );
            }
            if ( ( iSection == 1 ) || ( iSection == 6 ) )
            {
                text_file_writer.write_line( "loop_" );
                text_file_writer.write_line( "_atom_site_label" );
                text_file_writer.write_line( "_atom_site_fract_x" );
                text_file_writer.write_line( "_atom_site_fract_y" );
                text_file_writer.write_line( "_atom_site_fract_z" );
                text_file_writer.write_line( "_atom_site_adp_type" );
                text_file_writer.write_line( "_atom_site_U_iso_or_equiv" );
                // x y z U(eq)
                // N(1) 8211(8) 10638(7) 12233(1) 61(1)
                std::vector< std::string > data;
                while ( ( iLine != input_file.size() ) && ( input_file.line( iLine ).substr( 0, 5 ) != "TABLE" ) )
                {
                    words = split( input_file.line( iLine ) );
                    for ( size_t i( 0 ); i != words.size(); ++i )
                        data.push_back( words[i] );
                    ++iLine;
                }
                if ( ( data.size() % 5 ) != 0 )
                    throw std::runtime_error( "Number of data items inconsistent in Section 1 or 6." );
                for ( size_t i( 0 ); i != data.size(); i += 5 )
                {
                    std::string atom_label = data[i];
                    atom_label = remove( atom_label, '(' );
                    atom_label = remove( atom_label, ')' );
                    std::string x = clean_string( data[i+1], 4 );
                    std::string y = clean_string( data[i+2], 4 );
                    std::string z = clean_string( data[i+3], 4 );
                    std::string U = clean_string( data[i+4], 3 );
                    text_file_writer.write_line( pad( atom_label, 4 ) + " " + x + " " + y + " " + z + " " + "Uani" + " " + U );
                }
            }
            if ( ( iSection == 2 ) || ( iSection == 7 ) )
            {
                // x y z U(eq)
                // H(2A) 10149  8958 12367 80
                std::vector< std::string > data;
                while ( ( iLine != input_file.size() ) && ( input_file.line( iLine ).substr( 0, 5 ) != "TABLE" ) )
                {
                    words = split( input_file.line( iLine ) );
                    for ( size_t i( 0 ); i != words.size(); ++i )
                        data.push_back( words[i] );
                    ++iLine;
                }
                if ( ( data.size() % 5 ) != 0 )
                    throw std::runtime_error( "Number of data items inconsistent in Section 2 or 7." );
                for ( size_t i( 0 ); i != data.size(); i += 5 )
                {
                    std::string atom_label = data[i];
                    atom_label = remove( atom_label, '(' );
                    atom_label = remove( atom_label, ')' );
                    std::string x = clean_string( data[i+1], 4 );
                    std::string y = clean_string( data[i+2], 4 );
                    std::string z = clean_string( data[i+3], 4 );
                    std::string U = clean_string( data[i+4], 3 );
                    text_file_writer.write_line( pad( atom_label, 4 ) + " " + x + " " + y + " " + z + " " + "Uiso" + " " + U );
                }
            }
            if ( ( iSection == 3 ) || ( iSection == 8 ) )
            {
                text_file_writer.write_line( "loop_" );
                text_file_writer.write_line( "_atom_site_aniso_label" );
                text_file_writer.write_line( "_atom_site_aniso_U_11" );
                text_file_writer.write_line( "_atom_site_aniso_U_22" );
                text_file_writer.write_line( "_atom_site_aniso_U_33" );
                text_file_writer.write_line( "_atom_site_aniso_U_23" );
                text_file_writer.write_line( "_atom_site_aniso_U_13" );
                text_file_writer.write_line( "_atom_site_aniso_U_12" );
                // U11 U22 U33 U23 U13 U12
                // N(1) 63(4) 70(4) 50(3) 12(2) -2(3) 8(3)
                std::vector< std::string > data;
                while ( ( iLine != input_file.size() ) && ( input_file.line( iLine ).substr( 0, 5 ) != "TABLE" ) )
                {
                    words = split( input_file.line( iLine ) );
                    for ( size_t i( 0 ); i != words.size(); ++i )
                        data.push_back( words[i] );
                    ++iLine;
                }
                if ( ( data.size() % 7 ) != 0 )
                    throw std::runtime_error( "Number of data items inconsistent in Section 3 or 8." );
                for ( size_t i( 0 ); i != data.size(); i += 7 )
                {
                    std::string atom_label = data[i];
                    atom_label = remove( atom_label, '(' );
                    atom_label = remove( atom_label, ')' );
                    std::string U11 = clean_string( data[i+1], 3 );
                    std::string U22 = clean_string( data[i+2], 3 );
                    std::string U33 = clean_string( data[i+3], 3 );
                    std::string U23 = clean_string( data[i+4], 3 );
                    std::string U13 = clean_string( data[i+5], 3 );
                    std::string U12 = clean_string( data[i+6], 3 );
                    text_file_writer.write_line( pad( atom_label, 4 ) + " " + U11 + " " + U22 + " " + U33 + " " + U23 + " " + U13 + " " + U12 );
                }
            }
            if ( ( iSection == 4 ) || ( iSection == 9 ) )
            {
                // N(1)-C(2) 1.316(6)
                text_file_writer.write_line( "loop_" );
                text_file_writer.write_line( "  _geom_bond_atom_site_label_1" );
                text_file_writer.write_line( "  _geom_bond_atom_site_label_2" );
                text_file_writer.write_line( "  _geom_bond_distance" );
                std::vector< std::string > data;
                while ( ( iLine != input_file.size() ) && ( input_file.line( iLine ).substr( 0, 5 ) != "TABLE" ) )
                {
                    words = split( input_file.line( iLine ) );
                    for ( size_t i( 0 ); i != words.size(); ++i )
                        data.push_back( words[i] );
                    ++iLine;
                }
                if ( ( data.size() % 2 ) != 0 )
                    throw std::runtime_error( "Number of data items inconsistent in Section 4 or 9." );
                for ( size_t i( 0 ); i != data.size(); i += 2 )
                {
                    std::string atom_labels = data[i];
                    atom_labels = remove( atom_labels, '(' );
                    atom_labels = remove( atom_labels, ')' );
                    words = split( atom_labels, '-' );
                //    for ( size_t j( 0 ); j != words.size(); ++ j )
              //          std::cout << words[j] << " ";
                   // std::cout << std::endl;
                    if ( words.size() != 2 )
                        throw std::runtime_error( "Number of atom labels must be two." );
                    std::string length = data[i+1];
                    text_file_writer.write_line( pad( words[0], 4 ) + " " + pad( words[1], 4 ) + " " + length );
                }
            }
            if ( ( iSection == 5 ) || ( iSection == 10 ) )
            {
                // C(2)-N(1)-C(6) 115.0(5)
                text_file_writer.write_line( "loop_" );
                text_file_writer.write_line( "  _geom_angle_atom_site_label_1" );
                text_file_writer.write_line( "  _geom_angle_atom_site_label_2" );
                text_file_writer.write_line( "  _geom_angle_atom_site_label_3" );
                text_file_writer.write_line( "  _geom_angle" );
                std::vector< std::string > data;
                while ( ( iLine != input_file.size() ) && ( input_file.line( iLine ).substr( 0, 5 ) != "TABLE" ) )
                {
                    words = split( input_file.line( iLine ) );
                    for ( size_t i( 0 ); i != words.size(); ++i )
                        data.push_back( words[i] );
                    ++iLine;
                }
                if ( ( data.size() % 2 ) != 0 )
                    throw std::runtime_error( "Number of data items inconsistent in Section 5 or 10." );
                for ( size_t i( 0 ); i != data.size(); i += 2 )
                {
                    std::string atom_labels = data[i];
                    atom_labels = remove( atom_labels, '(' );
                    atom_labels = remove( atom_labels, ')' );
                    words = split( atom_labels, '-' );
                    if ( words.size() != 3 )
                        throw std::runtime_error( "Number of atom labels must be three." );
                    std::string length = data[i+1];
                    text_file_writer.write_line( pad( words[0], 4 ) + " " + pad( words[1], 4 ) + " " + pad( words[2], 4 ) + " " + length );
                }
            }
        }
    MACRO_END_GAME

    try // Insert one hydrogen atom between two atoms.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT

        if ( (true) )
        {
        std::string origin_atom_label( "O31" );
        std::string neighbour_atom_label( "O38" );
        Vector3D origin_atom_frac = crystal_structure.atom( crystal_structure.find_label( origin_atom_label ) ).position();
        Vector3D neighbour_atom_frac = crystal_structure.atom( crystal_structure.find_label( neighbour_atom_label ) ).position();
        // Find shortest disance between the two taking periodicity and space-group symmetry into account.
        double distance;
        Vector3D difference_frac;
        crystal_structure.shortest_distance( origin_atom_frac, neighbour_atom_frac, distance, difference_frac );
        if ( distance > 3.0 )
            std::cout << "Warning: distance > 3.0 A" << std::endl;
        // We need to be able to set the length of the X-H bond, so we must work in Cartesian coordinates.
        Vector3D difference_cart = crystal_structure.crystal_lattice().fractional_to_orthogonal( difference_frac );
        double target_bond_length( 1.0 );
        if ( crystal_structure.atom( crystal_structure.find_label( origin_atom_label ) ).element() == Element( "C" ) )
            target_bond_length = 1.089;
        else if ( crystal_structure.atom( crystal_structure.find_label( origin_atom_label ) ).element() == Element( "N" ) )
            target_bond_length = 1.015;
        else if ( crystal_structure.atom( crystal_structure.find_label( origin_atom_label ) ).element() == Element( "O" ) )
            target_bond_length = 0.993;
        difference_cart *= target_bond_length / distance;
        difference_frac = crystal_structure.crystal_lattice().orthogonal_to_fractional( difference_cart );
        Vector3D H_atom_frac = origin_atom_frac + difference_frac;
        crystal_structure.add_atom( Atom( Element( "H" ), H_atom_frac, "H" + size_t2string( crystal_structure.natoms() ) ) );
        }

        if ( (false) )
        {
        std::string origin_atom_label( "O41" );
        std::string neighbour_atom_label( "O40" );
        Vector3D origin_atom_frac = crystal_structure.atom( crystal_structure.find_label( origin_atom_label ) ).position();
        Vector3D neighbour_atom_frac = crystal_structure.atom( crystal_structure.find_label( neighbour_atom_label ) ).position();
        // Find shortest disance between the two taking periodicity and space-group symmetry into account.
        double distance;
        Vector3D difference_frac;
        crystal_structure.shortest_distance( origin_atom_frac, neighbour_atom_frac, distance, difference_frac );
        if ( distance > 3.0 )
            std::cout << "Warning: distance > 3.0 A" << std::endl;
        // We need to be able to set the length of the X-H bond, so we must work in Cartesian coordinates.
        Vector3D difference_cart = crystal_structure.crystal_lattice().fractional_to_orthogonal( difference_frac );
        double target_bond_length( 1.0 );
        if ( crystal_structure.atom( crystal_structure.find_label( origin_atom_label ) ).element() == Element( "C" ) )
            target_bond_length = 1.089;
        else if ( crystal_structure.atom( crystal_structure.find_label( origin_atom_label ) ).element() == Element( "N" ) )
            target_bond_length = 1.015;
        else if ( crystal_structure.atom( crystal_structure.find_label( origin_atom_label ) ).element() == Element( "O" ) )
            target_bond_length = 0.993;
        difference_cart *= target_bond_length / distance;
        difference_frac = crystal_structure.crystal_lattice().orthogonal_to_fractional( difference_cart );
        Vector3D H_atom_frac = origin_atom_frac + difference_frac;
        crystal_structure.add_atom( Atom( Element( "H" ), H_atom_frac, "H" + size_t2string( crystal_structure.natoms() ) ) );
        }

        if ( (false) )
        {
        std::string origin_atom_label( "O41" );
        std::string neighbour_atom_label( "O40" );
        Vector3D origin_atom_frac = crystal_structure.atom( crystal_structure.find_label( origin_atom_label ) ).position();
        Vector3D neighbour_atom_frac = crystal_structure.atom( crystal_structure.find_label( neighbour_atom_label ) ).position();
        // Find shortest disance between the two taking periodicity and space-group symmetry into account.
        double distance;
        Vector3D difference_frac;
        crystal_structure.second_shortest_distance( origin_atom_frac, neighbour_atom_frac, distance, difference_frac );
        if ( distance > 3.0 )
            std::cout << "Warning: distance > 3.0 A" << std::endl;
        // We need to be able to set the length of the X-H bond, so we must work in Cartesian coordinates.
        Vector3D difference_cart = crystal_structure.crystal_lattice().fractional_to_orthogonal( difference_frac );
        double target_bond_length( 1.0 );
        if ( crystal_structure.atom( crystal_structure.find_label( origin_atom_label ) ).element() == Element( "C" ) )
            target_bond_length = 1.089;
        else if ( crystal_structure.atom( crystal_structure.find_label( origin_atom_label ) ).element() == Element( "N" ) )
            target_bond_length = 1.015;
        else if ( crystal_structure.atom( crystal_structure.find_label( origin_atom_label ) ).element() == Element( "O" ) )
            target_bond_length = 0.993;
        difference_cart *= target_bond_length / distance;
        difference_frac = crystal_structure.crystal_lattice().orthogonal_to_fractional( difference_cart );
        Vector3D H_atom_frac = origin_atom_frac + difference_frac;
        crystal_structure.add_atom( Atom( Element( "H" ), H_atom_frac, "H" + size_t2string( crystal_structure.natoms() ) ) );
        }

        if ( (false) )
        {
        std::string origin_atom_label( "O2" );
        std::string neighbour_atom_label( "N2" );
        Vector3D origin_atom_frac = crystal_structure.atom( crystal_structure.find_label( origin_atom_label ) ).position();
        Vector3D neighbour_atom_frac = crystal_structure.atom( crystal_structure.find_label( neighbour_atom_label ) ).position();
        // Find shortest disance between the two taking periodicity and space-group symmetry into account.
        double distance;
        Vector3D difference_frac;
        crystal_structure.shortest_distance( origin_atom_frac, neighbour_atom_frac, distance, difference_frac );
        if ( distance > 3.0 )
            std::cout << "Warning: distance > 3.0 A" << std::endl;
        // We need to be able to set the length of the X-H bond, so we must work in Cartesian coordinates.
        Vector3D difference_cart = crystal_structure.crystal_lattice().fractional_to_orthogonal( difference_frac );
        double target_bond_length( 1.0 );
        if ( crystal_structure.atom( crystal_structure.find_label( origin_atom_label ) ).element() == Element( "C" ) )
            target_bond_length = 1.089;
        else if ( crystal_structure.atom( crystal_structure.find_label( origin_atom_label ) ).element() == Element( "N" ) )
            target_bond_length = 1.015;
        else if ( crystal_structure.atom( crystal_structure.find_label( origin_atom_label ) ).element() == Element( "O" ) )
            target_bond_length = 0.993;
        difference_cart *= target_bond_length / distance;
        difference_frac = crystal_structure.crystal_lattice().orthogonal_to_fractional( difference_cart );
        Vector3D H_atom_frac = origin_atom_frac + difference_frac;
        crystal_structure.add_atom( Atom( Element( "H" ), H_atom_frac, "H" + size_t2string( crystal_structure.natoms() ) ) );
        }

        crystal_structure.save_cif( append_to_file_name( input_file_name, "_H_inserted" ) );
    MACRO_END_GAME

    // add_hydrogen_atom_to_sp3_atom().
    try
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        std::string central_atom_label( "C35" );
        std::string neighbour_1_label( "O40" );
        std::string neighbour_2_label( "C36" );
        std::string neighbour_3_label( "C33" );
        Vector3D central_atom = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( central_atom_label ) ).position() );
        Vector3D neighbour_1 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( neighbour_1_label ) ).position() );
        Vector3D neighbour_2 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( neighbour_2_label ) ).position() );
        Vector3D neighbour_3 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( neighbour_3_label ) ).position() );
        Vector3D hydrogen_atom = add_hydrogen_atom_to_sp3_atom( central_atom, Element( "C" ), neighbour_1, neighbour_2, neighbour_3 );
        crystal_structure.add_atom( Atom( Element( "H" ), crystal_structure.crystal_lattice().orthogonal_to_fractional( hydrogen_atom ), "H" + size_t2string( crystal_structure.natoms() ) ) );
        crystal_structure.save_cif( append_to_file_name( input_file_name, "_Me_inserted" ) );
    MACRO_END_GAME

    try // Sudoku.
    {
        if ( false ) // Easy
        {
            std::vector< std::string > sudoku_string;
            sudoku_string.push_back( "900428000" );
            sudoku_string.push_back( "486000209" );
            sudoku_string.push_back( "302090854" );
            sudoku_string.push_back( "108000730" );
            sudoku_string.push_back( "009374010" );
            sudoku_string.push_back( "007850090" );
            sudoku_string.push_back( "820000906" );
            sudoku_string.push_back( "090006107" );
            sudoku_string.push_back( "060189000" );
            Sudoku sudoku( sudoku_string );
            std::cout << "###############################################" << std::endl;
            sudoku.show();
            Sudoku solved_sudoku = solve( sudoku );
            solved_sudoku.show();
        }
        if ( false ) // Medium
        {
            std::vector< std::string > sudoku_string;
            sudoku_string.push_back( "090000713" );
            sudoku_string.push_back( "000687000" );
            sudoku_string.push_back( "720039000" );
            sudoku_string.push_back( "060200537" );
            sudoku_string.push_back( "952700000" );
            sudoku_string.push_back( "000510020" );
            sudoku_string.push_back( "109046000" );
            sudoku_string.push_back( "000001409" );
            sudoku_string.push_back( "405000108" );
            Sudoku sudoku( sudoku_string );
            std::cout << "###############################################" << std::endl;
            Sudoku solved_sudoku = solve( sudoku );
            solved_sudoku.show();
        }
        if ( false ) // Difficult
        {
            std::vector< std::string > sudoku_string;
            sudoku_string.push_back( "090000867" );
            sudoku_string.push_back( "000010302" );
            sudoku_string.push_back( "603009000" );
            sudoku_string.push_back( "040700000" );
            sudoku_string.push_back( "000300081" );
            sudoku_string.push_back( "020805000" );
            sudoku_string.push_back( "908040070" );
            sudoku_string.push_back( "004080500" );
            sudoku_string.push_back( "006900010" );
            Sudoku sudoku( sudoku_string );
            std::cout << "###############################################" << std::endl;
            Sudoku solved_sudoku = solve( sudoku );
            solved_sudoku.show();
        }
        if ( false ) // Meister 1
        {
            std::vector< std::string > sudoku_string;
            sudoku_string.push_back( "200000008" );
            sudoku_string.push_back( "053020970" );
            sudoku_string.push_back( "070090030" );
            sudoku_string.push_back( "000007010" );
            sudoku_string.push_back( "700000004" );
            sudoku_string.push_back( "080600000" );
            sudoku_string.push_back( "020070050" );
            sudoku_string.push_back( "031050890" );
            sudoku_string.push_back( "500000002" );
            Sudoku sudoku( sudoku_string );
            std::cout << "###############################################" << std::endl;
            sudoku.show();
            Sudoku solved_sudoku = solve( sudoku );
            solved_sudoku.show();
        }
        if ( false ) // Meister 2
        {
            std::vector< std::string > sudoku_string;
            sudoku_string.push_back( "010753040" );
            sudoku_string.push_back( "070894050" );
            sudoku_string.push_back( "400612008" );
            sudoku_string.push_back( "007461920" );
            sudoku_string.push_back( "094527080" );
            sudoku_string.push_back( "201938400" );
            sudoku_string.push_back( "100275004" );
            sudoku_string.push_back( "020140030" );
            sudoku_string.push_back( "040380010" );
            Sudoku sudoku( sudoku_string );
            std::cout << "###############################################" << std::endl;
            Sudoku solved_sudoku = solve( sudoku );
            solved_sudoku.show();
        }
        if ( false ) // Schwer
        {
            std::vector< std::string > sudoku_string;
            sudoku_string.push_back( "006000082" );
            sudoku_string.push_back( "100400300" );
            sudoku_string.push_back( "003290106" );
            sudoku_string.push_back( "089000000" );
            sudoku_string.push_back( "054013000" );
            sudoku_string.push_back( "000000073" );
            sudoku_string.push_back( "260030800" );
            sudoku_string.push_back( "000005200" );
            sudoku_string.push_back( "000004705" );
            Sudoku sudoku( sudoku_string );
            std::cout << "###############################################" << std::endl;
            Sudoku solved_sudoku = solve( sudoku );
            solved_sudoku.show();
        }
        if ( false ) // Toughest
        {
            std::vector< std::string > sudoku_string;
            sudoku_string.push_back( "800000000" );
            sudoku_string.push_back( "003600000" );
            sudoku_string.push_back( "070090200" );
            sudoku_string.push_back( "050007000" );
            sudoku_string.push_back( "000045700" );
            sudoku_string.push_back( "000100030" );
            sudoku_string.push_back( "001000068" );
            sudoku_string.push_back( "008500010" );
            sudoku_string.push_back( "090000400" );
            Sudoku sudoku( sudoku_string );
            std::cout << "###############################################" << std::endl;
            Sudoku solved_sudoku = solve( sudoku );
            solved_sudoku.show();
        }
        if ( false ) // NYT
        {
            std::vector< std::string > sudoku_string;
            sudoku_string.push_back( "080020561" );
            sudoku_string.push_back( "000100007" );
            sudoku_string.push_back( "000500000" );
            sudoku_string.push_back( "050090408" );
            sudoku_string.push_back( "007850003" );
            sudoku_string.push_back( "090010050" );
            sudoku_string.push_back( "204001805" );
            sudoku_string.push_back( "060085000" );
            sudoku_string.push_back( "000200100" );
            Sudoku sudoku( sudoku_string );
            std::cout << "###############################################" << std::endl;
            Sudoku solved_sudoku = solve( sudoku );
            solved_sudoku.show();
        }
        if ( false ) // Learn something (Empty Rectangle)
        {
            std::vector< std::string > sudoku_string;
            sudoku_string.push_back( "800006305" );
            sudoku_string.push_back( "040000070" );
            sudoku_string.push_back( "000000000" );
            sudoku_string.push_back( "010038704" );
            sudoku_string.push_back( "000104000" );
            sudoku_string.push_back( "300070290" );
            sudoku_string.push_back( "000003000" );
            sudoku_string.push_back( "020000040" );
            sudoku_string.push_back( "506800002" );
            Sudoku sudoku( sudoku_string );
            std::cout << "###############################################" << std::endl;
            Sudoku solved_sudoku = solve( sudoku );
            solved_sudoku.show();
        }
        if ( false ) // Learn something (X-Wing)
        {
            std::vector< std::string > sudoku_string;
            sudoku_string.push_back( "600090007" );
            sudoku_string.push_back( "040007100" );
            sudoku_string.push_back( "002800050" );
            sudoku_string.push_back( "800000090" );
            sudoku_string.push_back( "000070000" );
            sudoku_string.push_back( "030000008" );
            sudoku_string.push_back( "050002300" );
            sudoku_string.push_back( "004500020" );
            sudoku_string.push_back( "900030004" );
            Sudoku sudoku( sudoku_string );
            std::cout << "###############################################" << std::endl;
            Sudoku solved_sudoku = solve( sudoku );
            solved_sudoku.show();
        }
        if ( true ) // AI Escargot 6313 guesses, stack pointer = 1
        {
            std::vector< std::string > sudoku_string;
            sudoku_string.push_back( "100007090" );
            sudoku_string.push_back( "030020008" );
            sudoku_string.push_back( "009600500" );
            sudoku_string.push_back( "005300900" );
            sudoku_string.push_back( "010080002" );
            sudoku_string.push_back( "600004000" );
            sudoku_string.push_back( "300000010" );
            sudoku_string.push_back( "040000007" );
            sudoku_string.push_back( "007000300" );
            Sudoku sudoku( sudoku_string );
            std::cout << "###############################################" << std::endl;
            Sudoku solved_sudoku = solve( sudoku );
            solved_sudoku.show();
        }
    MACRO_END_GAME

    try // BFDH.
    {

        //calculate_BFDH( crystal_lattice );
    MACRO_END_GAME

    try // RMSD
    {
        {
        std::vector< double > list_1;
        std::vector< double > list_2;
        list_1.push_back( 11.0 );
        list_2.push_back(  1.1 );
        list_1.push_back( 12.0 );
        list_2.push_back(  1.9 );
        list_1.push_back( 13.0 );
        list_2.push_back(  3.1 );
        list_1.push_back( 14.0 );
        list_2.push_back(  3.9 );
        double RMSD = calculate_sample_RMSD( list_1, list_2 );
        std::cout << RMSD << std::endl;
        CorrelationMatrix differences_matrix_1( list_1.size() );
        CorrelationMatrix differences_matrix_2( list_1.size() );
        differences_matrix_1.set_value_on_diagonal( 0.0 );
        differences_matrix_2.set_value_on_diagonal( 0.0 );
        for ( size_t i( 0 ); i != list_1.size(); ++i )
        {
            for ( size_t j( i ); j != list_1.size(); ++j )
            {
                differences_matrix_1.set_value( i, j, list_1[i] - list_1[j] );
                differences_matrix_2.set_value( i, j, list_2[i] - list_2[j] );
            }
        }
        RMSD = 0.0;
        size_t N( 0 );
        for ( size_t i( 0 ); i != list_1.size(); ++i )
        {
            for ( size_t j( i ); j != list_1.size(); ++j )
            {
                RMSD += square( differences_matrix_2.value( i, j ) - differences_matrix_1.value( i, j ) );
                ++N;
            }
        }
        RMSD /= N;
        RMSD = std::sqrt( RMSD );
        std::cout << RMSD << std::endl;
        }
        {
        std::vector< double > list_1;
        std::vector< double > list_2;
        RandomNumberGenerator_double RNG;
        for ( size_t i( 0 ); i != 1000; ++i )
        {
            list_1.push_back( RNG.next_number() - 145.0 + i );
            list_2.push_back( RNG.next_number() - 4.5   + i );
        }
        double RMSD = calculate_sample_RMSD( list_1, list_2 );
        std::cout << RMSD << std::endl;
        CorrelationMatrix differences_matrix_1( list_1.size() );
        CorrelationMatrix differences_matrix_2( list_1.size() );
        differences_matrix_1.set_value_on_diagonal( 0.0 );
        differences_matrix_2.set_value_on_diagonal( 0.0 );
        for ( size_t i( 0 ); i != list_1.size(); ++i )
        {
            for ( size_t j( i ); j != list_1.size(); ++j )
            {
                differences_matrix_1.set_value( i, j, list_1[i] - list_1[j] );
                differences_matrix_2.set_value( i, j, list_2[i] - list_2[j] );
            }
        }
        RMSD = 0.0;
        size_t N( 0 );
        for ( size_t i( 0 ); i != list_1.size(); ++i )
        {
            for ( size_t j( i ); j != list_1.size(); ++j )
            {
                RMSD += square( differences_matrix_2.value( i, j ) - differences_matrix_1.value( i, j ) );
                ++N;
            }
        }
        RMSD /= N;
        RMSD = std::sqrt( RMSD );
        std::cout << RMSD << std::endl;
        }
    MACRO_END_GAME

    try // Farey
    {
        if ( (false) )
        {
            RandomNumberGenerator_double RNG;
            double maximum_difference( 0.0 );
            for ( size_t i( 0 ); i != 1000000; ++i )
            {
                double target = RNG.next_number();
                Fraction result = Farey( target, 1000 );
                double difference = std::abs( result.to_double() - target );
                if ( difference > maximum_difference )
                    maximum_difference = difference;
              //  std::cout << "target = " << target << ", fraction = " << result.to_string() << " = " << result.to_double() << std::endl;
            }
            std::cout << "maximum_difference = " << maximum_difference << std::endl;
        }
        {
            double target = CONSTANT_PI;
            Fraction result = Farey_terminate_on_error( target, 0.0000005 );
            std::cout << "target = " << target << ", fraction = " << result.to_string() << " = " << result.to_double() << std::endl;
        }
        {
            double target = CONSTANT_PI;
            Fraction result = Farey( target, 1000 );
            std::cout << "target = " << target << ", fraction = " << result.to_string() << " = " << result.to_double() << std::endl;
        }
        {
            double target = CONSTANT_PI;
            Fraction result = Farey( target, 10000 );
            std::cout << "target = " << target << ", fraction = " << result.to_string() << " = " << result.to_double() << std::endl;
        }
        {
            double target = 0.167;
            Fraction result = Farey_terminate_on_error( target, 0.0005 );
            std::cout << "target = " << target << ", fraction = " << result.to_string() << " = " << result.to_double() << std::endl;
        }
        {
            double target = 0.333;
            Fraction result = Farey_terminate_on_error( target, 0.0005 );
            std::cout << "target = " << target << ", fraction = " << result.to_string() << " = " << result.to_double() << std::endl;
        }
        {
            double target = 0.605551;
            Fraction result = Farey_terminate_on_error( target, 0.0000005 );
            std::cout << "target = " << target << ", fraction = " << result.to_string() << " = " << result.to_double() << std::endl;
        }
        {
            double target = -0.605551;
            Fraction result = Farey( target, 100 ); // 3/5 is wrong, it should be 17/28
            std::cout << "target = " << target << ", fraction = " << result.to_string() << " = " << result.to_double() << std::endl;
        }
        {
            double target = 0.5;
            Fraction result = Farey( target, 1000 );
            std::cout << "target = " << target << ", fraction = " << result.to_string() << " = " << result.to_double() << std::endl;
        }
        {
            double target = 0.2;
            Fraction result = Farey( target, 1000 );
            std::cout << "target = " << target << ", fraction = " << result.to_string() << " = " << result.to_double() << std::endl;
        }
        {
            double target = 1.0/3.0;
            Fraction result = Farey( target, 1000 );
            std::cout << "target = " << target << ", fraction = " << result.to_string() << " = " << result.to_double() << std::endl;
        }
        {
            double target = 0.25;
            Fraction result = Farey( target, 1 );
            std::cout << "target = " << target << ", fraction = " << result.to_string() << " = " << result.to_double() << std::endl;
        }
        {
            double target = 0.75;
            Fraction result = Farey( target, 1 );
            std::cout << "target = " << target << ", fraction = " << result.to_string() << " = " << result.to_double() << std::endl;
        }
        {
            double target = 0.95;
            Fraction result = Farey( target, 1000 );
            std::cout << "target = " << target << ", fraction = " << result.to_string() << " = " << result.to_double() << std::endl;
        }

    MACRO_END_GAME

    try // Loop over all space groups
    {
        RandomNumberGenerator_integer rng;
        ReflectionList reflection_list;
        for ( size_t i( 0 ); i != 1; ++i )
            reflection_list.push_back( MillerIndices( rng.next_number( -100, 100 ), rng.next_number( -100, 100 ), rng.next_number( -100, 100 ) ), 1.0, 1.0, 2 );
        std::vector< MillerIndices > PO_directions;
        if ( (false) )
        {
            for ( size_t i( 0 ); i != 3; ++i )
            {
                MillerIndices miller_indices( rng.next_number( -5, 5 ), rng.next_number( -5, 5 ), rng.next_number( -5, 5 ) );
                miller_indices.make_relative_prime();
                PO_directions.push_back( miller_indices );
            }
        }
        else
        {
            int max_index( 1 );
            for ( int i( -max_index ); i != max_index+1; ++i )
            {
                for ( int j( -max_index ); j != max_index+1; ++j )
                {
                    for ( int k( -max_index ); k != max_index+1; ++k )
                    {
                        MillerIndices miller_indices( i, j, k );
                        if ( miller_indices.is_000() )
                            continue;
                        miller_indices.make_relative_prime();
                        PO_directions.push_back( miller_indices );
                    }
                }
            }
        }
        TextFileReader_2 input_file( FileName( "C:\\Users\\jacco\\Documents\\Data\\CCDC\\IT.cif" ) );
        if ( input_file.size() != 8795 )
            throw std::runtime_error( "read_cif(): symmetry line must have same number of items as specified in loop." );
        size_t iLine( 0 );
        std::vector< std::string > words;
        std::vector< double > r_values;
        r_values.push_back( 0.7 );
//        r_values.push_back( 0.8 );
//        r_values.push_back( 1.2 );
//        r_values.push_back( 1.3 );
        for ( size_t i( 0 ); i != 230; ++i )
        {
            std::string space_group_name = input_file.line( iLine ).substr( 5 );
            space_group_name = space_group_name.substr( 0, space_group_name.find_first_of( "(" ) );
            ++iLine;
            words = split( input_file.line( iLine ) );
            double a = string2double( words[1] );
            ++iLine;
            words = split( input_file.line( iLine ) );
            double b = string2double( words[1] );
            ++iLine;
            words = split( input_file.line( iLine ) );
            double c = string2double( words[1] );
            ++iLine;
            words = split( input_file.line( iLine ) );
            Angle alpha = Angle::from_degrees( string2double( words[1] ) );
            ++iLine;
            words = split( input_file.line( iLine ) );
            Angle beta = Angle::from_degrees( string2double( words[1] ) );
            ++iLine;
            words = split( input_file.line( iLine ) );
            Angle gamma = Angle::from_degrees( string2double( words[1] ) );
            CrystalLattice crystal_lattice( a, b, c, alpha, beta, gamma );
            iLine += 3;
            std::vector< SymmetryOperator > symmetry_operators;
            bool finished( false );
            do // Read the symmetry operators
            {
                words = split( input_file.line( iLine ) );
                if ( ( words[0][0] == '_' ) || ( words[0] == "loop_" ) )
                {
                    finished = true;
                }
                else
                {
                    if ( words.size() != 1 )
                        throw std::runtime_error( "read_cif(): symmetry line must have same number of items as specified in loop." );
                    SymmetryOperator symmetry_operator( words[0] );
                    symmetry_operators.push_back( symmetry_operator );
                    ++iLine;
                }
            } while ( ! finished );
            SpaceGroup space_group( symmetry_operators, space_group_name );
            PointGroup laue_class = space_group.laue_class();
            for ( size_t iPOD( 0 ); iPOD != PO_directions.size(); ++iPOD )
            {
                Vector3D PO_vector = reciprocal_lattice_point( PO_directions[iPOD], crystal_lattice );
                for ( size_t iMD_par( 0 ); iMD_par != r_values.size(); ++iMD_par )
                {
                    double r = r_values[iMD_par];
                    for ( size_t iReflection( 0 ); iReflection != reflection_list.size(); ++iReflection )
                    {
                        // Calculate equivalent reflections and check that the March-Dollase PO corrections are the same for all of them
                        // Now check that the March-Dollase PO corrections are the same for all of them
                        Vector3D H = reciprocal_lattice_point( reflection_list.miller_indices( iReflection ), crystal_lattice );
                        Angle alpha = angle( PO_vector, H );
                        double reference_PO = std::pow( square(r) * square(alpha.cosine()) + square(alpha.sine())/r, -3.0/2.0 );
                        for ( size_t j( 0 ); j != laue_class.nsymmetry_operators(); ++j )
                        {
                            MillerIndices equivalent_reflection = reflection_list.miller_indices( iReflection ) * laue_class.symmetry_operator( j );
                            // Now check that the March-Dollase PO corrections are the same for all of them
                            Vector3D H = reciprocal_lattice_point( equivalent_reflection, crystal_lattice );
                            Angle alpha = angle( PO_vector, H );
                            double current_PO = std::pow( square(r) * square(alpha.cosine()) + square(alpha.sine())/r, -3.0/2.0 );
                            if ( ! nearly_equal( current_PO, reference_PO ) )
                            {
                                std::cout << space_group.name() << " " << PO_directions[iPOD] << " " << reference_PO << " " << current_PO << " " << reflection_list.miller_indices( iReflection ) << " " << std::endl;
                            }

                        }
                    }
                }
            }
            iLine += 10;
        }
    MACRO_END_GAME

    try // .cell to .cif
    {
        if ( argc == 1 )
        {
            char a;
            std::cout << "Usage:" << std::endl;
            std::cout << std::endl;
            std::cout << "Cell2Cif.exe <CrystalStructure.cell>" << std::endl;
            std::cout << std::endl;
            std::cout << "Output: a cif file." << std::endl;
            std::cin >> a;
            return 0;
        }
        FileName file_name( argv[ 1 ] );
        if ( to_upper( file_name.extension() ) != "CELL" )
            throw std::runtime_error( "Please give the name of a .cell file that needs to be converted to a .cif file." );
        std::cout << "Now reading cell... " + file_name.full_name() << std::endl;
        CrystalStructure crystal_structure;
        read_cell( file_name, crystal_structure );
        crystal_structure.save_cif( replace_extension( file_name, "cif" ) );
    MACRO_END_GAME

    try // Calculate angles between two phenyl rings for a list of cifs.
    {
        MACRO_ONE_FILELISTNAME_AS_ARGUMENT
        TextFileWriter text_file_writer( FileName( file_list_file_name.directory(), "phenyl_angles", "txt" ) );
        std::vector< std::string > ring_1;
        std::vector< std::string > ring_2;
        ring_1.push_back( "C0" );
        ring_1.push_back( "C1" );
        ring_1.push_back( "C2" );
        ring_1.push_back( "C3" );
        ring_1.push_back( "C4" );
        ring_1.push_back( "C5" );
        ring_2.push_back( "C6" );
        ring_2.push_back( "C7" );
        ring_2.push_back( "C8" );
        ring_2.push_back( "C9" );
        ring_2.push_back( "C10" );
        ring_2.push_back( "C11" );
        // Note that for Z'=1 structures, we have e.g. "C1_0", for Z'=2 structures, we have e.g. "C1_0" and "C1_1".
        std::string header_str( "Rank " );
        // Write header
        text_file_writer.write_line( header_str );
        for ( size_t i( 0 ); i != file_list.size(); ++i )
//        for ( size_t i( 0 ); i != 10; ++i )
        {
            CrystalStructure crystal_structure;
            std::cout << "Now reading cif... " + file_list.value( i ).full_name() << std::endl;
            read_cif( file_list.value( i ), crystal_structure );
            CrystalLattice crystal_lattice = crystal_structure.crystal_lattice();
            bool Zprime_is_two( false );
            try
            {
                /*size_t index =*/ crystal_structure.atom( ring_1[0] + "_1" );
                Zprime_is_two = true;
            }
            catch ( std::exception & e ) {}
            std::string output_Z1_string( size_t2string( i+1 ) + " " );
            std::string output_Z2_string( size_t2string( i+1 ) + " " );
            { // Z' = 1
                std::vector< Vector3D > points_1;
                for ( size_t i( 0 ); i != ring_1.size(); ++i )
                {
                    size_t i1 = crystal_structure.atom( ring_1[i] + "_0" );
                    Atom atom = crystal_structure.atom( i1 );
                    points_1.push_back( crystal_lattice.fractional_to_orthogonal( atom.position() ) );
                }
                Plane plane_1( points_1 );
                std::vector< Vector3D > points_2;
                for ( size_t i( 0 ); i != ring_2.size(); ++i )
                {
                    size_t i1 = crystal_structure.atom( ring_2[i] + "_0" );
                    Atom atom = crystal_structure.atom( i1 );
                    points_2.push_back( crystal_lattice.fractional_to_orthogonal( atom.position() ) );
                }
                Plane plane_2( points_2 );
                Angle phi = angle( plane_1, plane_2 );
                output_Z1_string += double2string( phi.value_in_degrees() ) + " ";
                text_file_writer.write_line( output_Z1_string );
            }
            if ( Zprime_is_two )
            {
                std::vector< Vector3D > points_1;
                for ( size_t i( 0 ); i != ring_1.size(); ++i )
                {
                    size_t i1 = crystal_structure.atom( ring_1[i] + "_1" );
                    Atom atom = crystal_structure.atom( i1 );
                    points_1.push_back( crystal_lattice.fractional_to_orthogonal( atom.position() ) );
                }
                Plane plane_1( points_1 );
                std::vector< Vector3D > points_2;
                for ( size_t i( 0 ); i != ring_2.size(); ++i )
                {
                    size_t i1 = crystal_structure.atom( ring_2[i] + "_1" );
                    Atom atom = crystal_structure.atom( i1 );
                    points_2.push_back( crystal_lattice.fractional_to_orthogonal( atom.position() ) );
                }
                Plane plane_2( points_2 );
                Angle phi = angle( plane_1, plane_2 );
                output_Z2_string += double2string( phi.value_in_degrees() ) + " ";
                text_file_writer.write_line( output_Z2_string );
            }
        }
    MACRO_END_GAME

    try // Write out line number, R(, F( of TMFF file.
    {
        if ( argc != 2 )
        {
            std::cout << "Please give the name of a TMFF .ff file" << std::endl;
            throw std::runtime_error( "Please supply a <directory> and a <base name> as command-line arguments." );
        }
        FileName input_file_file_name( argv[ 1 ] );
        TextFileReader_2 input_file( input_file_file_name );
        TextFileWriter output_file( replace_extension( FileName( argv[ 1 ] ), "txt" )  );
        std::string line;
        for ( size_t iLine( 0 ); iLine != input_file.size(); ++iLine )
        {
            line = input_file.line( iLine );
            size_t iPos = 0;
            // Find all occurrences of R(
            while ( iPos != std::string::npos )
            {
                iPos = line.find( "R(", iPos );
                if ( iPos != std::string::npos )
                {
                    size_t iPos2 = line.find( ",", iPos );
                    output_file.write_line( size_t2string( iLine+1, 6, ' ' ) + " " + line.substr( iPos+2, iPos2-iPos-2 ) );
                    iPos = iPos2;
                }
            }
            // Find all occurrences of F(
            iPos = 0;
            while ( iPos != std::string::npos )
            {
                iPos = line.find( "F(", iPos );
                if ( iPos != std::string::npos )
                {
                    size_t iPos2 = line.find( ",", iPos );
                    output_file.write_line( size_t2string( iLine+1, 6, ' ' ) + " " + line.substr( iPos+2, iPos2-iPos-2 ) );
                    iPos = iPos2;
                }
            }
        }
    MACRO_END_GAME

    try // Generate R input file for Pawley or Loopstra-Rietveld plot.
        {
        if ( argc != 3 )
        {
            std::cout << "The following files are required:" << std::endl;
            std::cout << "<directory>\\<base name>.xye" << std::endl;
            std::cout << "<directory>\\<base name>_profile.txt" << std::endl;
            std::cout << "<directory>\\<base name>_tickmarks.txt" << std::endl;
            throw std::runtime_error( "Please supply a <directory> and a <base name> as command-line arguments." );
        }
        std::string directory( append_backslash( argv[ 1 ] ) );
        std::string base_name( argv[ 2 ] );
        GeneratePowderCIF generate_powder_cif( directory, base_name );
    //enum ZoomPolicy { ALWAYS_ZOOM, NEVER_ZOOM, ZOOM_OVER_40 };
        generate_powder_cif.generate_R_input_file( GeneratePowderCIF::NEVER_ZOOM );
    MACRO_END_GAME

    try // test
    {
        FileName input_file_name( argv[ 1 ] );
        PowderPattern powder_pattern;
        powder_pattern.read_xye( input_file_name );
        std::cout << "Cumulative intensity = " << powder_pattern.cumulative_intensity() << std::endl;
    MACRO_END_GAME

    try // Convert powder pattern in ASCII .raw format to .xye.
    {
        if ( argc != 2 )
            throw std::runtime_error( "Please give the name of a .raw file." );
        FileName input_file_name( argv[ 1 ] );
        PowderPattern powder_pattern;
        powder_pattern.read_raw( input_file_name );
        powder_pattern.save_xye( replace_extension( input_file_name, "xye" ), true );
    MACRO_END_GAME

    // Rotate group.
    try
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        // The group is attached to atom 1
        std::string atom_1_label( "C18" );
        std::string atom_2_label( "C13" );
        Vector3D C1 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_1_label ) ).position() );
        Vector3D C2 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_2_label ) ).position() );
        NormalisedVector3D n = normalised_vector( C1 - C2 );
        std::vector< std::string > labels_to_be_rotated;
        labels_to_be_rotated.push_back( "H21" );
        labels_to_be_rotated.push_back( "H22" );
        labels_to_be_rotated.push_back( "C21" );
        labels_to_be_rotated.push_back( "H24" );
        labels_to_be_rotated.push_back( "H25" );
        labels_to_be_rotated.push_back( "H26" );
        Angle additional_angle = Angle::from_degrees( 120.0 );
        for ( size_t i( 0 ); i != labels_to_be_rotated.size(); ++i )
        {
            Atom new_atom( crystal_structure.atom( crystal_structure.find_label( labels_to_be_rotated[i] ) ) );
            Vector3D tVector = crystal_structure.crystal_lattice().fractional_to_orthogonal( new_atom.position() );
            tVector = rotate_point_about_axis( tVector, C1, n, additional_angle );
            new_atom.set_position( crystal_structure.crystal_lattice().orthogonal_to_fractional( tVector ) );
            crystal_structure.set_atom( crystal_structure.find_label( labels_to_be_rotated[i] ), new_atom );
        }
        crystal_structure.save_cif( append_to_file_name( input_file_name, "_group_rotated" ) );
    MACRO_END_GAME

    // Add powder diffraction patterns and add noise and a background.
    try
    {
        PowderPattern powder_pattern_1;
        powder_pattern_1.read_xye( FileName( "W:\\c71\\p0024\\Form_I.xye" ) );
        bool add_second_pattern( false );
        PowderPattern powder_pattern_2;
        if ( add_second_pattern )
        {
            powder_pattern_2.read_xye( FileName( "\\\\Mac\\Home\\Documents\\Data_mac\\ContractResearch\\AMS\\Loratadine\\GWO30a.xye" ) );
        }
        PowderPattern background_1;
        background_1.read_xye( FileName( "W:\\c71\\p0024\\\\Form_I_BKGR.xye" ) );
        background_1.normalise_highest_peak( 1000.0 );
        PowderPattern background_2;
        if ( add_second_pattern )
        {
            background_2.read_xye( FileName( "\\\\Mac\\Home\\Documents\\Data_mac\\ContractResearch\\AMS\\Loratadine\\GWO30a_BKGR.xye" ) );
            background_2.normalise_highest_peak( 1000.0 );
        }
        powder_pattern_1.add_constant_background( 20.0 );
//        powder_pattern_1 += background_1;
        if ( add_second_pattern )
        {
            powder_pattern_1 += powder_pattern_2;
            powder_pattern_1 += background_2;
        }
        powder_pattern_1.add_Poisson_noise();
        powder_pattern_1.save_xye( FileName( "W:\\c71\\p0024\\Form_I_sum.xye" ), true );
    MACRO_END_GAME

    try // find_match( CrystalStructure, CrystalStructure ) crystal structure 2 (rhs) is the one that gets changed, so 1 is the target.
    {
        if ( argc < 3 )
            throw std::runtime_error( "Please give the name of two .cif files." );
        FileName file_name_1;
        FileName file_name_2;
        file_name_1 = FileName( argv[ 1 ] );
        file_name_2 = FileName( argv[ 2 ] );
        CrystalStructure crystal_structure_1;
        std::cout << "Now reading cif... " + file_name_1.full_name() << std::endl;
        read_cif( file_name_1, crystal_structure_1 );
        CrystalStructure crystal_structure_2;
        std::cout << "Now reading cif... " + file_name_2.full_name() << std::endl;
        read_cif( file_name_2, crystal_structure_2 );
        std::vector< int > integer_shifts;
        SymmetryOperator symmetry_operator = find_match( crystal_structure_1, crystal_structure_2, 4, integer_shifts, true, true );
        Matrix3D rotation = symmetry_operator.rotation();
        Vector3D shift = symmetry_operator.translation();
        shift.set_x( shift.x() + integer_shifts[0] );
        shift.set_y( shift.y() + integer_shifts[1] );
        shift.set_z( shift.z() + integer_shifts[2] );
        for ( size_t i( 0 ); i != crystal_structure_2.natoms(); ++i )
        {
            Atom new_atom( crystal_structure_2.atom( i ) );
            new_atom.set_position( ( rotation * crystal_structure_2.atom( i ).position() ) + shift );
            if ( new_atom.ADPs_type() == Atom::ANISOTROPIC )
                new_atom.set_anisotropic_displacement_parameters( rotate_adps( new_atom.anisotropic_displacement_parameters(), rotation, crystal_structure_2.crystal_lattice() ) );
            crystal_structure_2.set_atom( i, new_atom );
        }
        crystal_structure_2.save_cif( append_to_file_name( file_name_2, "_mtrans" ) );
    MACRO_END_GAME

    try // Find structure in FileList.txt using Rene de Gelder's normalised weighted cross correlations.
    {
        if ( argc != 3 )
            throw std::runtime_error( "Please give the name of a .cif file and a FileList.txt file." );
        FileName target_file_name( argv[ 1 ] );
        FileName file_list_file_name( argv[ 2 ] );
        if ( to_lower( file_list_file_name.extension() ) == "cif" )
            std::swap( target_file_name, file_list_file_name );
        CrystalStructure target_crystal_structure;
        std::cout << "Now reading cif... " + target_file_name.full_name() << std::endl;
        read_cif( target_file_name, target_crystal_structure );
        target_crystal_structure.apply_space_group_symmetry();
        Angle two_theta_start( 3.0, Angle::DEGREES );
        Angle two_theta_end(  35.0, Angle::DEGREES );
        Angle two_theta_step( 0.01, Angle::DEGREES );
        double FWHM( 0.1 );
        PowderPattern target_powder_pattern;
        {
            PowderPatternCalculator powder_pattern_calculator( target_crystal_structure );
            powder_pattern_calculator.set_wavelength( 1.54056 );
            powder_pattern_calculator.set_two_theta_start( two_theta_start );
            powder_pattern_calculator.set_two_theta_end( two_theta_end );
            powder_pattern_calculator.set_two_theta_step( two_theta_step );
            powder_pattern_calculator.set_FWHM( FWHM );
            powder_pattern_calculator.calculate( target_powder_pattern );
        }
        FileList file_list( file_list_file_name );
        if ( file_list.empty() )
            throw std::runtime_error( std::string( "No files in file list " ) + file_list_file_name.full_name() );
        double highest_correlation( 0.0 );
        size_t highest_correlation_index( 0 );
        TextFileWriter text_file_writer( FileName( "C:\\Data_Win\\matches.txt" ) );
        std::vector< std::string > water_labels;
        water_labels.push_back( "O0_1" );
        water_labels.push_back( "H0_1" );
        water_labels.push_back( "H1_1" );
        water_labels.push_back( "O0_2" );
        water_labels.push_back( "H0_2" );
        water_labels.push_back( "H1_2" );
        for ( size_t i( 0 ); i != file_list.size(); ++i )
        {
            CrystalStructure crystal_structure;
 //           std::cout << "Now reading cif... " + file_list.value( i ).full_name() << std::endl;
            read_cif( file_list.value( i ), crystal_structure );

            // ### REMOVE WATER ###
            for ( size_t k( 0 ); k != water_labels.size(); ++k )
            {
                size_t iAtom = crystal_structure.find_label( water_labels[k] );
                if ( iAtom != crystal_structure.natoms() )
                {
                    Atom new_atom = crystal_structure.atom( iAtom );
                    new_atom.set_occupancy( 0.0 );
                    crystal_structure.set_atom( iAtom, new_atom );
                }
            }
            crystal_structure.apply_space_group_symmetry();
//            std::cout << "Now calculating powder pattern... " + size_t2string( i, 4, '0' ) << std::endl;
            PowderPatternCalculator powder_pattern_calculator( crystal_structure );
            powder_pattern_calculator.set_wavelength( 1.54056 );
            powder_pattern_calculator.set_two_theta_start( two_theta_start );
            powder_pattern_calculator.set_two_theta_end( two_theta_end );
            powder_pattern_calculator.set_two_theta_step( two_theta_step );
            powder_pattern_calculator.set_FWHM( FWHM );
            PowderPattern powder_pattern;
            powder_pattern_calculator.calculate( powder_pattern );
            double correlation = normalised_weighted_cross_correlation( target_powder_pattern, powder_pattern, Angle( 3.0, Angle::DEGREES ) );
            if ( correlation > 0.95 )
            {
                text_file_writer.write_line( double2string( correlation ) + " " + size_t2string( i+1 ) );
            }
            if ( correlation > highest_correlation )
            {
                highest_correlation = correlation;
                highest_correlation_index = i;
                std::cout << "highest_correlation = " << highest_correlation << std::endl;
                std::cout << "highest_correlation_index = " << highest_correlation_index+1 << std::endl;
            }
        }
    MACRO_END_GAME

    try // Find unit-cell with all angles greater or smaller than 90 degrees (necessary to e.g. convert P1 to standard setting).
    // Also allows unit-cell axes to be ordered by length
    {
        if ( argc < 2 )
            throw std::runtime_error( "Please give the name of a .cif file." );
        FileName input_file_name( argv[ 1 ] );
        CrystalStructure crystal_structure;
        std::cout << "Now reading cif... " + input_file_name.full_name() << std::endl;
        read_cif( input_file_name, crystal_structure );
        CrystalLattice old_crystal_lattice = crystal_structure.crystal_lattice();
        Matrix3D best_transformation_matrix;
        double best_difference( 1080.0 );
        bool order_lengths_ascending( true );
        bool order_lengths_descending( false );
        if ( order_lengths_ascending && order_lengths_descending )
            throw std::runtime_error( "order_lengths_ascending and order_lengths_descending are both true." );
        bool angles_greater( true );
        int limit = 3;
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
                    CrystalLattice new_crystal_lattice( old_crystal_lattice );
                    // Transform
                    Matrix3D transformation_matrix( i1, i2, i3, j1, j2, j3, k1, k2, k3 );
                    if ( ! nearly_equal( transformation_matrix.determinant(), 1.0 ) )
                        continue;
                    new_crystal_lattice.transform( transformation_matrix );
                    Angle angle_sum;
                    if ( order_lengths_ascending )
                    {
                        if ( new_crystal_lattice.a() > new_crystal_lattice.b() )
                            continue;
                        if ( new_crystal_lattice.b() > new_crystal_lattice.c() )
                            continue;
                        if ( new_crystal_lattice.a() > new_crystal_lattice.c() )
                            continue;
                    }
                    if ( order_lengths_descending )
                    {
                        if ( new_crystal_lattice.a() < new_crystal_lattice.b() )
                            continue;
                        if ( new_crystal_lattice.b() < new_crystal_lattice.c() )
                            continue;
                        if ( new_crystal_lattice.a() < new_crystal_lattice.c() )
                            continue;
                    }
                    if ( angles_greater )
                    {

                        if ( new_crystal_lattice.alpha() < Angle::angle_90_degrees() )
                            continue;
                        if ( new_crystal_lattice.beta()  < Angle::angle_90_degrees() )
                            continue;
                        if ( new_crystal_lattice.gamma() < Angle::angle_90_degrees() )
                            continue;
                    }
                    else // smaller
                    {
                        if ( new_crystal_lattice.alpha() > Angle::angle_90_degrees() )
                            continue;
                        if ( new_crystal_lattice.beta()  > Angle::angle_90_degrees() )
                            continue;
                        if ( new_crystal_lattice.gamma() > Angle::angle_90_degrees() )
                            continue;
                    }
                    angle_sum = new_crystal_lattice.alpha() + new_crystal_lattice.beta() + new_crystal_lattice.gamma();
                    if ( std::abs( angle_sum.value_in_degrees() - 270.0 ) < best_difference )
                    {
                        best_difference = std::abs( angle_sum.value_in_degrees() - 270.0 );
                        best_transformation_matrix = transformation_matrix;
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
        best_transformation_matrix.show();
        CrystalLattice new_crystal_lattice( old_crystal_lattice );
        new_crystal_lattice.transform( best_transformation_matrix );
        new_crystal_lattice.print();
        std::cout << "Best difference = " << best_difference << std::endl;
        std::cout << std::endl;
    MACRO_END_GAME

    try // From C-centred to primitive, then find unit-cell angles close to 90 degrees.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        crystal_structure.transform( Matrix3D(  0.5,  0.5,  0.0,
                                               -0.5,  0.5,  0.0,
                                                0.0,  0.0,  1.0 ) );
        SpaceGroup space_group = crystal_structure.space_group();
        // In Mercury, if the space-group name and the set of symmetry operators do not match up,
        // the space-group name takes precedence, so we have to erase it to ensure that the
        // symmetry operators are used instead.
        space_group.set_name( "" );
        space_group.remove_duplicate_symmetry_operators();
        crystal_structure.set_space_group( space_group );
        CrystalLattice old_crystal_lattice = crystal_structure.crystal_lattice();
        Matrix3D identity_matrix;
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
                    // Make a copy
                    CrystalLattice new_lattice( old_crystal_lattice );
                    // Transform
                    Matrix3D transformation_matrix( i1, i2, i3, j1, j2, j3, k1, k2, k3 );
                    if ( ! nearly_equal( transformation_matrix.determinant(), 1.0 ) )
                        continue;
                    new_lattice.transform( transformation_matrix );
                    Angle tolerance = Angle::from_degrees( 70.0 );
                    if ( ( new_lattice.alpha() < tolerance ) || ( new_lattice.alpha() > ( Angle::angle_180_degrees() - tolerance ) ) )
                        continue;
                    if ( ( new_lattice.beta()  < tolerance ) || ( new_lattice.beta()  > ( Angle::angle_180_degrees() - tolerance ) ) )
                        continue;
                    if ( ( new_lattice.gamma() < tolerance ) || ( new_lattice.gamma() > ( Angle::angle_180_degrees() - tolerance ) ) )
                        continue;
                    {
                        transformation_matrix.show();
                        new_lattice.print();
                        Angle biggest_difference = absolute( new_lattice.alpha() - Angle::angle_90_degrees() );
                        biggest_difference = std::max( biggest_difference, absolute( new_lattice.beta() - Angle::angle_90_degrees() ) );
                        biggest_difference = std::max( biggest_difference, absolute( new_lattice.gamma() - Angle::angle_90_degrees() ) );
                        std::cout << "  " << (transformation_matrix-identity_matrix).sum_of_absolute_elements() << "   biggest diff = " << biggest_difference << std::endl;
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

    // From C-centred to primitive. Don't forget to make the unit-cell angles closer to 90 degrees.
    try
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        crystal_structure.transform( Matrix3D(  0.5,  0.5,  0.0,
                                               -0.5,  0.5,  0.0,
                                                0.0,  0.0,  1.0 ) );
        SpaceGroup space_group = crystal_structure.space_group();
        // In Mercury, if the space-group name and the set of symmetry operators do not match up,
        // the space-group name takes precedence, so we have to erase it to ensure that the
        // symmetry operators are used instead.
        space_group.set_name( "" );
        space_group.remove_duplicate_symmetry_operators();
        crystal_structure.set_space_group( space_group );
        Matrix3D rotation(  1.0,  0.0,  0.0,
                            0.0,  1.0,  0.0,
                            0.0,  0.0,  1.0 );
        Vector3D shift( 0.0, 0.0, 0.0 );
        for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
        {
            Atom new_atom( crystal_structure.atom( i ) );
            new_atom.set_position( ( rotation * crystal_structure.atom( i ).position() ) + shift );
            if ( new_atom.ADPs_type() == Atom::ANISOTROPIC )
            {
                SymmetricMatrix3D U_cif = new_atom.anisotropic_displacement_parameters().U_cif( crystal_structure.crystal_lattice() );
                U_cif = Matrix3D2SymmetricMatrix3D( rotation * U_cif );
                new_atom.set_anisotropic_displacement_parameters( U_cif_2_U_cart( U_cif, crystal_structure.crystal_lattice() ) );
            }
            crystal_structure.set_atom( i, new_atom );
        }
        crystal_structure.save_cif( append_to_file_name( input_file_name, "_C2P" ) );
    MACRO_END_GAME

    // Calculate molecular volume and packing coefficient.
    try
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        crystal_structure.apply_space_group_symmetry();
        double total_void_volume = void_volume( crystal_structure );
        std::cout << "Total void volume = " << double2string( total_void_volume ) << std::endl;
        std::cout << "Unit-cell volume = " <<  double2string( crystal_structure.crystal_lattice().volume() ) << std::endl;
        std::cout << "Molecular volume = " << double2string( ( crystal_structure.crystal_lattice().volume() - total_void_volume ) / crystal_structure.space_group().nsymmetry_operators() ) << std::endl;
        std::cout << "Packing coefficient = " << double2string( ( crystal_structure.crystal_lattice().volume() - total_void_volume ) / crystal_structure.crystal_lattice().volume() ) << std::endl;
    MACRO_END_GAME

    try // Background subtraction by Brueckner.
    {
        MACRO_ONE_XYEFILENAME_AS_ARGUMENT
        std::cout << powder_pattern.average_two_theta_step() << std::endl;
        PowderPattern background = calculate_Brueckner_background( powder_pattern,
                                                                   50, // niterations
                                                                   round_to_int( 50.0 * ( Angle::from_degrees( 0.015 ) / powder_pattern.average_two_theta_step() ) ), // window
                                                                   true, // apply_smoothing
                                                                   5 ); // smoothing_window
        background.save_xye( append_to_file_name( input_file_name, "_BKGR" ), true );
        powder_pattern -= background;
        powder_pattern.save_xye( append_to_file_name( input_file_name, "_NOBKGR" ), true );
    MACRO_END_GAME

    try // Generate a powder pattern based on a list of d-spacings and intensities.
    {
        // We need a dummy crystal stucture.
        CrystalStructure crystal_structure;
        crystal_structure.apply_space_group_symmetry();
        PowderPatternCalculator powder_pattern_calculator( crystal_structure );
        powder_pattern_calculator.set_two_theta_start( Angle::from_degrees( 0.0 ) );
        powder_pattern_calculator.set_two_theta_end( Angle::from_degrees( 47.0 ) );
        powder_pattern_calculator.set_two_theta_step( Angle::from_degrees( 0.02 ) );
        PowderPattern powder_pattern;
        double FWHM = 0.5;
        double lambda = 1.54056;
        std::vector< double > peak_positions;
        std::vector< double > F2_values;
        std::vector< double > FWHM_values;
        peak_positions.push_back( 2.0 * arcsine( lambda / ( 2.0 * 9.23 ) ).value_in_degrees() ); F2_values.push_back(   500.0 ); FWHM_values.push_back( FWHM );
        peak_positions.push_back( 2.0 * arcsine( lambda / ( 2.0 * 7.51 ) ).value_in_degrees() ); F2_values.push_back(   500.0 ); FWHM_values.push_back( FWHM );
        peak_positions.push_back( 2.0 * arcsine( lambda / ( 2.0 * 6.31 ) ).value_in_degrees() ); F2_values.push_back(  3000.0 ); FWHM_values.push_back( FWHM );
        peak_positions.push_back( 2.0 * arcsine( lambda / ( 2.0 * 5.80 ) ).value_in_degrees() ); F2_values.push_back(  2000.0 ); FWHM_values.push_back( FWHM );
        peak_positions.push_back( 2.0 * arcsine( lambda / ( 2.0 * 5.44 ) ).value_in_degrees() ); F2_values.push_back(  1000.0 ); FWHM_values.push_back( FWHM );
        peak_positions.push_back( 2.0 * arcsine( lambda / ( 2.0 * 5.22 ) ).value_in_degrees() ); F2_values.push_back(  4000.0 ); FWHM_values.push_back( FWHM );
        peak_positions.push_back( 2.0 * arcsine( lambda / ( 2.0 * 4.91 ) ).value_in_degrees() ); F2_values.push_back(   500.0 ); FWHM_values.push_back( FWHM );
        peak_positions.push_back( 2.0 * arcsine( lambda / ( 2.0 * 4.55 ) ).value_in_degrees() ); F2_values.push_back( 10000.0 ); FWHM_values.push_back( FWHM );
        peak_positions.push_back( 2.0 * arcsine( lambda / ( 2.0 * 4.31 ) ).value_in_degrees() ); F2_values.push_back(  4000.0 ); FWHM_values.push_back( FWHM );
        peak_positions.push_back( 2.0 * arcsine( lambda / ( 2.0 * 4.13 ) ).value_in_degrees() ); F2_values.push_back(  3000.0 ); FWHM_values.push_back( FWHM );
        peak_positions.push_back( 2.0 * arcsine( lambda / ( 2.0 * 3.99 ) ).value_in_degrees() ); F2_values.push_back(  9000.0 ); FWHM_values.push_back( FWHM );
        peak_positions.push_back( 2.0 * arcsine( lambda / ( 2.0 * 3.78 ) ).value_in_degrees() ); F2_values.push_back(  2000.0 ); FWHM_values.push_back( FWHM );
        peak_positions.push_back( 2.0 * arcsine( lambda / ( 2.0 * 3.63 ) ).value_in_degrees() ); F2_values.push_back( 10000.0 ); FWHM_values.push_back( FWHM );
        peak_positions.push_back( 2.0 * arcsine( lambda / ( 2.0 * 3.46 ) ).value_in_degrees() ); F2_values.push_back(  1000.0 ); FWHM_values.push_back( FWHM );
        peak_positions.push_back( 2.0 * arcsine( lambda / ( 2.0 * 3.24 ) ).value_in_degrees() ); F2_values.push_back(  4000.0 ); FWHM_values.push_back( FWHM );
        peak_positions.push_back( 2.0 * arcsine( lambda / ( 2.0 * 3.17 ) ).value_in_degrees() ); F2_values.push_back(  5000.0 ); FWHM_values.push_back( FWHM );
        peak_positions.push_back( 2.0 * arcsine( lambda / ( 2.0 * 3.04 ) ).value_in_degrees() ); F2_values.push_back(  1000.0 ); FWHM_values.push_back( FWHM );
        peak_positions.push_back( 2.0 * arcsine( lambda / ( 2.0 * 2.93 ) ).value_in_degrees() ); F2_values.push_back(  4000.0 ); FWHM_values.push_back( FWHM );
        peak_positions.push_back( 2.0 * arcsine( lambda / ( 2.0 * 2.87 ) ).value_in_degrees() ); F2_values.push_back(  5000.0 ); FWHM_values.push_back( FWHM );
        peak_positions.push_back( 2.0 * arcsine( lambda / ( 2.0 * 2.78 ) ).value_in_degrees() ); F2_values.push_back(  4000.0 ); FWHM_values.push_back( FWHM );
        peak_positions.push_back( 2.0 * arcsine( lambda / ( 2.0 * 2.72 ) ).value_in_degrees() ); F2_values.push_back(  5000.0 ); FWHM_values.push_back( FWHM );
        peak_positions.push_back( 2.0 * arcsine( lambda / ( 2.0 * 2.69 ) ).value_in_degrees() ); F2_values.push_back(  1000.0 ); FWHM_values.push_back( FWHM );
        peak_positions.push_back( 2.0 * arcsine( lambda / ( 2.0 * 2.58 ) ).value_in_degrees() ); F2_values.push_back(  2000.0 ); FWHM_values.push_back( FWHM );
        peak_positions.push_back( 2.0 * arcsine( lambda / ( 2.0 * 2.46 ) ).value_in_degrees() ); F2_values.push_back(  8000.0 ); FWHM_values.push_back( FWHM );
        peak_positions.push_back( 2.0 * arcsine( lambda / ( 2.0 * 2.44 ) ).value_in_degrees() ); F2_values.push_back(  3000.0 ); FWHM_values.push_back( FWHM );
        peak_positions.push_back( 2.0 * arcsine( lambda / ( 2.0 * 2.36 ) ).value_in_degrees() ); F2_values.push_back(  1000.0 ); FWHM_values.push_back( FWHM );
        peak_positions.push_back( 2.0 * arcsine( lambda / ( 2.0 * 2.31 ) ).value_in_degrees() ); F2_values.push_back(  3000.0 ); FWHM_values.push_back( FWHM );
        peak_positions.push_back( 2.0 * arcsine( lambda / ( 2.0 * 2.25 ) ).value_in_degrees() ); F2_values.push_back(  3000.0 ); FWHM_values.push_back( FWHM );
        peak_positions.push_back( 2.0 * arcsine( lambda / ( 2.0 * 2.12 ) ).value_in_degrees() ); F2_values.push_back(  3000.0 ); FWHM_values.push_back( FWHM );
        peak_positions.push_back( 2.0 * arcsine( lambda / ( 2.0 * 2.06 ) ).value_in_degrees() ); F2_values.push_back(  3000.0 ); FWHM_values.push_back( FWHM );
        peak_positions.push_back( 2.0 * arcsine( lambda / ( 2.0 * 1.99 ) ).value_in_degrees() ); F2_values.push_back(  3000.0 ); FWHM_values.push_back( FWHM );
     //   powder_pattern_calculator.calculate_for_testing( peak_positions, F2_values, FWHM_values, powder_pattern );
        powder_pattern.save_xye( FileName( "/Volumes/Staff/jvds/AMS_ModelSystems/EthylenediamineTartrate/powder_pattern.xye" ), false );
    MACRO_END_GAME

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

    // Calculate powder diffraction pattern.
    try
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        crystal_structure.apply_space_group_symmetry();
        Angle two_theta_start( 3.0, Angle::DEGREES );
        Angle two_theta_end(  35.0, Angle::DEGREES );
        Angle two_theta_step( 0.01, Angle::DEGREES );
        double FWHM( 0.1 );
        PowderPattern powder_pattern;
        PowderPatternCalculator powder_pattern_calculator( crystal_structure );
        powder_pattern_calculator.set_wavelength( 1.54056 );
        powder_pattern_calculator.set_two_theta_start( two_theta_start );
        powder_pattern_calculator.set_two_theta_end( two_theta_end );
        powder_pattern_calculator.set_two_theta_step( two_theta_step );
        powder_pattern_calculator.set_FWHM( FWHM );

        // ### PREFERRED ORIENTATION ###
      //  powder_pattern_calculator.set_preferred_orientation( MillerIndices( 0, 0, 1 ), 1.0 );

        powder_pattern_calculator.calculate( powder_pattern );
        powder_pattern.save_xye( replace_extension( input_file_name, "xye" ), false );
    MACRO_END_GAME

    // Find voids.
    try
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        crystal_structure.apply_space_group_symmetry();
        double probe_radius = 1.75;
        double volume = find_voids( crystal_structure, probe_radius );
        std::cout << double2string( volume ) + " " + double2string( volume / crystal_structure.space_group().nsymmetry_operators() ) << std::endl;
    MACRO_END_GAME

    try // Recalculate ESDs XRPD pattern.
    {
        MACRO_ONE_XYEFILENAME_AS_ARGUMENT
        powder_pattern.recalculate_estimated_standard_deviations();
        std::cout << powder_pattern.average_two_theta_step() << std::endl;
        powder_pattern.save_xye( append_to_file_name( input_file_name, "_new_ESDs" ), true );
    MACRO_END_GAME

    // Add OH hydrogen atom.
    try
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        // The hydrogen atom is attached to atom_1
        std::string atom_1_label( "O1" );
        std::string atom_2_label( "C13" );
        Vector3D atom_1 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_1_label ) ).position() );
        Vector3D atom_2 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_2_label ) ).position() );
        size_t nhydrogen_atoms = 12;
        std::vector< Vector3D > hydrogen_atoms = add_hydrogen_atoms( atom_1, atom_2, nhydrogen_atoms );
        for ( size_t i( 0 ); i != nhydrogen_atoms; ++i )
            crystal_structure.add_atom( Atom( Element( "H" ), crystal_structure.crystal_lattice().orthogonal_to_fractional( hydrogen_atoms[i] ), "H" + size_t2string( crystal_structure.natoms() ) ) );
        crystal_structure.save_cif( append_to_file_name( input_file_name, "_H_added" ) );
    MACRO_END_GAME

    try // Fake powder pattern consisting of two peaks.
    {
        TextFileWriter text_file_writer( FileName( "/Volumes/Staff/jvds/AMS_ModelSystems/EthylenediamineTartrate/powder_pattern.txt" ) );
        // We need a dummy crystal stucture.
        CrystalStructure crystal_structure;
        crystal_structure.apply_space_group_symmetry();
        PowderPatternCalculator powder_pattern_calculator( crystal_structure );
        powder_pattern_calculator.set_two_theta_start( Angle::from_degrees( 0.0 ) );
        powder_pattern_calculator.set_two_theta_end( Angle::from_degrees( 50.0 ) );
        powder_pattern_calculator.set_two_theta_step( Angle::from_degrees( 0.01 ) );
        PowderPattern reference_powder_pattern;
        {
        std::vector< double > peak_positions;
        std::vector< double > F2_values;
        std::vector< double > FWHM_values;
        peak_positions.push_back( 10.0 );
        F2_values.push_back( 100.0 );
        FWHM_values.push_back( 0.1 );
        peak_positions.push_back( 20.0 );
        F2_values.push_back( 10.0 );
        FWHM_values.push_back( 0.1 );
    //    powder_pattern_calculator.calculate_for_testing( peak_positions, F2_values, FWHM_values, reference_powder_pattern );
        }
        for ( size_t i( 0 ); i != 1000; ++i )
        {
            std::vector< double > peak_positions;
            std::vector< double > F2_values;
            std::vector< double > FWHM_values;
            peak_positions.push_back( 10.0 + ( i * powder_pattern_calculator.two_theta_step().value_in_degrees() ) );
            F2_values.push_back( 100.0 );
            FWHM_values.push_back( 0.15 );
            peak_positions.push_back( 20.0 - ( i * powder_pattern_calculator.two_theta_step().value_in_degrees() ) );
            F2_values.push_back( 10.0 );
            FWHM_values.push_back( 0.15 );
            PowderPattern powder_pattern;
    //        powder_pattern_calculator.calculate_for_testing( peak_positions, F2_values, FWHM_values, powder_pattern );
            double nwcc = normalised_weighted_cross_correlation( reference_powder_pattern, powder_pattern, Angle::from_degrees( 3.0 ) );
            // The first pattern is supposed to be the experimental pattern, and its ESDs are used as "the" weights. The ESDs of the second pattern are ignored.
            double R_wp = Rwp( reference_powder_pattern, powder_pattern );
            text_file_writer.write_line( size_t2string( i ) + " " + double2string( nwcc ) + " " + double2string( 1.0 - ( R_wp / 16.0 ) ) );
        }
//        powder_pattern.save_xye( FileName( "C:\\Data_Win\\TwoPeaks.xye" ), false );
    MACRO_END_GAME

    try // From R-centred to primitive, then find unit-cell angles close to 90 degrees.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        Matrix3D R_to_P( 2.0/3.0, 1.0/3.0, 1.0/3.0,
                         1.0/3.0, 2.0/3.0, 2.0/3.0,
                             0.0,     0.0,     1.0 );
        crystal_structure.transform( R_to_P );
        SpaceGroup space_group = crystal_structure.space_group();
        // In Mercury, if the space-group name and the set of symmetry operators do not match up,
        // the space-group name takes precedence, so we have to erase it to ensure that the
        // symmetry operators are used instead.
        space_group.set_name( "" );
        space_group.remove_duplicate_symmetry_operators();
        crystal_structure.set_space_group( space_group );
        CrystalLattice old_crystal_lattice = crystal_structure.crystal_lattice();
        Matrix3D identity_matrix;
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
                    // Make a copy
                    CrystalLattice new_lattice( old_crystal_lattice );
                    // Transform
                    Matrix3D transformation_matrix( i1, i2, i3, j1, j2, j3, k1, k2, k3 );
                    if ( ! nearly_equal( transformation_matrix.determinant(), 1.0 ) )
                        continue;
                    new_lattice.transform( transformation_matrix );
                    Angle tolerance = Angle::from_degrees( 68.0 );
                    if ( ( new_lattice.alpha() < tolerance ) || ( new_lattice.alpha() > ( Angle::angle_180_degrees() - tolerance ) ) )
                        continue;
                    if ( ( new_lattice.beta()  < tolerance ) || ( new_lattice.beta()  > ( Angle::angle_180_degrees() - tolerance ) ) )
                        continue;
                    if ( ( new_lattice.gamma() < tolerance ) || ( new_lattice.gamma() > ( Angle::angle_180_degrees() - tolerance ) ) )
                        continue;
                    {
                        transformation_matrix.show();
                        ( transformation_matrix * R_to_P ).show();
                        new_lattice.print();
                        Angle biggest_difference = absolute( new_lattice.alpha() - Angle::angle_90_degrees() );
                        biggest_difference = std::max( biggest_difference, absolute( new_lattice.beta() - Angle::angle_90_degrees() ) );
                        biggest_difference = std::max( biggest_difference, absolute( new_lattice.gamma() - Angle::angle_90_degrees() ) );
                        std::cout << "  " << (transformation_matrix-identity_matrix).sum_of_absolute_elements() << "   biggest diff = " << biggest_difference << std::endl;
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

    try //
    {
        PowderPattern powder_pattern;
        TextFileReader text_file_reader( FileName( "\\\\Mac\\Home\\Downloads\\crystals-09-00384-s001\\fosfomycin_xrpd_input.txt" ) );
        text_file_reader.set_skip_empty_lines( false );
        std::vector< std::string > words;
        Angle two_theta_start = Angle::from_degrees( 0.5 );
        Angle two_theta_step = Angle::from_degrees( 0.00100 );
        size_t i( 0 );
        while ( text_file_reader.get_next_line( words ) )
        {
            if ( words.size() != 4 )
                std::cout << "Warning: words.size() != 4" << std::endl;
            powder_pattern.push_back( ( i * two_theta_step ) + two_theta_start, string2double( words[0] ) );
            ++i;
        }
        powder_pattern.save_xye( FileName( "\\\\Mac\\Home\\Downloads\\crystals-09-00384-s001\\fosfomycin.xye" ), true );
    MACRO_END_GAME

    try // Apply zero-point error correction to a powder pattern.
    {
        MACRO_ONE_XYEFILENAME_AS_ARGUMENT
        powder_pattern.correct_zero_point_error( Angle::from_degrees( 0.16 ) );
        powder_pattern.save_xye( append_to_file_name( input_file_name, "_zp" ), true );
    MACRO_END_GAME

    try // Find transformation.
    {
        if ( argc < 2 )
            throw std::runtime_error( "Please give the name of a .cif file." );
        FileName input_file_name_1( argv[ 1 ] );
        CrystalStructure target_crystal_structure;
        std::cout << "Now reading cif... " + input_file_name_1.full_name() << std::endl;
        read_cif( input_file_name_1, target_crystal_structure );
        target_crystal_structure.apply_space_group_symmetry();

        FileName input_file_name_2( argv[ 2 ] );
        CrystalStructure crystal_structure;
        std::cout << "Now reading cif... " + input_file_name_2.full_name() << std::endl;
        read_cif( input_file_name_2, crystal_structure );

        Angle two_theta_start( 5.0, Angle::DEGREES );
        Angle two_theta_end(  35.0, Angle::DEGREES );
        Angle two_theta_step( 0.01, Angle::DEGREES );
        double FWHM( 0.1 );
        PowderPattern target_powder_pattern;
        {
        PowderPatternCalculator powder_pattern_calculator( target_crystal_structure );
        powder_pattern_calculator.set_wavelength( 1.54056 );
        powder_pattern_calculator.set_two_theta_start( two_theta_start );
        powder_pattern_calculator.set_two_theta_end( two_theta_end );
        powder_pattern_calculator.set_two_theta_step( two_theta_step );
        powder_pattern_calculator.set_FWHM( FWHM );
        powder_pattern_calculator.calculate( target_powder_pattern );
        }

        double greatest_similarity = -1.0;
        size_t shift_steps = 8;
        std::vector< Vector3D > shifts;
        if ( ( shift_steps == 0 ) || ( shift_steps == 1 ) )
            shifts.push_back( Vector3D() );
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

        for ( size_t iShifts( 0 ); iShifts != shifts.size(); ++iShifts )
        {
            std::cout << " ." << std::endl;
            for ( size_t k( 0 ); k != crystal_structure.space_group().nsymmetry_operators(); ++k )
            {
                // Make a copy
                CrystalStructure new_crystal_structure( crystal_structure );
                // Transform
                for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
                {
                    Atom new_atom( crystal_structure.atom( i ) );
                    Vector3D current_position = crystal_structure.space_group().symmetry_operator( k ) * ( crystal_structure.atom( i ).position() + shifts[iShifts] );
                    new_atom.set_position( current_position );
                    new_crystal_structure.set_atom( i, new_atom );
                }
                new_crystal_structure.apply_space_group_symmetry();
                PowderPatternCalculator powder_pattern_calculator( new_crystal_structure );
                powder_pattern_calculator.set_wavelength( 1.54056 );
                powder_pattern_calculator.set_two_theta_start( two_theta_start );
                powder_pattern_calculator.set_two_theta_end( two_theta_end );
                powder_pattern_calculator.set_two_theta_step( two_theta_step );
                powder_pattern_calculator.set_FWHM( FWHM );
                PowderPattern new_powder_pattern;
                powder_pattern_calculator.calculate( new_powder_pattern );
                double similarity = normalised_weighted_cross_correlation( target_powder_pattern, new_powder_pattern );
                if ( similarity > ( greatest_similarity + 0.0001 ) )
                {
                    greatest_similarity = similarity;
//                        transformation_matrix.show();
//                        new_lattice.print();
                    std::cout << "greatest_similarity = " << greatest_similarity << std::endl;
//                        std::cout << std::endl;
                }
            }
        }
    MACRO_END_GAME

    // Add centre of symmetry, not at origin.
    try
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        std::vector< SymmetryOperator > old_symmetry_operators = crystal_structure.space_group().symmetry_operators();
        std::vector< SymmetryOperator > new_symmetry_operators( old_symmetry_operators );
        SymmetryOperator additional_symmetry_operator( "-x+1/4,-y+1/4,-z+1/4" );
        for ( size_t i( 0 ); i != old_symmetry_operators.size(); ++i )
            new_symmetry_operators.push_back( additional_symmetry_operator * old_symmetry_operators[i] );
        SpaceGroup new_space_group( new_symmetry_operators );
        crystal_structure.set_space_group( new_space_group );
        crystal_structure.save_cif( append_to_file_name( input_file_name, "_newSpGr" ) );
    MACRO_END_GAME

    try
    {
        double epsilon = 1.0E-5;
        {
            double x = -1.0 - epsilon;
            std::cout << "x " << x << " " << static_cast<int>( ( x ) ) << std::endl;
        }
        {
            double x = -1.0 + epsilon;
            std::cout << "x " << x << " " << static_cast<int>( ( x ) ) << std::endl;
        }
        {
            double x =  0.0 - epsilon;
            std::cout << "x " << x << " " << static_cast<int>( ( x ) ) << std::endl;
        }
        {
            double x =  0.0 + epsilon;
            std::cout << "x " << x << " " << static_cast<int>( ( x ) ) << std::endl;
        }
        {
            double x =  1.0 - epsilon;
            std::cout << "x " << x << " " << static_cast<int>( ( x ) ) << std::endl;
        }
        {
            double x =  1.0 + epsilon;
            std::cout << "x " << x << " " << static_cast<int>( ( x ) ) << std::endl;
        }
        double x( 3.0 );
        std::cout << std::fixed;
        std::cout.precision(50);
        std::cout << x << std::endl;
    MACRO_END_GAME

    // Print c.o.m.
    try
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        Vector3D com = crystal_structure.centre_of_mass();
        std::cout << "Centre of mass = " << std::endl;
        com.show();
    MACRO_END_GAME

    try // Split cif with multiple crystal structures
    {
        FileName file_name( "C:\\Data_Win\\Research\\Refereeing\\GF.cif" );
        TextFileReader text_file_reader_1( file_name );
        text_file_reader_1.set_skip_empty_lines( true );
        std::string line;
        bool more_left = text_file_reader_1.get_next_line( line );
        text_file_reader_1.push_back_last_line();
        while ( more_left && text_file_reader_1.get_next_line( line ) )
        {
            size_t iPos = to_lower( line ).find( "data_" );
            if ( iPos == std::string::npos )
            {
                std::cout << line << std::endl;
                throw std::runtime_error( "First line does not contain data_ keyword." );
            }
            std::string identifier = strip( line.substr( iPos + 5 ) );
            FileName output_file_name;
            if ( identifier.empty() )
                output_file_name = generate_unique_file_name( file_name );
            else
                output_file_name = FileName( file_name.directory(), identifier, "cif" );
            TextFileWriter text_file_writer( output_file_name );
            text_file_writer.write_line( line ); // The line containing "data_"
            more_left = text_file_reader_1.get_next_line( line );
            while ( more_left && ( to_lower( line ).find( "data_" ) == std::string::npos ) )
            {
                text_file_writer.write_line( line );
                more_left = text_file_reader_1.get_next_line( line );
            }
            text_file_reader_1.push_back_last_line();
        }
    MACRO_END_GAME

    try // Generate FileList.txt .
    {
        if ( argc != 2 )
            throw std::runtime_error( "Please give the name of a directory." );
        TextFileWriter text_file_writer( FileName( argv[ 1 ], "FileList", "txt" ) );
        for ( size_t i( 0 ); i != 265; ++i )
            text_file_writer.write_line( "structure_" + size_t2string( i+1, 6, '0') + ".cif" );
    MACRO_END_GAME

    try // Calculate powder pattern.
    {
        if ( argc != 2 )
            throw std::runtime_error( "Please give the name of a .cif file." );
        FileName input_file_name( argv[ 1 ] );
        CrystalStructure crystal_structure;
        std::cout << "Now reading cif... " + input_file_name.full_name() << std::endl;
        read_cif( input_file_name, crystal_structure );
        crystal_structure.apply_space_group_symmetry();
        Angle two_theta_start( 5.0, Angle::DEGREES );
        Angle two_theta_end(  50.0, Angle::DEGREES );
        Angle two_theta_step( 0.02, Angle::DEGREES );
        double FWHM( 0.1 );
        std::cout << "Now calculating powder pattern... " << std::endl;
        PowderPatternCalculator powder_pattern_calculator( crystal_structure );
        powder_pattern_calculator.set_wavelength( 1.54056 );
        powder_pattern_calculator.set_two_theta_start( two_theta_start );
        powder_pattern_calculator.set_two_theta_end( two_theta_end );
        powder_pattern_calculator.set_two_theta_step( two_theta_step );
        powder_pattern_calculator.set_FWHM( FWHM );
//        powder_pattern_calculator.set_preferred_orientation( MillerIndices( 0, 0, 1 ), 0.75 );
        PowderPattern powder_pattern;
        powder_pattern_calculator.calculate( powder_pattern );
        powder_pattern.save_xye( FileName( input_file_name.directory(), input_file_name.file_name(), "xye" ), true );
    MACRO_END_GAME

    try // Show average 2theta step.
    {
        MACRO_ONE_XYEFILENAME_AS_ARGUMENT
        std::cout << "Average 2theta step: " << powder_pattern.average_two_theta_step() << std::endl;
    MACRO_END_GAME

    try // Convert powder pattern in .xrdml format to .xye.
    {
        if ( argc != 2 )
            throw std::runtime_error( "Please give the name of a .xrdml file." );
        FileName input_file_name( argv[ 1 ] );
        PowderPattern powder_pattern;
        powder_pattern.read_xrdml( input_file_name );
        std::cout << powder_pattern.average_two_theta_step() << std::endl;
        powder_pattern.save_xye( replace_extension( input_file_name, "xye" ), true );
    MACRO_END_GAME

    try // Repair XRPD pattern extracted from a .png
    {
        MACRO_ONE_XYEFILENAME_AS_ARGUMENT
        powder_pattern.sort_two_theta();
        powder_pattern.average_if_two_theta_equal();
        std::cout << powder_pattern.average_two_theta_step() << std::endl;
        powder_pattern.recalculate_estimated_standard_deviations();
        powder_pattern.save_xye( append_to_file_name( input_file_name, "_recal" ), false );
    MACRO_END_GAME

    try // Recalculate variable slit intensities to fixed slit intensities.
    {
        MACRO_ONE_XYEFILENAME_AS_ARGUMENT
        std::cout << powder_pattern.average_two_theta_step() << std::endl;
        powder_pattern.convert_to_fixed_slit();
        powder_pattern.save_xye( append_to_file_name( input_file_name, "_fixed_slits" ), true );
    MACRO_END_GAME

    try // Convert .gzmat to .zmatrix format (replace variables by their values).
    {
        if ( argc != 2 )
            throw std::runtime_error( "Please give the name of a .gzmat file." );
        FileName input_file_name( argv[ 1 ] );
        TextFileReader_2 text_file_reader( input_file_name );
        size_t iPos_Variables = text_file_reader.find_whole_word( "Variables:" );
        if ( iPos_Variables == 0 )
            throw std::runtime_error( "\"Variables:\" not found." );
        TextFileWriter text_file_writer( replace_extension( input_file_name, "zmatrix" ) );
        Splitter splitter( "=" );
        splitter.set_merge_delimiters( false );
        std::string line;
        std::vector< std::string > words;
//H  38  r75  3  a75  2  d75
//H  38  r76  3  a76  2  d76
//H  38  r77  3  a77  2  d77
//Variables:
//r2= 11.5889
//r3= 5.9129
//a3=  61.24
        for ( size_t iLine( 0 ); iLine != iPos_Variables; ++iLine )
        {
            words = split( text_file_reader.line( iLine ) );
            if ( words.size() == 0 )
            {
                std::cout << "WARNING: unexpected empty line." << std::endl;
                continue;
            }
            if ( is_even( words.size() ) )
                throw std::runtime_error( "Cannot interpret line with even number of words." );
            std::string new_line;
            new_line += words[0];
            if ( words.size() >= 3 )
                new_line += " " + words[1] + " " + extract_variable_value( text_file_reader.line( text_file_reader.find( words[2], iPos_Variables ) ) );
            if ( words.size() >= 5 )
                new_line += " " + words[3] + " " + extract_variable_value( text_file_reader.line( text_file_reader.find( words[4], iPos_Variables ) ) );
            if ( words.size() >= 7 )
                new_line += " " + words[5] + " " + extract_variable_value( text_file_reader.line( text_file_reader.find( words[6], iPos_Variables ) ) );
            if ( words.size() >= 9 )
                throw std::runtime_error( "Cannot interpret line with >7 words." );
            text_file_writer.write_line( new_line );
        }
    MACRO_END_GAME

    // Add two hydrogen atoms to an sp3 atom.
    try
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT

        {
        std::string neighbour_1_atom_label( "C20" );
        std::string origin_atom_label( "C21" ); // This is the atom to which the hydrogen atoms are added
        std::string neighbour_2_atom_label( "N4" );
        Vector3D neighbour_1_atom_frac = crystal_structure.atom( crystal_structure.find_label( neighbour_1_atom_label ) ).position();
        Vector3D origin_atom_frac = crystal_structure.atom( crystal_structure.find_label( origin_atom_label ) ).position();
        Vector3D neighbour_2_atom_frac = crystal_structure.atom( crystal_structure.find_label( neighbour_2_atom_label ) ).position();
        // Find shortest disance between the two taking periodicity and space-group symmetry into account.
        double distance;
        Vector3D r_1mino_frac;
        crystal_structure.shortest_distance( origin_atom_frac, neighbour_1_atom_frac, distance, r_1mino_frac );
        if ( distance > 3.0 )
            std::cout << "Warning: distance 1 > 3.0 A" << std::endl;
        Vector3D r_2mino_frac;
        crystal_structure.shortest_distance( origin_atom_frac, neighbour_2_atom_frac, distance, r_2mino_frac );
        if ( distance > 3.0 )
            std::cout << "Warning: distance 2 > 3.0 A" << std::endl;
        // We need to be able to set the length of the X-H bond, so we must work in Cartesian coordinates.
        Vector3D neighbour_1_atom_cart = crystal_structure.crystal_lattice().fractional_to_orthogonal( neighbour_1_atom_frac );
        Vector3D origin_atom_cart = crystal_structure.crystal_lattice().fractional_to_orthogonal( origin_atom_frac );
        Vector3D neighbour_2_atom_cart = crystal_structure.crystal_lattice().fractional_to_orthogonal( neighbour_2_atom_frac );
        Vector3D r_1mino_cart = crystal_structure.crystal_lattice().fractional_to_orthogonal( r_1mino_frac );
        NormalisedVector3D n_1mino_cart = normalised_vector( r_1mino_cart );
        Vector3D r_2mino_cart = crystal_structure.crystal_lattice().fractional_to_orthogonal( r_2mino_frac );
        NormalisedVector3D n_2mino_cart = normalised_vector( r_2mino_cart );
        NormalisedVector3D x = normalised_vector( -0.5*( n_1mino_cart + n_2mino_cart ) );
        Plane plane( neighbour_1_atom_cart, origin_atom_cart, neighbour_2_atom_cart );
        NormalisedVector3D y = plane.normal();
        double target_bond_length( 1.0 );
        if ( crystal_structure.atom( crystal_structure.find_label( origin_atom_label ) ).element() == Element( "C" ) )
            target_bond_length = 1.089;
        else if ( crystal_structure.atom( crystal_structure.find_label( origin_atom_label ) ).element() == Element( "N" ) )
            target_bond_length = 1.015;
        else if ( crystal_structure.atom( crystal_structure.find_label( origin_atom_label ) ).element() == Element( "O" ) )
            target_bond_length = 0.993;
        Angle angle = Angle::from_degrees( 106.0 / 2.0 );
        double o = target_bond_length * angle.sine();
        double a = target_bond_length * angle.cosine();
        // The y-axis is perpendicular to the plane of the three non-hydrogen atoms, so one hydrogen is added above (+y) and one below (-y) the plane.
        Vector3D H_atom_1_cart = origin_atom_cart + o*y + a*x;
        Vector3D H_atom_1_frac = crystal_structure.crystal_lattice().orthogonal_to_fractional( H_atom_1_cart );
        crystal_structure.add_atom( Atom( Element( "H" ), H_atom_1_frac, "H" + size_t2string( crystal_structure.natoms() ) ) );
        Vector3D H_atom_2_cart = origin_atom_cart - o*y + a*x;
        Vector3D H_atom_2_frac = crystal_structure.crystal_lattice().orthogonal_to_fractional( H_atom_2_cart );
        crystal_structure.add_atom( Atom( Element( "H" ), H_atom_2_frac, "H" + size_t2string( crystal_structure.natoms() ) ) );
        }
        {
        std::string neighbour_1_atom_label( "C18" );
        std::string origin_atom_label( "C19" ); // This is the atom to which the hydrogen atoms are added
        std::string neighbour_2_atom_label( "C20" );
        Vector3D neighbour_1_atom_frac = crystal_structure.atom( crystal_structure.find_label( neighbour_1_atom_label ) ).position();
        Vector3D origin_atom_frac = crystal_structure.atom( crystal_structure.find_label( origin_atom_label ) ).position();
        Vector3D neighbour_2_atom_frac = crystal_structure.atom( crystal_structure.find_label( neighbour_2_atom_label ) ).position();
        // Find shortest disance between the two taking periodicity and space-group symmetry into account.
        double distance;
        Vector3D r_1mino_frac;
        crystal_structure.shortest_distance( origin_atom_frac, neighbour_1_atom_frac, distance, r_1mino_frac );
        if ( distance > 3.0 )
            std::cout << "Warning: distance 1 > 3.0 A" << std::endl;
        Vector3D r_2mino_frac;
        crystal_structure.shortest_distance( origin_atom_frac, neighbour_2_atom_frac, distance, r_2mino_frac );
        if ( distance > 3.0 )
            std::cout << "Warning: distance 2 > 3.0 A" << std::endl;
        // We need to be able to set the length of the X-H bond, so we must work in Cartesian coordinates.
        Vector3D neighbour_1_atom_cart = crystal_structure.crystal_lattice().fractional_to_orthogonal( neighbour_1_atom_frac );
        Vector3D origin_atom_cart = crystal_structure.crystal_lattice().fractional_to_orthogonal( origin_atom_frac );
        Vector3D neighbour_2_atom_cart = crystal_structure.crystal_lattice().fractional_to_orthogonal( neighbour_2_atom_frac );
        Vector3D r_1mino_cart = crystal_structure.crystal_lattice().fractional_to_orthogonal( r_1mino_frac );
        NormalisedVector3D n_1mino_cart = normalised_vector( r_1mino_cart );
        Vector3D r_2mino_cart = crystal_structure.crystal_lattice().fractional_to_orthogonal( r_2mino_frac );
        NormalisedVector3D n_2mino_cart = normalised_vector( r_2mino_cart );
        NormalisedVector3D x = normalised_vector( -0.5*( n_1mino_cart + n_2mino_cart ) );
            // The order of the points determines the sign of the normal n.
        Plane plane( neighbour_1_atom_cart, origin_atom_cart, neighbour_2_atom_cart );
        NormalisedVector3D y = plane.normal();
        double target_bond_length( 1.0 );
        if ( crystal_structure.atom( crystal_structure.find_label( origin_atom_label ) ).element() == Element( "C" ) )
            target_bond_length = 1.089;
        else if ( crystal_structure.atom( crystal_structure.find_label( origin_atom_label ) ).element() == Element( "N" ) )
            target_bond_length = 1.015;
        else if ( crystal_structure.atom( crystal_structure.find_label( origin_atom_label ) ).element() == Element( "O" ) )
            target_bond_length = 0.993;
        Angle angle = Angle::from_degrees( 106.0 / 2.0 );
        double o = target_bond_length * angle.sine();
        double a = target_bond_length * angle.cosine();
        Vector3D H_atom_1_cart = origin_atom_cart + o*y + a*x;
        Vector3D H_atom_1_frac = crystal_structure.crystal_lattice().orthogonal_to_fractional( H_atom_1_cart );
        crystal_structure.add_atom( Atom( Element( "H" ), H_atom_1_frac, "H" + size_t2string( crystal_structure.natoms() ) ) );
        Vector3D H_atom_2_cart = origin_atom_cart - o*y + a*x;
        Vector3D H_atom_2_frac = crystal_structure.crystal_lattice().orthogonal_to_fractional( H_atom_2_cart );
        crystal_structure.add_atom( Atom( Element( "H" ), H_atom_2_frac, "H" + size_t2string( crystal_structure.natoms() ) ) );
        }

        crystal_structure.save_cif( append_to_file_name( input_file_name, "_H_added" ) );
    MACRO_END_GAME

    try // Write .inp for TLS from .cif + two _restraints.txt files.
    {
        if ( argc != 2 )
            throw std::runtime_error( "Please give the name of a .cif file that needs to be converted to a _TLS.inp file." );
        FileName input_file_name( argv[ 1 ] );
        TLSWriter( input_file_name );
     MACRO_END_GAME

    try // Write .inp for TLS from a .inp file.
    {
        if ( argc != 2 )
            throw std::runtime_error( "Please give the name of a .inp file that needs to be converted to _TLS.inp." );
        FileName input_file_name( argv[ 1 ] );
        TLSWriter_2( input_file_name );
     MACRO_END_GAME

    try // Write .inp from .cif + two _restraints.txt files + .xye file.
    {
        if ( argc != 3 )
            throw std::runtime_error( "Please give the name of a .cif file and a .xye file that need to be converted to a .inp file." );
        inp_writer( FileName( argv[ 1 ] ), FileName( argv[ 2 ] ) );
     MACRO_END_GAME

    try // Calculate all torsion angles and inversions for a list of cifs.
    {
        MACRO_ONE_FILELISTNAME_AS_ARGUMENT
        TextFileWriter text_file_writer( FileName( file_list_file_name.directory(), "torsions_and_inversions", "txt" ) );
        std::vector< std::vector< std::string > > torsions;

        std::vector< std::string > ring_labels;
        ring_labels.push_back( "C6" );
        ring_labels.push_back( "C15" );
        ring_labels.push_back( "C13" );
        ring_labels.push_back( "C12" );
        ring_labels.push_back( "C11" );
        ring_labels.push_back( "C10" );
        ring_labels.push_back( "N3" );
        ring_labels.push_back( "C7" );
        ring_labels.push_back( "C5" );
        ring_labels.push_back( "C2" );
        ring_labels.push_back( "C1" );
        ring_labels.push_back( "C0" );

        add_all_torsions_in_ring( ring_labels, torsions );

        // Note that for Z'=1 structures, we have e.g. "C1_0", for Z'=2 structures, we have e.g. "C1_0" and "C1_1".
//        std::vector< std::string > torsion;
//        torsion.push_back( "C5" );
//        torsion.push_back( "C6" );
//        torsion.push_back( "C7" );
//        torsion.push_back( "N3" );
//        torsions.push_back( torsion );
//        torsion.clear();
//        torsion.push_back( "C3" );
//        torsion.push_back( "S0" );
//        torsion.push_back( "N4" );
//        torsion.push_back( "C14" );
//        torsions.push_back( torsion );
//        torsion.clear();
//        torsion.push_back( "C0" );
//        torsion.push_back( "C1" );
//        torsion.push_back( "N4" );
//        torsion.push_back( "C13" );
//        torsions.push_back( torsion );
//
//        std::vector< std::string > ring_labels;
//        ring_labels.push_back( "N3" );
//        ring_labels.push_back( "C12" );
//        ring_labels.push_back( "C16" );
//        ring_labels.push_back( "N4" );
//        ring_labels.push_back( "C15" );
//        ring_labels.push_back( "C14" );
//        ring_labels.push_back( "C13" );
//        add_all_torsions_in_ring( ring_labels, torsions );

        std::vector< std::vector< std::string > > inversions;
        std::vector< std::string > inversion;
//        inversion.push_back( "N4" );
//        inversion.push_back( "S0" );
//        inversion.push_back( "C14" );
//        inversion.push_back( "C13" );
//        inversions.push_back( inversion );
//        inversion.clear();
//        inversion.push_back( "N3" );
//        inversion.push_back( "C11" );
//        inversion.push_back( "C12" );
//        inversion.push_back( "C14" );
//        inversions.push_back( inversion );
        std::string header_str( "Rank " );
        // Write header
        for ( size_t i( 0 ); i != torsions.size(); ++i )
        {
            header_str += "t_" + torsions[i][0] + "_";
            header_str += torsions[i][1] + "_";
            header_str += torsions[i][2] + "_";
            header_str += torsions[i][3] + " ";
        }
        for ( size_t i( 0 ); i != inversions.size(); ++i )
        {
            header_str += "i_" + inversions[i][0] + "_";
            header_str += inversions[i][1] + "_";
            header_str += inversions[i][2] + "_";
            header_str += inversions[i][3] + " ";
        }
        text_file_writer.write_line( header_str );
        for ( size_t i( 0 ); i != file_list.size(); ++i )
        {
            CrystalStructure crystal_structure;
            std::cout << "Now reading cif... " + file_list.value( i ).full_name() << std::endl;
            read_cif( file_list.value( i ), crystal_structure );
            CrystalLattice crystal_lattice = crystal_structure.crystal_lattice();
            bool Zprime_is_two( false );
            if ( torsions.empty() )
            {
                try
                {
                    /*size_t index =*/ crystal_structure.atom( inversions[0][0] + "_1" );
                    Zprime_is_two = true;
                }
                catch ( std::exception & e ) {}
            }
            else
            {
                try
                {
                    /* size_t index =*/ crystal_structure.atom( torsions[0][0] + "_1" );
                    Zprime_is_two = true;
                }
                catch ( std::exception & e ) {}
            }
            std::string output_Z1_string( size_t2string( i+1 ) + " " );
            std::string output_Z2_string( size_t2string( i+1 ) + " " );
            // Torsions
            for ( size_t i( 0 ); i != torsions.size(); ++i )
            {
                size_t i1 = crystal_structure.atom( torsions[i][0] + "_0" );
                Atom atom1 = crystal_structure.atom( i1 );
                Vector3D r1 = crystal_lattice.fractional_to_orthogonal( atom1.position() );
                size_t i2 = crystal_structure.atom( torsions[i][1] + "_0" );
                Atom atom2 = crystal_structure.atom( i2 );
                Vector3D r2 = crystal_lattice.fractional_to_orthogonal( atom2.position() );
                size_t i3 = crystal_structure.atom( torsions[i][2] + "_0" );
                Atom atom3 = crystal_structure.atom( i3 );
                Vector3D r3 = crystal_lattice.fractional_to_orthogonal( atom3.position() );
                size_t i4 = crystal_structure.atom( torsions[i][3] + "_0" );
                Atom atom4 = crystal_structure.atom( i4 );
                Vector3D r4 = crystal_lattice.fractional_to_orthogonal( atom4.position() );
                Angle st = signed_torsion( r1, r2, r3, r4 );
                output_Z1_string += double2string( st.value_in_degrees() ) + " ";
            }
            if ( Zprime_is_two )
            {
                for ( size_t i( 0 ); i != torsions.size(); ++i )
                {
                    size_t i1 = crystal_structure.atom( torsions[i][0] + "_1" );
                    Atom atom1 = crystal_structure.atom( i1 );
                    Vector3D r1 = crystal_lattice.fractional_to_orthogonal( atom1.position() );
                    size_t i2 = crystal_structure.atom( torsions[i][1] + "_1" );
                    Atom atom2 = crystal_structure.atom( i2 );
                    Vector3D r2 = crystal_lattice.fractional_to_orthogonal( atom2.position() );
                    size_t i3 = crystal_structure.atom( torsions[i][2] + "_1" );
                    Atom atom3 = crystal_structure.atom( i3 );
                    Vector3D r3 = crystal_lattice.fractional_to_orthogonal( atom3.position() );
                    size_t i4 = crystal_structure.atom( torsions[i][3] + "_1" );
                    Atom atom4 = crystal_structure.atom( i4 );
                    Vector3D r4 = crystal_lattice.fractional_to_orthogonal( atom4.position() );
                    Angle st = signed_torsion( r1, r2, r3, r4 );
                    output_Z2_string += double2string( st.value_in_degrees() ) + " ";
                }
            }
            for ( size_t i( 0 ); i != inversions.size(); ++i )
            {
                size_t i1 = crystal_structure.atom( inversions[i][0] + "_0" );
                Atom atom1 = crystal_structure.atom( i1 );
                Vector3D r1 = crystal_lattice.fractional_to_orthogonal( atom1.position() );
                size_t i2 = crystal_structure.atom( inversions[i][1] + "_0" );
                Atom atom2 = crystal_structure.atom( i2 );
                Vector3D r2 = crystal_lattice.fractional_to_orthogonal( atom2.position() );
                size_t i3 = crystal_structure.atom( inversions[i][2] + "_0" );
                Atom atom3 = crystal_structure.atom( i3 );
                Vector3D r3 = crystal_lattice.fractional_to_orthogonal( atom3.position() );
                size_t i4 = crystal_structure.atom( inversions[i][3] + "_0" );
                Atom atom4 = crystal_structure.atom( i4 );
                Vector3D r4 = crystal_lattice.fractional_to_orthogonal( atom4.position() );
                Plane plane( r2, r3, r4 );
                double sd = plane.signed_distance( r1 );
                output_Z1_string += double2string( sd ) + " ";
            }
            if ( Zprime_is_two )
            {
                for ( size_t i( 0 ); i != inversions.size(); ++i )
                {
                    size_t i1 = crystal_structure.atom( inversions[i][0] + "_1" );
                    Atom atom1 = crystal_structure.atom( i1 );
                    Vector3D r1 = crystal_lattice.fractional_to_orthogonal( atom1.position() );
                    size_t i2 = crystal_structure.atom( inversions[i][1] + "_1" );
                    Atom atom2 = crystal_structure.atom( i2 );
                    Vector3D r2 = crystal_lattice.fractional_to_orthogonal( atom2.position() );
                    size_t i3 = crystal_structure.atom( inversions[i][2] + "_1" );
                    Atom atom3 = crystal_structure.atom( i3 );
                    Vector3D r3 = crystal_lattice.fractional_to_orthogonal( atom3.position() );
                    size_t i4 = crystal_structure.atom( inversions[i][3] + "_1" );
                    Atom atom4 = crystal_structure.atom( i4 );
                    Vector3D r4 = crystal_lattice.fractional_to_orthogonal( atom4.position() );
                    Plane plane( r2, r3, r4 );
                    double sd = plane.signed_distance( r1 );
                    output_Z2_string += double2string( sd ) + " ";
                }
            }
            text_file_writer.write_line( output_Z1_string );
            if ( Zprime_is_two )
                text_file_writer.write_line( output_Z2_string );
        }
    MACRO_END_GAME

    try // CrystalStructure::supercell().
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        size_t u( 1 );
        size_t v( 1 );
        size_t w( 2 );
        crystal_structure.supercell( u, v, w );
        crystal_structure.make_atom_labels_unique();
        crystal_structure.save_cif( append_to_file_name( input_file_name, "_" + size_t2string( u ) + "_" + size_t2string( v ) + "_" + size_t2string( w ) ) );
    MACRO_END_GAME

    try // Calculate similarity matrix.
    {
        MACRO_ONE_FILELISTNAME_AS_ARGUMENT
        CorrelationMatrix similarity_matrix = calculate_correlation_matrix( file_list );
        similarity_matrix.save( FileName( file_list_file_name.directory(), "SimilarityMatrix", "txt" ) );
    MACRO_END_GAME

    try // Split reflections from SHELX .hkl file into +(hkl) and -(hkl).
    {
        // We need a CrystalLattice and a space group
        FileName input_file_name( "C:\\Data_Win\\Research\\Refereeing\\bm5114sup5_2.cif" );
        CrystalStructure crystal_structure;
        std::cout << "Now reading cif... " + input_file_name.full_name() << std::endl;
        read_cif( input_file_name, crystal_structure );
        crystal_structure.apply_space_group_symmetry();
        ReflectionList reflection_list;
        reflection_list.read_hkl( FileName( "C:\\Data_Win\\Research\\Refereeing\\bm5114sup6.hkl" ) );
        std::cout << ".hkl file has been read." << std::endl;
        std::vector< bool > done( reflection_list.size(), false );
        TextFileWriter text_file_writer( FileName( "C:\\Data_Win\\Research\\Refereeing\\bm5114sup6_p_m.txt" ) );
        for ( size_t i( 0 ); i != reflection_list.size(); ++i )
        {
            if ( done[i] )
                continue;
            MillerIndices miller_indices_p = reflection_list.miller_indices( i );
            MillerIndices miller_indices_m = -miller_indices_p;
            double F_squared_p = 0.0;
            double F_squared_m = 0.0;
            size_t nreflections_p( 0 );
            size_t nreflections_m( 0 );
            // It is possible that the same reflection has been measured multiple times, so we cannot simply find "a" reflection based on
            // a set of hkl indices, as the list may contain the same set of hkl indices multiple times. Therefore, the only way is to go through the
            // entire list from start to finish every single time.
            for ( size_t j( 0 ); j != reflection_list.size(); ++j )
            {
                if ( done[j] )
                    continue;
                if ( reflection_list.miller_indices( j ) == miller_indices_p )
                {
                    ++nreflections_p;
                    done[j] = true;
                    F_squared_p += reflection_list.F_squared( j );
                }
                else if ( reflection_list.miller_indices( j ) == miller_indices_m )
                {
                    ++nreflections_m;
                    done[j] = true;
                    F_squared_m += reflection_list.F_squared( j );
                }
            }
            if ( ( nreflections_p == 0 ) || ( nreflections_m == 0 ) )
            {
                std::cout << "WARNING: nreflections_m == 0 for " + miller_indices_m.to_string() << std::endl;
                continue;
            }
            F_squared_p /= nreflections_p;
            if ( F_squared_p < 0.0 )
                F_squared_p = 0.0;
            F_squared_m /= nreflections_m;
            if ( F_squared_m < 0.0 )
                F_squared_m = 0.0;
            if ( miller_indices_m < miller_indices_p )
                text_file_writer.write_line( miller_indices_p.to_string() + " " + double2string( F_squared_p ) + " " + double2string( F_squared_m ) );
            else
                text_file_writer.write_line( miller_indices_m.to_string() + " " + double2string( F_squared_m ) + " " + double2string( F_squared_p ) );
        }
    MACRO_END_GAME

    // Change bond length.
    try
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        {
        std::string origin_atom_label( "C9_0" );
        std::string neighbour_atom_label( "H10_0" );
        Vector3D origin_atom_frac = crystal_structure.atom( crystal_structure.find_label( origin_atom_label ) ).position();
        Vector3D neighbour_atom_frac = crystal_structure.atom( crystal_structure.find_label( neighbour_atom_label ) ).position();
        // Find shortest disance between the two taking periodicity and space-group symmetry into account.
        double distance;
        Vector3D difference_frac;
        crystal_structure.shortest_distance( origin_atom_frac, neighbour_atom_frac, distance, difference_frac );
        if ( distance > 3.0 )
            std::cout << "Warning: distance > 3.0 A" << std::endl;
        // We need to be able to set the length of the X-H bond, so we must work in Cartesian coordinates.
        Vector3D difference_cart = crystal_structure.crystal_lattice().fractional_to_orthogonal( difference_frac );
        // C-H 1.089
        // N-H 1.015
        // O-H 0.993
        // C-F 1.37
        double target_bond_length( 1.37 );
        difference_cart *= target_bond_length / distance;
        difference_frac = crystal_structure.crystal_lattice().orthogonal_to_fractional( difference_cart );
        Atom new_atom( crystal_structure.atom( crystal_structure.find_label( neighbour_atom_label ) ) );
        new_atom.set_position( origin_atom_frac + difference_frac );
    //    new_atom.set_element( Element( "F" ) );
        crystal_structure.set_atom( crystal_structure.find_label( neighbour_atom_label ), new_atom );
        }
        {
        std::string origin_atom_label( "C11_0" );
        std::string neighbour_atom_label( "H12_0" );
        Vector3D origin_atom_frac = crystal_structure.atom( crystal_structure.find_label( origin_atom_label ) ).position();
        Vector3D neighbour_atom_frac = crystal_structure.atom( crystal_structure.find_label( neighbour_atom_label ) ).position();
        // Find shortest disance between the two taking periodicity and space-group symmetry into account.
        double distance;
        Vector3D difference_frac;
        crystal_structure.shortest_distance( origin_atom_frac, neighbour_atom_frac, distance, difference_frac );
        if ( distance > 3.0 )
            std::cout << "Warning: distance > 3.0 A" << std::endl;
        // We need to be able to set the length of the X-H bond, so we must work in Cartesian coordinates.
        Vector3D difference_cart = crystal_structure.crystal_lattice().fractional_to_orthogonal( difference_frac );
        // C-H 1.089
        // N-H 1.015
        // O-H 0.993
        // C-F 1.37
        double target_bond_length( 1.37 );
        difference_cart *= target_bond_length / distance;
        difference_frac = crystal_structure.crystal_lattice().orthogonal_to_fractional( difference_cart );
        Atom new_atom( crystal_structure.atom( crystal_structure.find_label( neighbour_atom_label ) ) );
        new_atom.set_position( origin_atom_frac + difference_frac );
    //    new_atom.set_element( Element( "F" ) );
        crystal_structure.set_atom( crystal_structure.find_label( neighbour_atom_label ), new_atom );
        }
        crystal_structure.save_cif( append_to_file_name( input_file_name, "_blc" ) );
    MACRO_END_GAME

    try // Generate a powder diffraction pattern from SHELX .hkl file.
    {
        // We need a CrystalLattice and a space group
        FileName input_file_name( "C:\\Data_Win\\Research\\Refereeing\\bm5114sup5_2.cif" );
        CrystalStructure crystal_structure;
        std::cout << "Now reading cif... " + input_file_name.full_name() << std::endl;
        read_cif( input_file_name, crystal_structure );
        crystal_structure.apply_space_group_symmetry();
        ReflectionList reflection_list;
        reflection_list.read_hkl( FileName( "C:\\Data_Win\\Research\\Refereeing\\bm5114sup6.hkl" ) );
        std::cout << ".hkl file has been read." << std::endl;
        PointGroup laue_class = crystal_structure.space_group().laue_class();
        ReflectionList reflection_list_final;
        std::vector< bool > done( reflection_list.size(), false );
        for ( size_t i( 0 ); i != reflection_list.size(); ++i )
        {
            if ( done[i] )
                continue;
            MillerIndices miller_indices = reflection_list.miller_indices( i );
            // Calculate equivalent reflections
            std::set< MillerIndices > equivalent_reflections;
            equivalent_reflections.insert( miller_indices );
            for ( size_t j( 0 ); j != laue_class.nsymmetry_operators(); ++j )
                equivalent_reflections.insert( miller_indices * laue_class.symmetry_operator( j ) );
            double F_squared = 0.0;
            size_t number_of_equivalent_reflections_that_have_actually_been_measured( 0 );
            // It is possible that the same reflection has been measured multiple times, so we cannot simply find "a" reflection based on
            // a set of hkl indices, as the list may contain the same set of hkl indices multiple times. Therefore, the only way is to go through the
            // entire list from start to finish every single time.
            for ( size_t j( 0 ); j != reflection_list.size(); ++j )
            {
                if ( done[j] )
                    continue;
                // Is this reflection one of the equivalent reflections?
                for ( std::set< MillerIndices >::const_iterator it = equivalent_reflections.begin(); it != equivalent_reflections.end(); ++it )
                {
                    if ( *it == reflection_list.miller_indices( j ) )
                    {
                        ++number_of_equivalent_reflections_that_have_actually_been_measured;
                        done[j] = true;
                        F_squared += reflection_list.F_squared( j );
                        break;
                    }
                }
            }
            F_squared /= number_of_equivalent_reflections_that_have_actually_been_measured;
            if ( F_squared < 0.0 )
                F_squared = 0.0;
            Vector3D H = miller_indices.h() * crystal_structure.crystal_lattice().a_star_vector() +
                         miller_indices.k() * crystal_structure.crystal_lattice().b_star_vector() +
                         miller_indices.l() * crystal_structure.crystal_lattice().c_star_vector();
            double d = 1.0 / ( H.length() );
            reflection_list_final.push_back( miller_indices, F_squared, d, equivalent_reflections.size() );
        }
        reflection_list_final.save( FileName( "C:\\Data_Win\\Research\\Refereeing\\bm5114sup6_cal.hkl" ) );
        PowderPatternCalculator powder_pattern_calculator( crystal_structure );
        powder_pattern_calculator.set_FWHM( 0.05 );
        powder_pattern_calculator.set_two_theta_end( Angle::from_degrees( 40.0 ) );
        powder_pattern_calculator.set_two_theta_step( Angle::from_degrees( 0.01 ) );
        PowderPattern powder_pattern;
        powder_pattern_calculator.calculate( reflection_list_final, powder_pattern );
        powder_pattern.save_xye( FileName( "C:\\Data_Win\\Research\\Refereeing\\bm5114sup6_cal.xye" ), true );
    MACRO_END_GAME

    try // Calculate RMSCD for two FileList.txt files without matching.
    {
        if ( argc != 3 )
            throw std::runtime_error( "Please give the names of two FileList.txt files." );
        FileName file_list_file_name_1( argv[ 1 ] );
        FileList file_list_1( file_list_file_name_1 );
        if ( file_list_1.empty() )
            throw std::runtime_error( std::string( "Error: No files in file list " ) + file_list_file_name_1.full_name() );
        FileName file_list_file_name_2( argv[ 2 ] );
        FileList file_list_2( file_list_file_name_2 );
        if ( file_list_2.empty() )
            throw std::runtime_error( std::string( "Error: No files in file list " ) + file_list_file_name_2.full_name() );
        if ( file_list_1.size() != file_list_2.size() )
            throw std::runtime_error( "Error: FileList.txt files do not contain the same number of entries." );
        TextFileWriter text_file_writer( FileName( file_list_file_name_1.directory(), "RMSCDs", "txt" ) );
        text_file_writer.write_line( file_list_file_name_1.full_name() );
        text_file_writer.write_line( file_list_file_name_2.full_name() );
        for ( size_t fi( 0 ); fi != file_list_1.size(); ++fi )
        {
            FileName file_name_1( file_list_1.value( fi ) );
            CrystalStructure crystal_structure_1;
            std::cout << "Now reading cif... " + file_name_1.full_name() << std::endl;
            read_cif( file_name_1, crystal_structure_1 );
            FileName file_name_2( file_list_2.value( fi ) );
            CrystalStructure crystal_structure_2;
            std::cout << "Now reading cif... " + file_name_2.full_name() << std::endl;
            read_cif( file_name_2, crystal_structure_2 );
//            double result = RMSCD_with_matching( crystal_structure_1, crystal_structure_2, false );
            double result = root_mean_square_Cartesian_displacement( crystal_structure_1, crystal_structure_2, false );
            text_file_writer.write_line( file_list_1.value( fi ).file_name() + " " + file_list_2.value( fi ).file_name() + " " + double2string( result ) );
        }
    MACRO_END_GAME

    try // *** doubles *** Generate coefficients for Chebyshev polynomial of the first kind.
    {
        size_t max_order = 50;
        std::vector< std::vector< double > > vector_of_coefficients;
        std::vector< double > coefficients_0;   // All zeros
        std::vector< double > coefficients_n;   // Order n
        std::vector< double > coefficients_n_1; // Order n - 1
        std::vector< double > coefficients_n_2; // Order n - 2
        for ( size_t j( 0 ); j != max_order+1; ++j )
            coefficients_0.push_back( 0 );
        coefficients_n = coefficients_0;
        coefficients_n[0] = 1;
        vector_of_coefficients.push_back( coefficients_n );
        coefficients_n = coefficients_0;
        coefficients_n[1] = 1;
        vector_of_coefficients.push_back( coefficients_n );
        for ( size_t i( 2 ); i != max_order+1; ++i )
        {
            coefficients_n = coefficients_0;
            coefficients_n_1 = vector_of_coefficients[ i - 1 ];
            coefficients_n_2 = vector_of_coefficients[ i - 2 ];
            for ( size_t j( 0 ); j != max_order; ++j )
            {
                coefficients_n[j+1] += 2*coefficients_n_1[j];
            }
            for ( size_t j( 0 ); j != max_order+1; ++j )
            {
                coefficients_n[j] -= coefficients_n_2[j];
            }
            vector_of_coefficients.push_back( coefficients_n );
        }

        std::cout << "double Chebyshev( const size_t order, const double x )" << std::endl;
        std::cout << "{" << std::endl;
        std::cout << "    switch( order )" << std::endl;
        std::cout << "    {" << std::endl;
        for ( size_t i( 0 ); i != max_order+1; ++i )
        {
            std::cout << "        case  " << size_t2string( i , 2, ' ' ) << " : return";
            coefficients_n = vector_of_coefficients[i];
            for ( size_t j( 0 ); j != max_order+1; ++j )
            {
                if ( coefficients_n[j] == 0 )
                    continue;
                if ( coefficients_n[j] < 0 )
                    std::cout << " - ";
                else
                    std::cout << " + ";
                std::cout << std::abs( coefficients_n[j] ) << ".0";
                if ( j == 1 )
                    std::cout << " * x";
                else if ( j > 1 )
                    std::cout << " * std::pow( x, " << j << " )";
            }
            std::cout << ";" << std::endl;
        }
        std::cout << "    }" << std::endl;
        std::cout << "}" << std::endl;
    MACRO_END_GAME

    try // Generate coefficients for Chebyshev polynomial of the first kind.
    {
        size_t max_order = 50;
        std::vector< std::vector< long int > > vector_of_coefficients;
        std::vector< long int > coefficients_0;   // All zeros
        std::vector< long int > coefficients_n;   // Order n
        std::vector< long int > coefficients_n_1; // Order n - 1
        std::vector< long int > coefficients_n_2; // Order n - 2
        for ( size_t j( 0 ); j != max_order+1; ++j )
            coefficients_0.push_back( 0 );
        coefficients_n = coefficients_0;
        coefficients_n[0] = 1;
        vector_of_coefficients.push_back( coefficients_n );
        coefficients_n = coefficients_0;
        coefficients_n[1] = 1;
        vector_of_coefficients.push_back( coefficients_n );
        for ( size_t i( 2 ); i != max_order+1; ++i )
        {
            coefficients_n = coefficients_0;
            coefficients_n_1 = vector_of_coefficients[ i - 1 ];
            coefficients_n_2 = vector_of_coefficients[ i - 2 ];
            for ( size_t j( 0 ); j != max_order; ++j )
            {
                coefficients_n[j+1] += 2*coefficients_n_1[j];
            }
            for ( size_t j( 0 ); j != max_order+1; ++j )
            {
                coefficients_n[j] -= coefficients_n_2[j];
            }
            vector_of_coefficients.push_back( coefficients_n );
        }

        std::cout << "double Chebyshev( const size_t order, const double x )" << std::endl;
        std::cout << "{" << std::endl;
        std::cout << "    switch( order )" << std::endl;
        std::cout << "    {" << std::endl;
        for ( size_t i( 0 ); i != max_order+1; ++i )
        {
            std::cout << "        case  " << size_t2string( i , 2, ' ' ) << " : return";
            coefficients_n = vector_of_coefficients[i];
            for ( size_t j( 0 ); j != max_order+1; ++j )
            {
                if ( coefficients_n[j] == 0 )
                    continue;
                if ( coefficients_n[j] < 0 )
                    std::cout << " - ";
                else
                    std::cout << " + ";
                std::cout << std::abs( coefficients_n[j] ) << ".0";
                if ( j == 1 )
                    std::cout << " * x";
                else if ( j > 1 )
                    std::cout << " * std::pow( x, " << j << " )";
            }
            std::cout << ";" << std::endl;
        }
        std::cout << "    }" << std::endl;
        std::cout << "}" << std::endl;
    MACRO_END_GAME

    // Apply space-group symmetry.
    try
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        crystal_structure.apply_space_group_symmetry();
        crystal_structure.save_cif( append_to_file_name( input_file_name, "_asgs" ) );
    MACRO_END_GAME

    try // Invert structure by multiplying all coordinates by -1.0.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        Matrix3D rotation( -1.0,  0.0,  0.0,
                            0.0, -1.0,  0.0,
                            0.0,  0.0, -1.0 );
        Vector3D shift( 1.0, 1.0, 1.0 );
        for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
        {
            Atom new_atom( crystal_structure.atom( i ) );
            new_atom.set_position( ( rotation * crystal_structure.atom( i ).position() ) + shift );
            crystal_structure.set_atom( i, new_atom );
        }
        crystal_structure.save_cif( append_to_file_name( input_file_name, "_inverted" ) );
    MACRO_END_GAME

    // Add hydrogen atom to sp2 ring carbon.
    try
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        {
        std::string atom_C_label( "C114" );
        std::string neighbour_1_label( "C60" );
        std::string neighbour_2_label( "N88" );
        Vector3D atom_C = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_C_label ) ).position() );
        Vector3D neighbour_1 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( neighbour_1_label ) ).position() );
        Vector3D neighbour_2 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( neighbour_2_label ) ).position() );
        Vector3D hydrogen_atom = add_hydrogen_atom_to_sp2_atom( atom_C, Element( "C" ), neighbour_1, neighbour_2 );
        crystal_structure.add_atom( Atom( Element( "H" ), crystal_structure.crystal_lattice().orthogonal_to_fractional( hydrogen_atom ), "H" + size_t2string( crystal_structure.natoms() ) ) );
        }
        crystal_structure.save_cif( append_to_file_name( input_file_name, "_CH_added" ) );
    MACRO_END_GAME

    // Correct carbon atom in aromatic ring.
    try
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        {
        std::string atom_1_label( "C60" );
        std::string atom_2_label( "C65" );
        std::string atom_3_label( "C74" );
        std::string atom_4_label( "N88" );
        std::string atom_5_label( "C97" );

        std::string atom_0_label( "C66" );

        Vector3D atom_0 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_0_label ) ).position() );

        Vector3D atom_1 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_1_label ) ).position() );
        Vector3D atom_2 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_2_label ) ).position() );
        Vector3D atom_3 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_3_label ) ).position() );
        Vector3D atom_4 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_4_label ) ).position() );
        Vector3D atom_5 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_5_label ) ).position() );
        std::vector< Vector3D > points;
        points.push_back( atom_1 );
        points.push_back( atom_2 );
        points.push_back( atom_3 );
        points.push_back( atom_4 );
        points.push_back( atom_5 );
        Plane plane( points );

        double distance = plane.signed_distance( atom_0 );
        Vector3D new_position = atom_0 + distance * plane.normal();
        crystal_structure.add_atom( Atom( Element( "C" ), crystal_structure.crystal_lattice().orthogonal_to_fractional( new_position ), "C" + size_t2string( crystal_structure.natoms() ) ) );
        }
        crystal_structure.save_cif( append_to_file_name( input_file_name, "_C_added" ) );
    MACRO_END_GAME

    // Add hydrogen atoms to sp3 nitrogen with hydrogen bond.
    try
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        {
        std::string atom_N_label( "N08" );
        std::string atom_2_label( "C00Q" );
        std::string atom_H_label( "H114" );
        Vector3D atom_N = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_N_label ) ).position() );
        Vector3D atom_2 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_2_label ) ).position() );
        double distance;
        Vector3D difference_vector;
        crystal_structure.shortest_distance( crystal_structure.atom( crystal_structure.find_label( atom_N_label ) ).position(), crystal_structure.atom( crystal_structure.find_label( atom_H_label ) ).position(), distance, difference_vector );
        Vector3D atom_H = crystal_structure.atom( crystal_structure.find_label( atom_N_label ) ).position() + difference_vector; // Fractional coordinates
        atom_H = crystal_structure.crystal_lattice().fractional_to_orthogonal( atom_H );
        std::vector< Vector3D > hydrogen_atoms = add_2_hydrogen_atoms_to_sp3_nitrogen( atom_N, atom_2, atom_H );
        crystal_structure.add_atom( Atom( Element( "H" ), crystal_structure.crystal_lattice().orthogonal_to_fractional( hydrogen_atoms[0] ), "H" + size_t2string( crystal_structure.natoms() ) ) );
        crystal_structure.add_atom( Atom( Element( "H" ), crystal_structure.crystal_lattice().orthogonal_to_fractional( hydrogen_atoms[1] ), "H" + size_t2string( crystal_structure.natoms() ) ) );
        }
        {
        std::string atom_N_label( "N03" );
        std::string atom_2_label( "C01A" );
        std::string atom_H_label( "H14" );
        Vector3D atom_N = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_N_label ) ).position() );
        Vector3D atom_2 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_2_label ) ).position() );
        double distance;
        Vector3D difference_vector;
        crystal_structure.shortest_distance( crystal_structure.atom( crystal_structure.find_label( atom_N_label ) ).position(), crystal_structure.atom( crystal_structure.find_label( atom_H_label ) ).position(), distance, difference_vector );
        Vector3D atom_H = crystal_structure.atom( crystal_structure.find_label( atom_N_label ) ).position() + difference_vector; // Fractional coordinates
        atom_H = crystal_structure.crystal_lattice().fractional_to_orthogonal( atom_H );
        std::vector< Vector3D > hydrogen_atoms = add_2_hydrogen_atoms_to_sp3_nitrogen( atom_N, atom_2, atom_H );
        crystal_structure.add_atom( Atom( Element( "H" ), crystal_structure.crystal_lattice().orthogonal_to_fractional( hydrogen_atoms[0] ), "H" + size_t2string( crystal_structure.natoms() ) ) );
        crystal_structure.add_atom( Atom( Element( "H" ), crystal_structure.crystal_lattice().orthogonal_to_fractional( hydrogen_atoms[1] ), "H" + size_t2string( crystal_structure.natoms() ) ) );
        }
        crystal_structure.save_cif( append_to_file_name( input_file_name, "_NH2_added" ) );
    MACRO_END_GAME

    try // Take a disordered atom that has been modelled as a large ADP and change it into a split-atom model.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        // C01O
        {
        Atom atom = crystal_structure.atom( crystal_structure.find_label( "C01O" ) );
        SymmetricMatrix3D C01O_Ucart = atom.anisotropic_displacement_parameters().U_cart();
        std::vector< double > eigenvalues;
        std::vector< NormalisedVector3D > eigenvectors;
        calculate_eigenvalues( C01O_Ucart, eigenvalues, eigenvectors );
        Vector3D new_position_1 = crystal_structure.crystal_lattice().fractional_to_orthogonal( atom.position() ) + 0.5 * eigenvectors[2];
        new_position_1 = crystal_structure.crystal_lattice().orthogonal_to_fractional( new_position_1 );
        Atom new_atom_1( Element( "C" ), new_position_1, "C1a" );
        crystal_structure.add_atom( new_atom_1 );
        Vector3D new_position_2 = crystal_structure.crystal_lattice().fractional_to_orthogonal( atom.position() ) - 0.5 * eigenvectors[2];
        new_position_2 = crystal_structure.crystal_lattice().orthogonal_to_fractional( new_position_2 );
        Atom new_atom_2( Element( "C" ), new_position_2, "C1b" );
        crystal_structure.add_atom( new_atom_2 );
        }
        // C01N
        {
        Atom atom = crystal_structure.atom( crystal_structure.find_label( "O1" ) );
        SymmetricMatrix3D C01N_Ucart = atom.anisotropic_displacement_parameters().U_cart();
        std::vector< double > eigenvalues;
        std::vector< NormalisedVector3D > eigenvectors;
        calculate_eigenvalues( C01N_Ucart, eigenvalues, eigenvectors );
        Vector3D new_position_1 = crystal_structure.crystal_lattice().fractional_to_orthogonal( atom.position() ) + 0.5 * eigenvectors[2];
        new_position_1 = crystal_structure.crystal_lattice().orthogonal_to_fractional( new_position_1 );
        Atom new_atom_1( Element( "O" ), new_position_1, "O1a" );
        crystal_structure.add_atom( new_atom_1 );
        Vector3D new_position_2 = crystal_structure.crystal_lattice().fractional_to_orthogonal( atom.position() ) - 0.5 * eigenvectors[2];
        new_position_2 = crystal_structure.crystal_lattice().orthogonal_to_fractional( new_position_2 );
        Atom new_atom_2( Element( "O" ), new_position_2, "O1b" );
        crystal_structure.add_atom( new_atom_2 );
        }
        crystal_structure.save_cif( append_to_file_name( input_file_name, "_split" ) );
    MACRO_END_GAME

    try // Analyse flexible five- and six-membered rings for a list of cifs.
    {
        MACRO_ONE_FILELISTNAME_AS_ARGUMENT
        TextFileWriter text_file_writer( FileName( file_list_file_name.directory(), file_list_file_name.file_name() + "_flexible_rings", "txt" ) );
        std::vector< std::vector< std::string > > rings_5;
        // Note that for Z'=1 structures, we have e.g. "C1_0", for Z'=2 structures, we have e.g. "C1_0" and "C1_1".
        std::vector< std::string > ring;
        ring.push_back( "C18" );
        ring.push_back( "C17" );
        ring.push_back( "C13" );
        ring.push_back( "C16" );
        ring.push_back( "O1" );
        rings_5.push_back( ring );
//        ring.clear();
//        ring.push_back( "C2" );
//        ring.push_back( "C12" );
//        ring.push_back( "C2" );
//        ring.push_back( "N1" );
//        ring.push_back( "C9" );
//        rings_5.push_back( ring );
        ring.clear();
        std::vector< std::vector< std::string > > rings_6;
//        ring.push_back( "N2" );
//        ring.push_back( "C10" );
//        ring.push_back( "C12" );
//        ring.push_back( "N3" );
//        ring.push_back( "C11" );
//        ring.push_back( "C9" );
//        rings_6.push_back( ring );
//        ring.clear();
//        ring.push_back( "C5" );
//        ring.push_back( "C8" );
//        ring.push_back( "N1" );
//        ring.push_back( "C9" );
//        rings_6.push_back( ring );
//        ring.clear();
        if ( rings_5.empty() && rings_6.empty() )
            throw std::runtime_error( "Analyse flexible five- and six-membered rings for a list of cifs: nothing to do." );
        std::string header_str;
        // Write header
        for ( size_t i( 0 ); i != rings_5.size(); ++i )
        {
            header_str += "r_" + rings_5[i][0] + "_";
            header_str += rings_5[i][1] + "_";
            header_str += rings_5[i][2] + "_";
            header_str += rings_5[i][3] + "_";
            header_str += rings_5[i][4] + " ";
        }
        for ( size_t i( 0 ); i != rings_6.size(); ++i )
        {
            header_str += "r_" + rings_6[i][0] + "_";
            header_str += rings_6[i][1] + "_";
            header_str += rings_6[i][2] + "_";
            header_str += rings_6[i][3] + "_";
            header_str += rings_6[i][4] + "_";
            header_str += rings_6[i][5] + " ";
        }
        text_file_writer.write_line( header_str );

        TextFileWriter output_file( FileName( "C:\\Data_Win\\ring_analyser_testing.txt" ) );
        TextFileWriter distance_from_plane_file( FileName( "C:\\Data_Win\\distance_from_plane.txt" ) );
        TextFileWriter planarity_file( FileName( "C:\\Data_Win\\ring_analyser_planarity.txt" ) );
        TextFileWriter diff_1st_2_file( FileName( "C:\\Data_Win\\diff_1st_2.txt" ) );

        for ( size_t fi( 0 ); fi != file_list.size(); ++fi )
        {
            CrystalStructure crystal_structure;
            std::cout << "Now reading cif... " + file_list.value( fi ).full_name() << std::endl;
            std::string analysis_01;
            std::string analysis_02;
            std::string output_string = file_list.value( fi ).file_name();
            read_cif( file_list.value( fi ), crystal_structure );
            CrystalLattice crystal_lattice = crystal_structure.crystal_lattice();
            bool Zprime_is_two( false );
            if ( rings_5.empty() )
            {
                try
                {
                    /*size_t index =*/ crystal_structure.atom( rings_6[0][0] + "_1" );
                    Zprime_is_two = true;
                }
                catch ( std::exception & e ) {}
            }
            else
            {
                try
                {
                    /*size_t index =*/ crystal_structure.atom( rings_5[0][0] + "_1" );
                    Zprime_is_two = true;
                }
                catch ( std::exception & e ) {}
            }
            FiveMemberedRingAnalyser five_membered_ring_analyser;
            for ( size_t j( 0 ); j != rings_5.size(); ++j )
            {
                std::vector< Vector3D > points;
                for ( size_t k( 0 ); k != rings_5[j].size(); ++k )
                {
                    size_t index = crystal_structure.atom( rings_5[j][k] + "_0" );
                    Atom atom = crystal_structure.atom( index );
                    Vector3D point = crystal_lattice.fractional_to_orthogonal( atom.position() );
                    points.push_back( point );
                }
                five_membered_ring_analyser.analyse( points );

                if ( ! five_membered_ring_analyser.is_planar() )
                {
                    output_file.write( file_list.value( fi ).file_name() + " " );
                    for ( size_t k( 0 ); k != 5; ++k )
                    {
                        output_file.write( double2string( five_membered_ring_analyser.distances_from_plane_[ five_membered_ring_analyser.sorted_map_[ 4-k ] ] ) + " " );
                    }
                    output_file.write_line();
                    output_file.write( file_list.value( fi ).file_name() + " " );
                    for ( size_t k( 0 ); k != 5; ++k )
                        output_file.write( double2string( five_membered_ring_analyser.rmsds_from_mean_plane_[ five_membered_ring_analyser.sorted_map_[ 4-k ] ] ) + " " );
                    output_file.write_line();
                    distance_from_plane_file.write( file_list.value( fi ).file_name() + " " );
                    distance_from_plane_file.write_line( double2string( five_membered_ring_analyser.distances_from_plane_[ five_membered_ring_analyser.sorted_map_[ 4 ] ] ) );
                    planarity_file.write( file_list.value( fi ).file_name() + " " );
                    planarity_file.write_line( double2string( five_membered_ring_analyser.rmsds_from_mean_plane_[ five_membered_ring_analyser.sorted_map_[ 4 ] ] ) );
                    diff_1st_2_file.write( file_list.value( fi ).file_name() + " " );
                    diff_1st_2_file.write_line( double2string( five_membered_ring_analyser.distances_from_plane_[ five_membered_ring_analyser.sorted_map_[ 4 ] ] - five_membered_ring_analyser.distances_from_plane_[ five_membered_ring_analyser.sorted_map_[ 3 ] ] ) );
                }

                if ( five_membered_ring_analyser.is_planar() )
                    analysis_01 += " flat";
                else if ( five_membered_ring_analyser.is_envelope() )
                {
                    analysis_01 += " envelope " + rings_5[j][five_membered_ring_analyser.unique_envelope_point_1()];
                    std::string neighbour_atom_label( "C19" );
                    size_t index = crystal_structure.atom( neighbour_atom_label + "_0" );
                    Atom atom = crystal_structure.atom( index );
                    Vector3D point = crystal_lattice.fractional_to_orthogonal( atom.position() );
                    FiveMemberedRingAnalyser::GeometryType geometry_type = five_membered_ring_analyser.axial_or_equatorial( point );
                    if ( geometry_type == FiveMemberedRingAnalyser::EQUATORIAL )
                        analysis_01 += " " + neighbour_atom_label + " = equatorial";
                    else if ( geometry_type == FiveMemberedRingAnalyser::AXIAL )
                        analysis_01 += " " + neighbour_atom_label + " = axial";
                }
                else if ( five_membered_ring_analyser.is_double_envelope() )
                    analysis_01 += " double_envelope ( " + rings_5[j][five_membered_ring_analyser.unique_envelope_point_1()] + " " +
                                                             rings_5[j][five_membered_ring_analyser.unique_envelope_point_2()] + " )";
                output_string += analysis_01;
                if ( Zprime_is_two )
                {
                    std::vector< Vector3D > points;
                    for ( size_t k( 0 ); k != 5; ++k )
                    {
                        size_t index = crystal_structure.atom( rings_5[j][k] + "_1" );
                        Atom atom = crystal_structure.atom( index );
                        Vector3D point = crystal_lattice.fractional_to_orthogonal( atom.position() );
                        points.push_back( point );
                    }
                    five_membered_ring_analyser.analyse( points );

                    if ( ! five_membered_ring_analyser.is_planar() )
                    {
                        output_file.write( file_list.value( fi ).file_name() + " " );
                        for ( size_t k( 0 ); k != 5; ++k )
                            output_file.write( double2string( five_membered_ring_analyser.distances_from_plane_[ five_membered_ring_analyser.sorted_map_[ 4-k ] ] ) + " " );
                        output_file.write_line();
                        output_file.write( file_list.value( fi ).file_name() + " " );
                        for ( size_t k( 0 ); k != 5; ++k )
                            output_file.write( double2string( five_membered_ring_analyser.rmsds_from_mean_plane_[ five_membered_ring_analyser.sorted_map_[ 4-k ] ] ) + " " );
                        output_file.write_line();
                        distance_from_plane_file.write( file_list.value( fi ).file_name() + " " );
                        distance_from_plane_file.write_line( double2string( five_membered_ring_analyser.distances_from_plane_[ five_membered_ring_analyser.sorted_map_[ 4 ] ] ) );
                        planarity_file.write( file_list.value( fi ).file_name() + " " );
                        planarity_file.write_line( double2string( five_membered_ring_analyser.rmsds_from_mean_plane_[ five_membered_ring_analyser.sorted_map_[ 4 ] ] ) );
                        diff_1st_2_file.write( file_list.value( fi ).file_name() + " " );
                        diff_1st_2_file.write_line( double2string( five_membered_ring_analyser.distances_from_plane_[ five_membered_ring_analyser.sorted_map_[ 4 ] ] - five_membered_ring_analyser.distances_from_plane_[ five_membered_ring_analyser.sorted_map_[ 3 ] ] ) );
                    }

                    if ( five_membered_ring_analyser.is_planar() )
                        analysis_02 += " flat";
                    else if ( five_membered_ring_analyser.is_envelope() )
                    {
                        analysis_02 += " envelope " + rings_5[j][five_membered_ring_analyser.unique_envelope_point_1()];
                        std::string neighbour_atom_label( "C19" );
                        size_t index = crystal_structure.atom( neighbour_atom_label + "_1" );
                        Atom atom = crystal_structure.atom( index );
                        Vector3D point = crystal_lattice.fractional_to_orthogonal( atom.position() );
                        FiveMemberedRingAnalyser::GeometryType geometry_type = five_membered_ring_analyser.axial_or_equatorial( point );
                        if ( geometry_type == FiveMemberedRingAnalyser::EQUATORIAL )
                            analysis_02 += " " + neighbour_atom_label + " = equatorial";
                        else if ( geometry_type == FiveMemberedRingAnalyser::AXIAL )
                            analysis_02 += " " + neighbour_atom_label + " = axial";
                    }
                    else if ( five_membered_ring_analyser.is_double_envelope() )
                        analysis_02 += " double_envelope ( " + rings_5[j][five_membered_ring_analyser.unique_envelope_point_1()] + " " +
                                                                 rings_5[j][five_membered_ring_analyser.unique_envelope_point_2()] + " )";
                    if ( analysis_01 == analysis_02 )
                        output_string += " x2";
                    else
                        output_string += analysis_02;
                }
                else // Z' == 1
                    output_string += "   ";
            }
            analysis_01.clear();
            analysis_02.clear();

            SixMemberedRingAnalyser six_membered_ring_analyser;
            for ( size_t j( 0 ); j != rings_6.size(); ++j )
            {
                std::vector< Vector3D > points;
                for ( size_t k( 0 ); k != rings_6[j].size(); ++k )
                {
                    size_t index = crystal_structure.atom( rings_6[j][k] + "_0" );
                    Atom atom = crystal_structure.atom( index );
                    Vector3D point = crystal_lattice.fractional_to_orthogonal( atom.position() );
                    points.push_back( point );
                }
                six_membered_ring_analyser.analyse( points );
                if ( six_membered_ring_analyser.is_planar() )
                    analysis_01 += " flat";
                else if ( six_membered_ring_analyser.is_chair() )
                {
                    analysis_01 += " chair";
                    std::string neighbour_atom_label( "C13" );
                    size_t index = crystal_structure.atom( neighbour_atom_label + "_0" );
                    Atom atom = crystal_structure.atom( index );
                    Vector3D point = crystal_lattice.fractional_to_orthogonal( atom.position() );
                    SixMemberedRingAnalyser::GeometryType geometry_type = six_membered_ring_analyser.axial_or_equatorial( point );
                    if ( geometry_type == SixMemberedRingAnalyser::EQUATORIAL )
                        analysis_01 += " " + neighbour_atom_label + " = equatorial";
                    else if ( geometry_type == SixMemberedRingAnalyser::AXIAL )
                        analysis_01 += " " + neighbour_atom_label + " = axial";
                }
                // The difference between boat and twisted boat is a gliding scale.
                else if ( six_membered_ring_analyser.is_boat() )
                    analysis_01 += " boat";
                else if ( six_membered_ring_analyser.is_twisted_boat() )
                    analysis_01 += " twisted boat";

                output_string += analysis_01;
                if ( Zprime_is_two )
                {
                    std::vector< Vector3D > points;
                    for ( size_t k( 0 ); k != rings_6[j].size(); ++k )
                    {
                        size_t index = crystal_structure.atom( rings_6[j][k] + "_1" );
                        Atom atom = crystal_structure.atom( index );
                        Vector3D point = crystal_lattice.fractional_to_orthogonal( atom.position() );
                        points.push_back( point );
                    }
                    six_membered_ring_analyser.analyse( points );
                    if ( six_membered_ring_analyser.is_planar() )
                        analysis_02 += " flat";
                    else if ( six_membered_ring_analyser.is_chair() )
                    {
                        analysis_02 += " chair";
                        std::string neighbour_atom_label( "C13" );
                        size_t index = crystal_structure.atom( neighbour_atom_label + "_1" );
                        Atom atom = crystal_structure.atom( index );
                        Vector3D point = crystal_lattice.fractional_to_orthogonal( atom.position() );
                        SixMemberedRingAnalyser::GeometryType geometry_type = six_membered_ring_analyser.axial_or_equatorial( point );
                        if ( geometry_type == SixMemberedRingAnalyser::EQUATORIAL )
                            analysis_02 += " " + neighbour_atom_label + " = equatorial";
                        else if ( geometry_type == SixMemberedRingAnalyser::AXIAL )
                            analysis_02 += " " + neighbour_atom_label + " = axial";
                    }
                    // The difference between boat and twisted boat is a gliding scale.
                    else if ( six_membered_ring_analyser.is_boat() )
                        analysis_02 += " boat";
                    else if ( six_membered_ring_analyser.is_twisted_boat() )
                        analysis_02 += " twisted boat";
                    if ( analysis_01 == analysis_02 )
                        output_string += " x2";
                    else
                        output_string += analysis_02;
                }
                else // Z' == 1
                    output_string += "   ";

            }
            text_file_writer.write_line( output_string );
        }
    MACRO_END_GAME

    try // Rebin powder pattern
    {
        if ( argc != 2 )
            throw std::runtime_error( "Please give the name of a .xye file." );
        FileName input_file_name( argv[ 1 ] );
        PowderPattern powder_pattern;
        powder_pattern.read_xye( input_file_name );
        size_t bin_size = 10;
        powder_pattern.rebin( bin_size );
        powder_pattern.save_xye( append_to_file_name( input_file_name, "_rebin_" + size_t2string( bin_size ) ), true );
    MACRO_END_GAME

    try // Remove H atoms from a set of cif files.
    {
        MACRO_ONE_FILELISTNAME_AS_ARGUMENT
        for ( size_t i( 0 ); i != file_list.size(); ++i )
            remove_hydrogen_atoms( file_list.value( i ) );
    MACRO_END_GAME

    try // Make atom labels unique.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
    //    crystal_structure.set_space_group( SpaceGroup() );
        crystal_structure.make_atom_labels_unique();
        crystal_structure.save_cif( append_to_file_name( input_file_name, "_ual" ) );
    MACRO_END_GAME

    try // Check space group
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        check_if_closed( crystal_structure.space_group().symmetry_operators() );
    MACRO_END_GAME

    try // Convert powder pattern in .brml format to .xye.
    {
        // First, the user must change .brml to .zip and extract all files
        if ( argc != 2 )
            throw std::runtime_error( "Please give the name and path of the RawData0.xml file." );
        FileName input_file_name( argv[ 1 ] );
        PowderPattern powder_pattern;
        powder_pattern.read_brml( input_file_name );
        powder_pattern.save_xye( replace_extension( input_file_name, "xye" ), true );
    MACRO_END_GAME

    try // Read Norman's file.
    {
        TextFileReader text_file_reader( FileName( "C:\\Data_Win\\iwrgbihwrbghwirb.txt" ) );
        std::string line;
        size_t iLine( 0 );
        RunningAverageAndESD<double> running_average;
        while ( text_file_reader.get_next_line( line ) )
        {
            ++iLine;
            if ( ( iLine % 8 ) == 0 )
            {
                double RMSCD = string2double( line );
                std::cout << RMSCD << std::endl;
                if ( RMSCD > 0.19 )
                    RMSCD = 0.042;
                running_average.add_value( RMSCD );
            }
        }
        std::cout << running_average.average() << std::endl;
        std::cout << running_average.estimated_standard_deviation() << std::endl;
    MACRO_END_GAME

    try // Make atom labels unique.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
    //    crystal_structure.set_space_group( SpaceGroup() );
        crystal_structure.make_atom_labels_unique();
        crystal_structure.save_cif( append_to_file_name( input_file_name, "_ual" ) );
    MACRO_END_GAME

    try // Generate a powder diffraction pattern from SHELX .hkl file.
    {
        ReflectionList reflection_list;
        reflection_list.read_hkl( FileName( "C:\\Data_Win\\ContractResearch\\GC\\Molecule05\\Experimental_data\\GP.hkl" ) );
        std::cout << ".hkl file has been read." << std::endl;
        // We need a CrystalLattice
        CrystalLattice crystal_lattice(  4.843,
                                         6.118,
                                        29.72,
                                        Angle::from_degrees( 87.384 ),
                                        Angle::from_degrees( 85.734 ),
                                        Angle::from_degrees( 70.335 ) );
        SpaceGroup space_group;
        PointGroup laue_class = space_group.laue_class();
        ReflectionList reflection_list_final;
        std::vector< bool > done( reflection_list.size(), false );
        for ( size_t i(0); i != reflection_list.size(); ++i )
        {
            if ( done[i] )
                continue;
            MillerIndices miller_indices = reflection_list.miller_indices( i );
            // Calculate equivalent reflections
            std::set< MillerIndices > equivalent_reflections;
            equivalent_reflections.insert( miller_indices );
            for ( size_t j( 0 ); j != laue_class.nsymmetry_operators(); ++j )
                equivalent_reflections.insert( miller_indices * laue_class.symmetry_operator( j ) );
            double F_squared = 0.0;
            size_t number_of_equivalent_reflections_that_have_actually_been_measured( 0 );
            for ( std::set< MillerIndices >::const_iterator it = equivalent_reflections.begin(); it != equivalent_reflections.end(); ++it )
            {
                size_t index = reflection_list.index( *it );
                if ( index == reflection_list.size() )
                {
                    std::cout << "Oops..." << std::endl;
                }
                else
                {
                    ++number_of_equivalent_reflections_that_have_actually_been_measured;
                    done[index] = true;
                    F_squared += reflection_list.F_squared( index );
                }
            }
            F_squared /= number_of_equivalent_reflections_that_have_actually_been_measured;
            if ( F_squared < 0.0 )
                F_squared = 0.0;
            Vector3D H = miller_indices.h() * crystal_lattice.a_star_vector() +
                         miller_indices.k() * crystal_lattice.b_star_vector() +
                         miller_indices.l() * crystal_lattice.c_star_vector();
            double d = 1.0 / ( H.length() );
            reflection_list_final.push_back( miller_indices, F_squared, d, equivalent_reflections.size() );
        }
        CrystalStructure crystal_structure;
        crystal_structure.set_space_group( space_group );
        PowderPatternCalculator powder_pattern_calculator( crystal_structure );
        powder_pattern_calculator.set_FWHM( 0.05 );
        powder_pattern_calculator.set_two_theta_end( Angle::from_degrees( 40.0 ) );
        powder_pattern_calculator.set_two_theta_step( Angle::from_degrees( 0.01 ) );
        PowderPattern powder_pattern;
        powder_pattern_calculator.calculate( reflection_list_final, powder_pattern );
        powder_pattern.save_xye( FileName( "C:\\Data_Win\\ContractResearch\\GC\\Molecule05\\Experimental_data\\GP_new.xye" ), true );
    MACRO_END_GAME

    try // Subtract two powder patterns.
    {
        FileName input_file_name_1( "C:\\Data_Win\\ContractResearch\\GC\\Molecule06\\GP_profile.xye" );
        PowderPattern powder_pattern_1;
        powder_pattern_1.read_xye( input_file_name_1 );
        FileName input_file_name_2( "C:\\Data_Win\\ContractResearch\\GC\\Molecule06\\GP_BKGR.xye" );
        PowderPattern powder_pattern_2;
        powder_pattern_2.read_xye( input_file_name_2 );
        powder_pattern_1 -= powder_pattern_2;

        powder_pattern_1.correct_zero_point_error( Angle::from_degrees( 0.14068 ) );
        powder_pattern_1.save_xye( append_to_file_name( input_file_name_1, "_zp" ), true );
    MACRO_END_GAME

    try // After DASH 6 finish .inp file.
    {
        if ( argc != 2 )
            throw std::runtime_error( "Please give the name of a .inp file that needs to be completed as argument." );
        FileName input_file_name( argv[ 1 ] );
        finish_inp( input_file_name );
    MACRO_END_GAME

    try // Convert powder pattern in .MDI format to .xye.
    {
        if ( argc != 2 )
            throw std::runtime_error( "Please give the name of a .mdi file." );
        FileName input_file_name( argv[ 1 ] );
        PowderPattern powder_pattern;
        powder_pattern.read_mdi( input_file_name );
        powder_pattern.save_xye( replace_extension( input_file_name, "xye" ), true );
    MACRO_END_GAME

    try // Calculate dipole moment for FileList.txt.
    {
        MACRO_ONE_FILELISTNAME_AS_ARGUMENT
        TextFileWriter text_file_writer( FileName( file_list_file_name.directory(), "dipole_results", "txt" ) );
        for ( size_t i( 0 ); i != file_list.size(); ++i )
        {
            CrystalStructure crystal_structure;
            std::cout << "Now reading cif... " + file_list.value( i ).full_name() << std::endl;
            read_cif( file_list.value( i ), crystal_structure );
            double dipole_moment_molecule = crystal_structure.dipole_moment();
            crystal_structure.apply_space_group_symmetry();
            double dipole_moment_crystal = crystal_structure.dipole_moment();
            if ( dipole_moment_crystal < 1E-10 )
                dipole_moment_crystal = 0.0;
            text_file_writer.write_line( FileName( "", file_list.value( i ).file_name(), file_list.value( i ).extension() ).full_name() + " " + double2string( dipole_moment_molecule ) + " " + double2string( dipole_moment_crystal ) );
        }
    MACRO_END_GAME

    try // Calculate dipole moment.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
//        SpaceGroup space_group;
//        space_group.add_inversion_at_origin();
//        crystal_structure.set_space_group( space_group );
        std::cout << crystal_structure.dipole_moment() << std::endl;
        crystal_structure.apply_space_group_symmetry();
        std::cout << crystal_structure.dipole_moment() << std::endl;
    MACRO_END_GAME

    try // Generate R input file for Pawley or Loopstra-Rietveld plot.
        {
        if ( argc != 3 )
        {
            std::cout << "The following files are required:" << std::endl;
            std::cout << "<directory>\\<base name>.xye" << std::endl;
            std::cout << "<directory>\\<base name>_profile.txt" << std::endl;
            std::cout << "<directory>\\<base name>_tickmarks.txt" << std::endl;
            throw std::runtime_error( "Please supply a <directory> and a <base name> as command-line arguments." );
        }
        std::string directory( append_backslash( argv[ 1 ] ) );
        std::string base_name( argv[ 2 ] );
        GeneratePowderCIF generate_powder_cif( directory, base_name );
    //enum ZoomPolicy { ALWAYS_ZOOM, NEVER_ZOOM, ZOOM_OVER_40 };
        generate_powder_cif.generate_R_input_file( GeneratePowderCIF::ZOOM_OVER_40 );
    MACRO_END_GAME

    try // Transform unit cell, NOT the atomic coordinates.
    {
        CrystalLattice crystal_lattice( 19.94392,
                                        14.69558,
                                         8.50626,
                                        Angle::from_degrees( 93.506 ),
                                        Angle::from_degrees( 91.852 ),
                                        Angle::from_degrees( 101.390 ) );

        Matrix3D transformation_matrix(  0.0,  1.0,  0.0,
                                         0.0,  0.0,  1.0,
                                         1.0,  0.0,  0.0 );
        crystal_lattice.transform( transformation_matrix );
     //   transformation_matrix.invert();
      //  transformation_matrix.show();

        std::cout << crystal_lattice.a() << " " <<
                     crystal_lattice.b() << " " <<
                     crystal_lattice.c() << " " <<
                     crystal_lattice.alpha() << " " <<
                     crystal_lattice.beta() << " " <<
                     crystal_lattice.gamma() << std::endl;
        crystal_lattice.print();
    MACRO_END_GAME

   // Move group of atoms.
    try
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        std::vector< std::string > atoms_to_move;
//        atoms_to_move.push_back( "O3" );
//        atoms_to_move.push_back( "P1" );
//        atoms_to_move.push_back( "O2" );
//        atoms_to_move.push_back( "O4" );
//        atoms_to_move.push_back( "H1" );

        // Vector to move along
        // NB direction
        // NB fractional vs orthogonal coordinates

        Vector3D difference_frac;
        if ( (false) )
        {
            // Move along C3-O1 in the direction of C3.
            std::string from_atom_label( "O1" );
            std::string to_atom_label( "C3" );

            Vector3D from_atom_frac = crystal_structure.atom( crystal_structure.find_label( from_atom_label ) ).position();
            Vector3D to_atom_frac = crystal_structure.atom( crystal_structure.find_label( to_atom_label ) ).position();
            // Find shortest disance between the two taking periodicity and space-group symmetry into account.
            double distance;
            crystal_structure.shortest_distance( from_atom_frac, to_atom_frac, distance, difference_frac );

            // We need to be able to set the distance, so we must work in Cartesian coordinates.
            Vector3D difference_cart = crystal_structure.crystal_lattice().fractional_to_orthogonal( difference_frac );
            double target_distance( 0.25 );
            difference_cart *= target_distance / distance;
            difference_frac = crystal_structure.crystal_lattice().orthogonal_to_fractional( difference_cart );
        }
        else
        {
            difference_frac = Vector3D( 0.0, -0.047, 0.0 );
        }

        if ( atoms_to_move.empty() )
        {
            for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
            {
                Atom new_atom( crystal_structure.atom( i ) );
                Vector3D new_position = new_atom.position() + difference_frac;
                new_atom.set_position( new_position );
                crystal_structure.set_atom( i, new_atom );
            }
        }
        else
        {
            for ( size_t i( 0 ); i != atoms_to_move.size(); ++i )
            {
                Atom new_atom( crystal_structure.atom( crystal_structure.find_label( atoms_to_move[i] ) ) );
                Vector3D new_position = new_atom.position() + difference_frac;
                new_atom.set_position( new_position );
                crystal_structure.set_atom( crystal_structure.find_label( atoms_to_move[i] ), new_atom );
            }
        }

        crystal_structure.save_cif( append_to_file_name( input_file_name, "_atoms_moved" ) );
    MACRO_END_GAME

    try // Add powder patterns, assumes VCT. Only lab data, for now.
    {
        std::vector< PowderPattern > powder_patterns;
        std::vector< double > noscp2ts;

        {
        PowderPattern powder_pattern;
        powder_pattern.read_xye( FileName( "C:\\Data_Win\\KOFHOC\\VCT\\NewSample\\range_01.xye" ) );
        powder_patterns.push_back( powder_pattern );
        noscp2ts.push_back( 1.0 );
        }
        {
        PowderPattern powder_pattern;
        powder_pattern.read_xye( FileName( "C:\\Data_Win\\KOFHOC\\VCT\\NewSample\\range_02.xye" ) );
        powder_patterns.push_back( powder_pattern );
        noscp2ts.push_back( 1.0 );
        }
        {
        PowderPattern powder_pattern;
        powder_pattern.read_xye( FileName( "C:\\Data_Win\\KOFHOC\\VCT\\NewSample\\range_03.xye" ) );
        powder_patterns.push_back( powder_pattern );
        noscp2ts.push_back( 2.0 );
        }
        {
        PowderPattern powder_pattern;
        powder_pattern.read_xye( FileName( "C:\\Data_Win\\KOFHOC\\VCT\\NewSample\\range_04.xye" ) );
        powder_patterns.push_back( powder_pattern );
        noscp2ts.push_back( 4.0 );
        }
        PowderPattern sum = add_powder_patterns( powder_patterns, noscp2ts );
        sum.save_xye( FileName( "C:\\Data_Win\\KOFHOC\\VCT\\NewSample\\VCT_sum.xye" ), true );
    MACRO_END_GAME

    try // Write FileList.txt file.
    {
        TextFileWriter text_file_writer( FileName( "C:\\Data_Win\\\\FileList.txt" ) );
        for ( size_t i( 0 ); i != 7252; ++i )
        {
            text_file_writer.write_line( "structure_" + size_t2string( i+1, 6, '0' ) + ".cif" );
        }
    MACRO_END_GAME

    try // Calculate ADPs from a set of frames (as cif files).
    {
        MACRO_ONE_FILELISTNAME_AS_ARGUMENT
        SpaceGroup space_group; // = SpaceGroup::P21c();
//        SpaceGroup space_group = SpaceGroup::P21c();
//        space_group.add_inversion_at_origin(); // Space group is now P-1
        // Stores results as .cif
//        Matrix3D transformation_matrix(  0.25,  0.0,  0.25,
//                                         0.0,  1.0,  0.0,
//                                        -1.0,  0.0,  0.0 );
        Matrix3D transformation_matrix(  1.0,  0.0,  0.0,
                                         0.0,  1.0,  0.0,
                                         0.0,  0.0,  1.0 );
        AnalyseTrajectory analyse_trajectory( file_list, 6, 12, 6, space_group, transformation_matrix );
  //      analyse_trajectory.save_centres_of_mass();
    MACRO_END_GAME

    try // Calculate all delta_AB for ADPs for all pairs of atoms to check if rigid-body approximation for TLS is valid.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        if ( crystal_structure.natoms() == 0 )
            return 0;
        std::vector< std::string > labels;
        std::vector< double > deltas;
        double max( 0.0 );
        for ( size_t i( 0 ); i != crystal_structure.natoms()-1; ++i )
        {
            if ( crystal_structure.atom( i ).ADPs_type() != Atom::ANISOTROPIC )
                continue;
            for ( size_t j( i+1 ); j != crystal_structure.natoms(); ++j )
            {
                if ( crystal_structure.atom( j ).ADPs_type() != Atom::ANISOTROPIC )
                    continue;
                Vector3D v = crystal_structure.atom( i ).position() - crystal_structure.atom( j ).position();
                double mu_2_v_i = crystal_structure.atom( i ).anisotropic_displacement_parameters().average_displacement_squared( v, crystal_structure.crystal_lattice() );
                double mu_2_v_j = crystal_structure.atom( j ).anisotropic_displacement_parameters().average_displacement_squared( v, crystal_structure.crystal_lattice() );
                labels.push_back( crystal_structure.atom( i ).label() + " " + crystal_structure.atom( j ).label() );
                deltas.push_back( std::abs( mu_2_v_i - mu_2_v_j ) );
                if ( std::abs( mu_2_v_i - mu_2_v_j ) > max )
                    max = std::abs( mu_2_v_i - mu_2_v_j );
//                std::cout << crystal_structure.atom( i ).label() << " " << crystal_structure.atom( j ).label() << " " << std::abs( mu_2_v_i - mu_2_v_j ) << std::endl;
            }
        }
        for ( size_t i( 0 ); i != deltas.size(); ++i )
            std::cout << labels[i] << " " << ASCII_histogram( 0, max, deltas[i], 25, ' ' ) << deltas[i] << std::endl;

    MACRO_END_GAME

    try // R. T. Downs tests.
    {
        CrystalLattice crystal_lattice( 14.333418,
                                         5.562235,
                                         8.530211,
                                        Angle::from_degrees( 90.00000 ),
                                        Angle::from_degrees( 92.02191 ),
                                        Angle::from_degrees( 90.00000 ) );

//    prm C01_x   0.01598
//    prm C01_y   0.09382
//    prm C01_z   0.05908
//    prm C02_x  -0.02925
//    prm C02_y   0.06049
//    prm C02_z   0.21568
//
//    site C1  x =  C01_x; y =  C01_y; z =  C01_z; occ C  0.5 ADPs_Keep_PD u11=u11C1; u22=u22C1; u33=u33C1; u12=u12C1; u13=u13C1; u23=u23C1;
//    site C2  x =  C02_x; y =  C02_y; z =  C02_z; occ C  0.5 ADPs_Keep_PD u11=u11C2; u22=u22C2; u33=u33C2; u12=u12C2; u13=u13C2; u23=u23C2;

        Vector3D v_C1(  0.01598, 0.09382, 0.05908 );
        Vector3D v_C2( -0.02925, 0.06049, 0.21568 );

//    prm ru11C1 = l22*rzC1^2 + l33*ryC1^2 - 2*ryC1*rzC1*l23 - 2*ryC1*s31 + 2*rzC1*s21 + t11; :  0.07362
//    prm ru22C1 = l11*rzC1^2 + l33*rxC1^2 - 2*rxC1*rzC1*l13 - 2*rzC1*s12 + 2*rxC1*s32 + t22; :  0.02469
//    prm ru33C1 = l11*ryC1^2 + l22*rxC1^2 - 2*rxC1*ryC1*l12 - 2*rxC1*s23 + 2*ryC1*s13 + t33; :  0.05492
//    prm ru12C1 = -rxC1*ryC1*l33 + rxC1*rzC1*l23 + ryC1*rzC1*l13 - rzC1^2*l12 - rzC1*s11 + rzC1*s22 + rxC1*s31 - ryC1*s32 + t12; :  0.02061
//    prm ru13C1 = -rxC1*rzC1*l22 + rxC1*ryC1*l23 - ryC1^2*l13 + ryC1*rzC1*l12 + ryC1*s11 - ryC1*s33 + rzC1*s23 - rxC1*s21 + t13; : -0.00637
//    prm ru23C1 = -ryC1*rzC1*l11 - rxC1^2*l23 + rxC1*ryC1*l13 + rxC1*rzC1*l12 - rxC1*s22 + rxC1*s33 + ryC1*s12 - rzC1*s13 + t23; :  0.00335
//
//    prm ru11C2 = l22*rzC2^2 + l33*ryC2^2 - 2*ryC2*rzC2*l23 - 2*ryC2*s31 + 2*rzC2*s21 + t11; :  0.11731
//    prm ru22C2 = l11*rzC2^2 + l33*rxC2^2 - 2*rxC2*rzC2*l13 - 2*rzC2*s12 + 2*rxC2*s32 + t22; :  0.06507
//    prm ru33C2 = l11*ryC2^2 + l22*rxC2^2 - 2*rxC2*ryC2*l12 - 2*rxC2*s23 + 2*ryC2*s13 + t33; :  0.05481
//    prm ru12C2 = -rxC2*ryC2*l33 + rxC2*rzC2*l23 + ryC2*rzC2*l13 - rzC2^2*l12 - rzC2*s11 + rzC2*s22 + rxC2*s31 - ryC2*s32 + t12; :  0.03667
//    prm ru13C2 = -rxC2*rzC2*l22 + rxC2*ryC2*l23 - ryC2^2*l13 + ryC2*rzC2*l12 + ryC2*s11 - ryC2*s33 + rzC2*s23 - rxC2*s21 + t13; :  0.00828
//    prm ru23C2 = -ryC2*rzC2*l11 - rxC2^2*l23 + rxC2*ryC2*l13 + rxC2*rzC2*l12 - rxC2*s22 + rxC2*s33 + ryC2*s12 - rzC2*s13 + t23; :  0.00183

    // Same order as u11 etc. in .cif file
    //SymmetricMatrix3D( const double a00, const double a11, const double a22,
    //                   const double a01, const double a02, const double a12 );

        // Expects U_cart
        AnisotropicDisplacementParameters ADPs_C1( SymmetricMatrix3D( 0.07362, 0.02469, 0.05492, 0.02061, -0.00637, 0.00335 ) );
        AnisotropicDisplacementParameters ADPs_C2( SymmetricMatrix3D( 0.11731, 0.06507, 0.05481, 0.03667,  0.00828, 0.00183 ) );

        Vector3D v = v_C2 - v_C1;

        double mu2_v_1 = ( v * ( transpose( crystal_lattice.Downs_G() ) * ( 2.0 * square( CONSTANT_PI ) * ADPs_C1.U_star( crystal_lattice ) ) * crystal_lattice.Downs_G() ) * v ) / ( 2.0 * square( CONSTANT_PI ) * v * crystal_lattice.Downs_G() * v );
        double mu2_v_2 = ( v * ( transpose( crystal_lattice.Downs_G() ) * ( 2.0 * square( CONSTANT_PI ) * ADPs_C2.U_star( crystal_lattice ) ) * crystal_lattice.Downs_G() ) * v ) / ( 2.0 * square( CONSTANT_PI ) * v * crystal_lattice.Downs_G() * v );
        std::cout << "delta_AB = " << mu2_v_2 - mu2_v_1 << std::endl;

        double delta_AB = v * ( crystal_lattice.Downs_D() * ( ADPs_C2.U_cart() - ADPs_C1.U_cart() ) * crystal_lattice.Downs_D() ) * v;
        std::cout << "delta_AB = " << delta_AB << std::endl;

    MACRO_END_GAME

    // Save all disorder combinations.
    try
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        std::vector< DisorderGroup > disorder_groups;
        DisorderGroup disorder_group;
        disorder_group.major_occupancy_labels_.push_back( "N202" );

        disorder_group.minor_occupancy_labels_.push_back( "N206" );
        disorder_groups.push_back( disorder_group );
        disorder_group.major_occupancy_labels_.clear();
        disorder_group.minor_occupancy_labels_.clear();
        disorder_group.major_occupancy_labels_.push_back( "C221" );
        disorder_group.major_occupancy_labels_.push_back( "H05A" );
        disorder_group.major_occupancy_labels_.push_back( "H05B" );

        disorder_group.minor_occupancy_labels_.push_back( "C223" );
        disorder_group.minor_occupancy_labels_.push_back( "H99A" );
        disorder_group.minor_occupancy_labels_.push_back( "H99B" );
        disorder_groups.push_back( disorder_group );
        disorder_group.major_occupancy_labels_.clear();
        disorder_group.minor_occupancy_labels_.clear();
        disorder_group.major_occupancy_labels_.push_back( "C121" );
        disorder_group.major_occupancy_labels_.push_back( "H14A" );
        disorder_group.major_occupancy_labels_.push_back( "H14B" );
        disorder_group.major_occupancy_labels_.push_back( "F101" );
        disorder_group.major_occupancy_labels_.push_back( "F102" );
        disorder_group.major_occupancy_labels_.push_back( "F103" );

        disorder_group.minor_occupancy_labels_.push_back( "C134" );
        disorder_group.minor_occupancy_labels_.push_back( "H98A" );
        disorder_group.minor_occupancy_labels_.push_back( "H98B" );
        disorder_group.minor_occupancy_labels_.push_back( "F1" );
        disorder_group.minor_occupancy_labels_.push_back( "F2" );
        disorder_group.minor_occupancy_labels_.push_back( "F3" );
        disorder_groups.push_back( disorder_group );
        disorder_group.major_occupancy_labels_.clear();
        disorder_group.minor_occupancy_labels_.clear();
        disorder_group.major_occupancy_labels_.push_back( "N105" );
        disorder_group.major_occupancy_labels_.push_back( "C118" );
        disorder_group.major_occupancy_labels_.push_back( "H9" );
        disorder_group.major_occupancy_labels_.push_back( "C120" );
        disorder_group.major_occupancy_labels_.push_back( "H3" );
        disorder_group.major_occupancy_labels_.push_back( "C117" );
        disorder_group.major_occupancy_labels_.push_back( "H8" );

        disorder_group.minor_occupancy_labels_.push_back( "N106" );
        disorder_group.minor_occupancy_labels_.push_back( "C132" );
        disorder_group.minor_occupancy_labels_.push_back( "H80" );
        disorder_group.minor_occupancy_labels_.push_back( "C133" );
        disorder_group.minor_occupancy_labels_.push_back( "H82" );
        disorder_group.minor_occupancy_labels_.push_back( "C131" );
        disorder_group.minor_occupancy_labels_.push_back( "H81" );
        disorder_groups.push_back( disorder_group );

        size_t npermutations = 1 << disorder_groups.size();
        std::cout << "npermutations = " << npermutations<< std::endl;
        for ( size_t i( 0 ); i != npermutations; ++i )
        {
            for ( size_t k( 0 ); k != crystal_structure.natoms(); ++k )
                crystal_structure.set_suppressed( k, false );
            for ( size_t j( 0 ); j != disorder_groups.size(); ++j )
            {
                if ( ( i & ( 1 << j ) ) == 0 )
                {
                    // Major occupancy
                    for ( size_t l( 0 ); l != disorder_groups[j].minor_occupancy_labels_.size(); ++l )
                    {
                        size_t index = crystal_structure.find_label( disorder_groups[j].minor_occupancy_labels_[l] );
                        if ( index == crystal_structure.natoms() )
                            throw std::runtime_error( "label >" + disorder_groups[j].minor_occupancy_labels_[l] + "< not found." );
                      //  std::cout << "disorder_groups[j].minor_occupancy_labels_[l] 0 = " << disorder_groups[j].minor_occupancy_labels_[l] << std::endl;
                      //  std::cout << "index 0 = " << index << std::endl;
                        crystal_structure.set_suppressed( index, true );
                    }
                }
                else
                {
                    // Minor occupancy
                    for ( size_t l( 0 ); l != disorder_groups[j].major_occupancy_labels_.size(); ++l )
                    {
                        size_t index = crystal_structure.find_label( disorder_groups[j].major_occupancy_labels_[l] );
                        if ( index == crystal_structure.natoms() )
                            throw std::runtime_error( "label >" + disorder_groups[j].minor_occupancy_labels_[l] + "< not found." );
                      //  std::cout << "disorder_groups[j].minor_occupancy_labels_[l] 1 = " << disorder_groups[j].minor_occupancy_labels_[l] << std::endl;
                      //  std::cout << "index 1 = " << index << std::endl;
                        crystal_structure.set_suppressed( index, true );
                    }
                }
            }
            crystal_structure.save_cif( append_to_file_name( input_file_name, size_t2string( i, 2, '0' ) ) );
        }
    MACRO_END_GAME

    try // Convert Materials Studio .xsd to .cif.
    {
        if ( argc != 2 ) \
            throw std::runtime_error( "Please give the name of a .xsd file." ); \
        FileName input_file_name( argv[ 1 ] ); \
        CrystalStructure crystal_structure; \
        std::cout << "Now reading xsd... " + input_file_name.full_name() << std::endl; \
        read_xsd( input_file_name, crystal_structure );
        crystal_structure.save_cif( replace_extension( input_file_name, "cif" ) );
    MACRO_END_GAME

    try // Make atom labels unique for a list of cif files.
    {
        MACRO_ONE_FILELISTNAME_AS_ARGUMENT
        for ( size_t i( 0 ); i != file_list.size(); ++i )
        {
            CrystalStructure crystal_structure;
            std::cout << "Now reading cif... " + file_list.value(i).full_name() << std::endl;
            read_cif( file_list.value(i), crystal_structure );
            crystal_structure.make_atom_labels_unique();
            crystal_structure.save_cif( append_to_file_name( file_list.value(i), "_ual" ) );
        }
    MACRO_END_GAME

    try // Find unit-cell transformation with a space-group setting as the target.
    {
        if ( argc < 2 )
            throw std::runtime_error( "Please give the name of a .cif file." );
        FileName input_file_name( argv[ 1 ] );
        CrystalStructure crystal_structure;
        std::cout << "Now reading cif... " + input_file_name.full_name() << std::endl;
        read_cif( input_file_name, crystal_structure );
        CrystalLattice old_crystal_lattice = crystal_structure.crystal_lattice();
        SpaceGroup old_space_group = crystal_structure.space_group();
        SpaceGroup target_space_group;
        if ( argc == 3 )
        {
            FileName file_name_2( argv[ 2 ] );
            CrystalStructure crystal_structure_2;
            std::cout << "Now reading cif... " + file_name_2.full_name() << std::endl;
            read_cif( file_name_2, crystal_structure_2 );
            target_space_group = crystal_structure_2.space_group();
        }
        else
        {
            std::vector< SymmetryOperator > symmetry_operators;
            symmetry_operators.push_back( SymmetryOperator( "x,y,z" ) );
            symmetry_operators.push_back( SymmetryOperator( "-x,y,-z" ) );
            symmetry_operators.push_back( SymmetryOperator( "1/2+x,1/2+y,z" ) );
            symmetry_operators.push_back( SymmetryOperator( "1/2-x,1/2+y,-z" ) );
            target_space_group = SpaceGroup( symmetry_operators, "C2" );
        }
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
                    Matrix3D transformation_matrix( i1, i2, i3, j1, j2, j3, k1, k2, k3 );
                    if ( ! nearly_equal( transformation_matrix.determinant(), 1.0 ) )
                        continue;
                    // Make a copy
                    SpaceGroup new_space_group( old_space_group );
                    // Transform
                    Matrix3D transformation_matrix_inverse_transpose( transformation_matrix );
                    transformation_matrix_inverse_transpose.invert();
                    transformation_matrix_inverse_transpose.transpose();
                    new_space_group.apply_similarity_transformation( SymmetryOperator( transformation_matrix_inverse_transpose, Vector3D() ) );
                    if ( same_symmetry_operators( new_space_group, target_space_group ) )
                    {
                        transformation_matrix.show();
                        CrystalLattice new_lattice( old_crystal_lattice );
                        new_lattice.transform( transformation_matrix );
                        new_lattice.print();
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

    try // Complete cif by merging .inp, _Hmi.cif, .cif and lists with bond lengths and valence angles.
    {
        if ( argc != 3 )
        {
            std::cout << "The following files are required:" << std::endl;
            std::cout << "<directory>\\<base name>.cif" << std::endl;
            std::cout << "<directory>\\<base name>.inp" << std::endl;
            std::cout << "<directory>\\<base name>.xye" << std::endl;
            std::cout << "<directory>\\<base name>_profile.txt" << std::endl;
            std::cout << "<directory>\\<base name>_tickmarks.txt" << std::endl;
            std::cout << "<directory>\\<base name>_additional_info.txt" << std::endl;
            std::cout << "<directory>\\<base name>-Bonds.tsv" << std::endl;
            std::cout << "<directory>\\<base name>-AllAngles.tsv" << std::endl;
            std::cout << "The following files are optional:" << std::endl;
            std::cout << "<directory>\\<base name>_Hmi.cif (from GRACE, Mercury or Materials Studio)" << std::endl;
            std::cout << "<directory>\\<base name>_relabel.txt" << std::endl;
            throw std::runtime_error( "Please supply a <directory> and a <base name> as command-line arguments." );
        }
        std::string directory( append_backslash( argv[ 1 ] ) );
        std::string base_name( argv[ 2 ] );
        GeneratePowderCIF generate_powder_cif( directory, base_name );
        generate_powder_cif.check_new_labels();
        generate_powder_cif.generate();
        generate_powder_cif.check_new_labels();
    MACRO_END_GAME

    try // Drunkard's walk.
    {
        int size = 100;
        RunningAverageAndESD<double> static_average_nmoves;
        RunningAverageAndESD<double> moving_average_nmoves;

        RandomNumberGenerator_integer seed_generator;
        for ( size_t i( 0 ); i != 1000; ++i )
        {
            if ( i < 10 )
                std::cout << "Starting simulation number " << i << std::endl;
            else if ( ( i < 100 ) && ( ( i % 10 ) == 0 ) )
                std::cout << "Starting simulation number " << i << std::endl;
            else if ( ( ( i % 100 ) == 0 ) )
                std::cout << "Starting simulation number " << i << std::endl;
            DrunkardsWalk drunkards_walk_1;
            DrunkardsWalk drunkards_walk_2;
            bool use_random_starting_positions( true );
            if ( use_random_starting_positions )
            {
                BagOfNumbers bag_of_numbers( size+1, seed_generator.next_number( 1, 1000000 ) );
                drunkards_walk_1 = DrunkardsWalk( GridPoint2D( 0, 0 ), GridPoint2D( size, size ), GridPoint2D( bag_of_numbers.draw_with_replace(), bag_of_numbers.draw_with_replace() ), seed_generator.next_number( 1, 1000000 ) );
                drunkards_walk_2 = DrunkardsWalk( GridPoint2D( 0, 0 ), GridPoint2D( size, size ), GridPoint2D( bag_of_numbers.draw_with_replace(), bag_of_numbers.draw_with_replace() ), seed_generator.next_number( 1, 1000000 ) );
            }
            else
            {
                drunkards_walk_1 = DrunkardsWalk( GridPoint2D( 0, 0 ), GridPoint2D( size, size ), GridPoint2D( 28, 66 ), seed_generator.next_number( 1, 1000000 ) );
                drunkards_walk_2 = DrunkardsWalk( GridPoint2D( 0, 0 ), GridPoint2D( size, size ), GridPoint2D( 96, 37 ), seed_generator.next_number( 1, 1000000 ) );
            }
            GridPoint2D current_position_1 = drunkards_walk_1.current_position();
            GridPoint2D current_position_2 = drunkards_walk_2.current_position();
            GridPoint2D static_starting_position = drunkards_walk_2.current_position();
            bool static_person_found( false );
            size_t static_person_found_nmoves( 0 );
            bool moving_person_found( false );
            size_t moving_person_found_nmoves( 0 );
            size_t nmoves( 0 );

            if ( ( current_position_1.x_ == static_starting_position.x_ ) && ( current_position_1.y_ == static_starting_position.y_ ) )
            {
                static_person_found = true;
                static_person_found_nmoves = nmoves;
            }
            if ( ( current_position_1.x_ == current_position_2.x_ ) && ( current_position_1.y_ == current_position_2.y_ ) )
            {
                moving_person_found = true;
                moving_person_found_nmoves = nmoves;
            }
            ++nmoves;
            current_position_1 = drunkards_walk_1.next_step();
            if ( ( ( current_position_1.x_ + current_position_1.y_ ) % 2 ) != ( ( current_position_2.x_ + current_position_2.y_ ) % 2 ) )
                current_position_2 = drunkards_walk_2.next_step();

            while ( ( ! static_person_found ) || ( ! moving_person_found ) )
            {

                if ( ! static_person_found )
                {
                    if ( ( current_position_1.x_ == static_starting_position.x_ ) && ( current_position_1.y_ == static_starting_position.y_ ) )
                    {
                        static_person_found = true;
                        static_person_found_nmoves = nmoves;
                    }
                }
                if ( ! moving_person_found )
                {
                    if ( ( current_position_1.x_ == current_position_2.x_ ) && ( current_position_1.y_ == current_position_2.y_ ) )
                    {
                        moving_person_found = true;
                        moving_person_found_nmoves = nmoves;
                    }
                }
                ++nmoves;
                if ( nmoves > ( std::pow( size, 3 ) + 100000 ) )
                {
                    std::cout << "i = " << i << std::endl;
                    throw std::runtime_error( "Number of moves > " + int2string( std::pow( size, 3 ) ) );
                }
                current_position_1 = drunkards_walk_1.next_step();
                current_position_2 = drunkards_walk_2.next_step();
            }
            static_average_nmoves.add_value( static_cast<double>(static_person_found_nmoves) );
            moving_average_nmoves.add_value( static_cast<double>(moving_person_found_nmoves) );
        }
        std::cout << "Static average = " << static_average_nmoves.average() << ", moving average = " << moving_average_nmoves.average() << std::endl;
        std::cout << "Static ESD = " << static_average_nmoves.estimated_standard_deviation() << ", moving ESD = " << moving_average_nmoves.estimated_standard_deviation() << std::endl;
    MACRO_END_GAME

    try // CrystalStructure::convert_to_P1().
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        crystal_structure.convert_to_P1();
        crystal_structure.make_atom_labels_unique();
        crystal_structure.save_cif( append_to_file_name( input_file_name, "_P1" ) );
    MACRO_END_GAME

    try // Rescale unit-cell volume.
    {
        MACRO_ONE_FILELISTNAME_AS_ARGUMENT
        for ( size_t i( 0 ); i != file_list.size(); ++i )
        {
            CrystalStructure crystal_structure;
            std::cout << "Now reading cif... " + file_list.value(i).full_name() << std::endl;
            read_cif( file_list.value(i), crystal_structure );
            CrystalLattice crystal_lattice = crystal_structure.crystal_lattice();
            crystal_lattice.rescale_volume( 2586.833, 4 );
            crystal_structure.set_crystal_lattice( crystal_lattice );
            crystal_structure.save_cif( append_to_file_name( file_list.value(i), "_rv" ) );
        }
    MACRO_END_GAME

    try // Renumber cif files.
    {
        MACRO_ONE_FILELISTNAME_AS_ARGUMENT
        for ( int i( file_list.size() - 1 ); i > 0; --i )
        {
            if ( i < 8 )
                continue;
            size_t diff( 4 );
            if ( i < 23 )
                --diff;
            if ( i < 19 )
                diff -= 2;
            TextFileReader_2 input_file( file_list.value( i ) );
            if ( file_list.value( i ).extension() != "cif" )
                std::cout << "Warning: structure does not have .cif extension." << std::endl;
            TextFileWriter output_file( FileName( file_list.value( i ).directory(), "structure_" + size_t2string( i + diff + 1, 6, '0' ), file_list.value( i ).extension() ) );
            for ( size_t iLine( 0 ); iLine != input_file.size(); ++iLine )
            {
                if ( input_file.line( iLine ).substr( 0, 5 ) == "data_" )
                    output_file.write_line( "data_" + size_t2string( i + diff + 1, 6, '0' ) );
                else
                    output_file.write_line( input_file.line( iLine ) );
            }
        }
    MACRO_END_GAME

    try // Calculate dipole moment.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
//        SpaceGroup space_group;
//        space_group.add_inversion_at_origin();
//        crystal_structure.set_space_group( space_group );
        std::cout << crystal_structure.dipole_moment() << std::endl;
        // @@ The dipole moment of a periodic structure is not uniquely defined, the following generates nonsense.
        crystal_structure.apply_space_group_symmetry();
        std::cout << crystal_structure.dipole_moment() << std::endl;
    MACRO_END_GAME

    try // Calculate density for .cif file.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        crystal_structure.apply_space_group_symmetry();
        std::cout << crystal_structure.density() << std::endl;
    MACRO_END_GAME

    try // Calculate density for FileList.txt.
    {
        MACRO_ONE_FILELISTNAME_AS_ARGUMENT
        TextFileWriter text_file_writer( FileName( file_list_file_name.directory(), "densities", "txt" ) );
        for ( size_t i( 0 ); i != file_list.size(); ++i )
        {
            CrystalStructure crystal_structure;
            std::cout << "Now reading cif... " + file_list.value( i ).full_name() << std::endl;
            read_cif( file_list.value( i ), crystal_structure );
            crystal_structure.apply_space_group_symmetry();
            double density = crystal_structure.density();
            text_file_writer.write_line( FileName( "", file_list.value( i ).file_name(), file_list.value( i ).extension() ).full_name() + " " + double2string( density ) );
        }
    MACRO_END_GAME

    try // Extract SGE identifiers.
    {
        FileName input_file_name( "C:\\Data_Win\\SGE.txt" );
        TextFileReader text_file_reader( input_file_name );
        TextFileWriter output_file( append_to_file_name( input_file_name, "_ids" ) );
        std::vector< std::string > words;
        while ( text_file_reader.get_next_line( words ) )
        {
            if ( words.size() == 8 )
            {
                if ( ( words[ 3 ] == "research" ) && ( words[ 4 ] == "r" ) && ( words[ 7 ] == "1" ) )
                {
                    output_file.write( words[ 0 ] + " " );
                }
            }
        }
    MACRO_END_GAME

    try // Convert powder patterns in .txt format to .xye.
    {
        MACRO_ONE_FILELISTNAME_AS_ARGUMENT
        for ( size_t i( 0 ); i != file_list.size(); ++i )
        {
            FileName input_file_name = file_list.value(i);
            PowderPattern powder_pattern;
            powder_pattern.read_txt( input_file_name );
            powder_pattern.save_xye( replace_extension( input_file_name, "xye" ), true );
        }
    MACRO_END_GAME

    try // XIJKIK
    {
        FileName input_file_name( "C:\\GD\\XIJKIK.cif" );
        for ( size_t i( 1 ); i != 7; ++i )
        {
            for ( size_t j( 1 ); j != 7; ++j )
            {
                TextFileWriter text_file_writer( append_to_file_name( input_file_name, "_" + size_t2string(i) + "_" + size_t2string(j) ) );
                TextFileReader text_file_reader( input_file_name );
                std::string line;
                while ( text_file_reader.get_next_line( line ) )
                {
                    size_t iPos = line.find( "H_1" );
                    if ( iPos!= std::string::npos )
                    {
                        std::vector< std::string > words = split( line );
                        if ( ( words[0].substr( 3, 1 ) == size_t2string(i) )  )
                            text_file_writer.write_line( line );
                    }
                    else
                    {
                        iPos = line.find( "H_2" );
                        if ( iPos!= std::string::npos )
                        {
                            std::vector< std::string > words = split( line );
                            if ( ( words[0].substr( 3, 1 ) == size_t2string(j) )  )
                                text_file_writer.write_line( line );
                        }
                        else
                            text_file_writer.write_line( line );
                    }
                }
            }
        }
    MACRO_END_GAME

    try // Unit-cell transformation
    {
        CrystalLattice crystal_lattice(  9.912606,
                                        14.130999,
                                        18.496520,
                                        Angle::from_degrees( 90.064381 ),
                                        Angle::from_degrees( 72.695258 ),
                                        Angle::from_degrees( 71.861569 ) );
        // Make a copy
        CrystalLattice new_lattice( crystal_lattice );
        // Transform
        Matrix3D transformation_matrix( -1, 2, 0, -1, 0, 0, 0, -1, 1 );
        new_lattice.transform( transformation_matrix );
        new_lattice.print();
        std::cout << new_lattice.volume() << std::endl;
    MACRO_END_GAME

    try // Calculate unit-cell volume
    {
        CrystalLattice crystal_lattice( 14.872602,
                                        32.557974,
                                         5.085214,
                                        Angle::from_degrees( 90.0 ),
                                        Angle::from_degrees( 98.994 ),
                                        Angle::from_degrees( 90.0 ) );
        std::cout << crystal_lattice.volume() << std::endl;
    MACRO_END_GAME

    try // Analyse unit-cell parameter files from Materials Studio
    {
        // This should be a class and the class should also collect:
        // The temperatures
        // The experimental unit cell
        // The energy-minimised unit cell (= unit cell at T = 0 K)
        MACRO_ONE_FILELISTNAME_AS_ARGUMENT
        double start_time( 0.0 );
        TextFileWriter text_file_writer( FileName( file_list.base_directory(), "Analysis", "txt" ) );
        for ( size_t i( 0 ); i != file_list.size(); ++i )
        {
            RunningAverageAndESD<double> average_a;
            RunningAverageAndESD<double> average_b;
            RunningAverageAndESD<double> average_c;
            RunningAverageAndESD<Angle> average_alpha;
            RunningAverageAndESD<Angle> average_beta;
            RunningAverageAndESD<Angle> average_gamma;
            RunningAverageAndESD<double> average_volume;
            std::cout << "Now reading file... " + file_list.value( i ).full_name() << std::endl;
            TextFileReader text_file_reader( file_list.value( i ) );
            std::vector< std::string > words;
            size_t u = 6;
            size_t v = 2;
            size_t w = 12;
            text_file_reader.get_next_line( words );
            while ( text_file_reader.get_next_line( words ) )
            {
                if ( string2double( words[0] ) < start_time )
                    continue;
                // t, a, t, b, t, c, t, alpha, t, beta, t, gamma
                // 0  1  2  3  4  5  6    7    8    9   10   11
                //  "Length A" "Length B" "Length C" "alpha" "beta" "gamma"
                CrystalLattice crystal_lattice( string2double( words[ 0] ) / u,
                                                string2double( words[ 1] ) / v,
                                                string2double( words[ 2] ) / w,
                                                Angle::from_degrees( string2double( words[ 3] ) ),
                                                Angle::from_degrees( string2double( words[ 4] ) ),
                                                Angle::from_degrees( string2double( words[ 5] ) ) );
                // Transform the lattice
                crystal_lattice.transform( Matrix3D(  1.0,  0.0,  0.0,
                                                      0.0,  1.0,  0.0,
                                                      0.0,  0.0,  1.0 ) );
                average_a.add_value( crystal_lattice.a() );
                average_b.add_value( crystal_lattice.b() );
                average_c.add_value( crystal_lattice.c() );
                average_alpha.add_value( crystal_lattice.alpha() );
                average_beta.add_value( crystal_lattice.beta() );
                average_gamma.add_value( crystal_lattice.gamma() );
                average_volume.add_value( crystal_lattice.volume() );
            }
            text_file_writer.write_line( double2string( average_a.average() ) + " " +
                                         double2string( average_a.estimated_standard_deviation() ) + " " +
                                         double2string( average_b.average() ) + " " +
                                         double2string( average_b.estimated_standard_deviation() ) + " " +
                                         double2string( average_c.average() ) + " " +
                                         double2string( average_c.estimated_standard_deviation() ) + " " +
                                         double2string( average_alpha.average().value_in_degrees() ) + " " +
                                         double2string( average_alpha.estimated_standard_deviation().value_in_degrees() ) + " " +
                                         double2string( average_beta.average().value_in_degrees() ) + " " +
                                         double2string( average_beta.estimated_standard_deviation().value_in_degrees() ) + " " +
                                         double2string( average_gamma.average().value_in_degrees() ) + " " +
                                         double2string( average_gamma.estimated_standard_deviation().value_in_degrees() ) + " " +
                                         double2string( average_volume.average() ) + " " +
                                         double2string( average_volume.estimated_standard_deviation() )
                                       );
        }
    MACRO_END_GAME

    try // Materials Studio cell.xcd file
    {
        std::vector< FileName > file_names;
        file_names.push_back( FileName( "C:\\Data\\GD\\MD\\VI_6x2x12_50K_production Cell.xcd"  ) );
        file_names.push_back( FileName( "C:\\Data\\GD\\MD\\VI_6x2x12_100K_production Cell.xcd" ) );
        file_names.push_back( FileName( "C:\\Data\\GD\\MD\\VI_6x2x12_150K_production Cell.xcd" ) );
        file_names.push_back( FileName( "C:\\Data\\GD\\MD\\VI_6x2x12_200K_production Cell.xcd" ) );
        file_names.push_back( FileName( "C:\\Data\\GD\\MD\\VI_6x2x12_250K_production Cell.xcd" ) );
        file_names.push_back( FileName( "C:\\Data\\GD\\MD\\VI_6x2x12_300K_production Cell.xcd" ) );
        file_names.push_back( FileName( "C:\\Data\\GD\\MD\\VI_6x2x12_350K_production Cell.xcd" ) );
        file_names.push_back( FileName( "C:\\Data\\GD\\MD\\VI_6x2x12_400K_production Cell.xcd" ) );
        file_names.push_back( FileName( "C:\\Data\\GD\\MD\\VI_6x2x12_450K_production Cell.xcd" ) );
        file_names.push_back( FileName( "C:\\Data\\GD\\MD\\VI_6x2x12_500K_production Cell.xcd" ) );
        file_names.push_back( FileName( "C:\\Data\\GD\\MD\\VI_6x2x12_550K_production Cell.xcd" ) );
        file_names.push_back( FileName( "C:\\Data\\GD\\MD\\VI_6x2x12_600K_production Cell.xcd" ) );

        for ( size_t f( 0 ); f != file_names.size(); ++f )
        {
            TextFileReader text_file_reader( file_names[f] );
            TextFileWriter output_file( replace_extension( file_names[f], "txt" ) );
            std::vector< std::string > words;
            Splitter splitter_comma( "," );
            Splitter splitter_space( " " );
            std::vector< std::string > names;
            std::vector< std::vector< std::string > > values;
            size_t i( 0 );
            std::vector< std::string > value;
//<?xml version="1.0" encoding="latin1"?>
//<!DOCTYPE XCD []>
//<XCD Version="6.0" NumCharts="2">
//	<CHART_2D>
//		<DATA_2D NumSeries="3">
//			<SERIES_2D UniqueID="0" Name="Length A" NumPoints="863">
//				<POINT_2D XY="0.029,92.22462986173539"/>
//				<POINT_2D XY="0.058,92.2201391946457"/>
//				<POINT_2D XY="0.08699999999999999,92.2148028929542"/>
//				<POINT_2D XY="0.116,92.21203319627358"/>
//				<POINT_2D XY="0.145,92.20473408978023"/>
            bool first( true );
            while ( text_file_reader.get_next_line( words ) )
            {
                if ( words[0] == "<SERIES_2D" )
                {
                    std::string fragment = extract_delimited_text( text_file_reader.get_line(), "Name=\"", "\" NumPoints" );
                    names.push_back( fragment );
                    if ( ! first )
                    {
                        values.push_back( value );
                        value.clear();
                    }
                    first = false;
                }
                if ( words[0] == "<POINT_2D" )
                {
                    std::string fragment = extract_delimited_text( words[1], "\"", "\"" );
                    words = splitter_comma.split( fragment );
                    value.push_back( words[1] );
                }
            }
            values.push_back( value );
            output_file.write( "\"" + names[i] + "\"" );
            for ( size_t i( 1 ); i != names.size(); ++i )
            {
                output_file.write( std::string( " " ) + "\"" + names[i] + "\"" );
            }
            output_file.write_line();
            for ( size_t j( 0 ); j != values[i].size(); ++j )
            {
                for ( size_t i( 0 ); i != values.size(); ++i )
                {
                    output_file.write( " " + values[i][j] );
                }
                output_file.write_line();
            }

         }
    MACRO_END_GAME

    try // Write CASTEP input files for optimising H atoms (unit cell fixed, non-H fixed).
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        WriteCASTEPFile write_CASTEP_file( crystal_structure, input_file_name.directory(), input_file_name.file_name() );
        write_CASTEP_file.set_job_type( WriteCASTEPFile::H_ATOMS_ONLY );
        write_CASTEP_file.write();
    MACRO_END_GAME

    try // Run tests
    {
        run_tests();
    MACRO_END_GAME

    try // Write lean
    {
        MACRO_ONE_FILELISTNAME_AS_ARGUMENT
        file_list.set_prepend_file_name_with_basedirectory( true );
        for ( size_t i( 0 ); i != file_list.size(); ++i )
        {
            CrystalStructure crystal_structure;
            std::cout << "Now reading cif... " + file_list.value( i ).full_name() << std::endl;
            read_cif( file_list.value( i ), crystal_structure );
            crystal_structure.save_cif( append_to_file_name( file_list.value( i ), "_lean" ) );
        }
        TextFileWriter text_file_writer( file_list_file_name );
        for ( size_t i( 0 ); i != file_list.size(); ++i )
            text_file_writer.write_line( append_to_file_name( file_list.value( i ), "_lean" ).full_name() );
    MACRO_END_GAME

    try // Test DoubleWithESD
    {
        DoubleWithESD dwe_1( "10.81(8)" );
        double conversion_factor = 1.0/( 8.0 *CONSTANT_PI * CONSTANT_PI );
        DoubleWithESD dwe_2( conversion_factor * dwe_1.value(), conversion_factor * dwe_1.estimated_standard_deviation() );
        std::cout << dwe_2.crystallographic_style() << std::endl;
    MACRO_END_GAME

    try // Transform multiple ADPs
    {
        CrystalLattice crystal_lattice(  8.58,
                                         5.451,
                                        16.788,
                                        Angle::from_degrees( 90.0 ),
                                        Angle::from_degrees( 122.20 ),
                                        Angle::from_degrees( 90.0 ) );

        Matrix3D transformation_matrix(  1.0, 0.0, 1.0,
                                         0.0, 1.0, 0.0,
                                        -1.0, 0.0, 0.0  );

        std::vector< SymmetricMatrix3D > adps;
        adps.push_back( SymmetricMatrix3D(  0.14836,  0.05714,  0.08416, -0.00810,  0.08528, -0.01550 ) );
        adps.push_back( SymmetricMatrix3D(  0.14486,  0.08051,  0.12858,  0.02659,  0.11302,  0.01117 ) );
        adps.push_back( SymmetricMatrix3D(  0.09215,  0.07667,  0.12268,  0.03455,  0.07809,  0.03040 ) );
        adps.push_back( SymmetricMatrix3D(  0.05079,  0.05597,  0.06506,  0.00711,  0.03321,  0.01489 ) );
        adps.push_back( SymmetricMatrix3D(  0.05110,  0.03355,  0.04095, -0.00008,  0.02688,  0.00289 ) );
        adps.push_back( SymmetricMatrix3D(  0.07744,  0.04618,  0.04276, -0.01079,  0.03542, -0.00437 ) );
        adps.push_back( SymmetricMatrix3D(  0.07888,  0.08503,  0.04439, -0.03099,  0.02012, -0.00729 ) );
        adps.push_back( SymmetricMatrix3D(  0.05077,  0.09590,  0.06198, -0.00029,  0.01456,  0.02235 ) );
        adps.push_back( SymmetricMatrix3D(  0.20650,  0.10503,  0.08812, -0.02209,  0.10228, -0.03379 ) );
        adps.push_back( SymmetricMatrix3D(  0.19514,  0.05036,  0.12764, -0.01507,  0.11394, -0.02106 ) );
        adps.push_back( SymmetricMatrix3D(  0.22592,  0.12908,  0.19566,  0.05927,  0.17318,  0.00695 ) );
        adps.push_back( SymmetricMatrix3D(  0.14260,  0.11355,  0.13768,  0.01337,  0.11522,  0.02547 ) );
        adps.push_back( SymmetricMatrix3D(  0.13620,  0.07307,  0.14836,  0.04827,  0.09898,  0.04142 ) );
        adps.push_back( SymmetricMatrix3D(  0.09139,  0.13883,  0.17878,  0.04792,  0.09232,  0.04670 ) );
        adps.push_back( SymmetricMatrix3D(  0.06001,  0.06176,  0.07683, -0.00417,  0.04025,  0.01221 ) );
        adps.push_back( SymmetricMatrix3D(  0.07132,  0.03514,  0.05639, -0.00033,  0.03916,  0.00588 ) );
        adps.push_back( SymmetricMatrix3D(  0.08405,  0.05377,  0.04823, -0.00747,  0.03865,  0.00529 ) );
        adps.push_back( SymmetricMatrix3D(  0.11598,  0.15346,  0.04522, -0.05110,  0.01819, -0.01577 ) );
        adps.push_back( SymmetricMatrix3D(  0.11185,  0.08549,  0.08166, -0.05075,  0.04542, -0.01993 ) );
        adps.push_back( SymmetricMatrix3D(  0.08412,  0.10229,  0.07905,  0.01981,  0.02885,  0.04349 ) );
        adps.push_back( SymmetricMatrix3D(  0.05065,  0.17294,  0.09855, -0.00467,  0.01174,  0.03027 ) );
        TextFileWriter output_file( FileName( "C:\\Data\\GD\\adps.txt" ) );

        for ( size_t i( 0 ); i != adps.size(); ++i )
        {
            AnisotropicDisplacementParameters U_cif_new = transform_adps( adps[i], transformation_matrix, crystal_lattice );
            output_file.write_line( double2string( U_cif_new.value( 0, 0 ) ) + " " +
                                    double2string( U_cif_new.value( 1, 1 ) ) + " " +
                                    double2string( U_cif_new.value( 2, 2 ) ) + " " +
                                    double2string( U_cif_new.value( 0, 1 ) ) + " " +
                                    double2string( U_cif_new.value( 0, 2 ) ) + " " +
                                    double2string( U_cif_new.value( 1, 2 ) ) );
        }
    MACRO_END_GAME

    try // Test PowderPatternCalculator::calculate_equivalent_reflections()
    {
        FileName input_file_name( "C:\\Data\\for_testing\\P32_No_145.cif" );
        CrystalStructure crystal_structure;
        read_cif( input_file_name, crystal_structure );
        PowderPatternCalculator powder_pattern_calculator( crystal_structure );
        powder_pattern_calculator.calculate_reflection_list();
        ReflectionList reflection_list = powder_pattern_calculator.reflection_list();
        reflection_list.save( replace_extension( input_file_name, "txt" ) );
    MACRO_END_GAME

    try // Transform unit cell, NOT the atomic coordinates
    {
//        FileName input_file_name( "C:\\Data\\Tatiana\\Benzo\\261\\261_1.cif" );
//        CrystalStructure crystal_structure;
//        read_cif( input_file_name, crystal_structure );
//        CrystalLattice crystal_lattice = crystal_structure.crystal_lattice();

        CrystalLattice crystal_lattice( 8.5792,
                                        5.4508,
                                        16.4090,
                                        Angle::from_degrees( 90.0 ),
                                        Angle::from_degrees( 120.0311 ),
                                        Angle::from_degrees( 90.0 ) );

        Matrix3D transformation_matrix(  1.0,  0.0,  1.0,
                                         0.0,  1.0,  0.0,
                                        -1.0,  0.0,  0.0 );
        crystal_lattice.transform( transformation_matrix );
        transformation_matrix.invert();
        transformation_matrix.show();

        std::cout << crystal_lattice.a() << " " <<
                     crystal_lattice.b() << " " <<
                     crystal_lattice.c() << " " <<
                     crystal_lattice.alpha() << " " <<
                     crystal_lattice.beta() << " " <<
                     crystal_lattice.gamma() << std::endl;
    MACRO_END_GAME

    try // Unit cell matches
    {
        Angle two_theta_start( 5.0, Angle::DEGREES );
        Angle two_theta_end(  50.0, Angle::DEGREES );
        Angle two_theta_step( 0.01, Angle::DEGREES );
        double FWHM( 0.1 );
        // Create a dummy XRPD pattern
        CrystalLattice crystal_lattice_1( 14.07915,
                                          8.0043,
                                          11.820425,
                                          Angle::from_degrees( 90.0 ),
                                          Angle::from_degrees( 116.1443 ),
                                          Angle::from_degrees( 90.0 ) );
        std::vector< SymmetryOperator > symmetry_operators_1;
        symmetry_operators_1.push_back( SymmetryOperator( "x,y,z" ) );
        symmetry_operators_1.push_back( SymmetryOperator( "x,-y,1/2+z" ) );
        SpaceGroup space_group_1( symmetry_operators_1, "Pc" );
        CrystalStructure crystal_structure_1;
        crystal_structure_1.set_crystal_lattice( crystal_lattice_1 );
        crystal_structure_1.set_space_group( space_group_1 );
        PowderPatternCalculator powder_pattern_calculator_1( crystal_structure_1 );
        powder_pattern_calculator_1.set_wavelength( 1.54056 );
        powder_pattern_calculator_1.set_two_theta_start( two_theta_start );
        powder_pattern_calculator_1.set_two_theta_end( two_theta_end );
        powder_pattern_calculator_1.set_two_theta_step( two_theta_step );
        powder_pattern_calculator_1.set_FWHM( FWHM );
        powder_pattern_calculator_1.calculate_reflection_list();
        powder_pattern_calculator_1.set_structure_factors_to_1();
        PowderPattern powder_pattern_1;
        powder_pattern_calculator_1.calculate( powder_pattern_calculator_1.reflection_list(), powder_pattern_1 );
//        powder_pattern_calculator_1.reflection_list().save( FileName( "C:\\Data\\GF\\dummy_exp_pattern.hkl" ) );
//        powder_pattern_1.save( FileName( "C:\\Data\\Tatiana\\ADN\\dummy_exp_pattern.xye" ) );
        TextFileWriter text_file_writer( FileName( "C:\\Data\\GF\\matches.txt" ) );
        FileList file_list( FileName( "C:\\Data\\GF\\FileList.txt" ) );
        for ( size_t i( 0 ); i != file_list.size(); ++i )
        {
            CrystalStructure crystal_structure_2;
            std::cout << "Now reading cif... " + file_list.value( i ).full_name() << std::endl;
            read_cif( file_list.value( i ), crystal_structure_2 );
            PowderPatternCalculator powder_pattern_calculator_2( crystal_structure_2 );
            powder_pattern_calculator_2.set_wavelength( 1.54056 );
            powder_pattern_calculator_2.set_two_theta_start( two_theta_start );
            powder_pattern_calculator_2.set_two_theta_end( two_theta_end );
            powder_pattern_calculator_2.set_two_theta_step( two_theta_step );
            powder_pattern_calculator_2.set_FWHM( FWHM );
            powder_pattern_calculator_2.calculate_reflection_list();
            powder_pattern_calculator_2.set_structure_factors_to_1();
            PowderPattern powder_pattern_2;
            powder_pattern_calculator_2.calculate( powder_pattern_calculator_2.reflection_list(), powder_pattern_2 );
            double nwcc = normalised_weighted_cross_correlation( powder_pattern_1, powder_pattern_2, Angle( 3.0, Angle::DEGREES ) );
            text_file_writer.write_line( double2string( nwcc ) + " " + file_list.value( i ).full_name() );
        }
    MACRO_END_GAME

    try // Transform ADPs
    {
        SymmetricMatrix3D U_cif( 0.04682, 0.02399, 0.05938, 0.00258, -0.00735, -0.00165 );

        Matrix3D transformation_matrix(  1.0,  0.0,  0.0,
                                         0.0, -1.0,  0.0,
                                        -1.0,  0.0, -1.0 );

        CrystalLattice crystal_lattice( 13.8593,
                                        10.5242,
                                        14.8044,
                                        Angle::from_degrees( 90.0 ),
                                        Angle::from_degrees( 92.0014 ),
                                        Angle::from_degrees( 90.0 ) );

        AnisotropicDisplacementParameters U_cif_new = transform_adps( U_cif, transformation_matrix, crystal_lattice );
        U_cif_new.show();
    MACRO_END_GAME

    try // Analyse volumes in .cif file before and after minimisation
    {
        MACRO_ONE_FILELISTNAME_AS_ARGUMENT
        TextFileWriter text_file_writer( FileName( file_list_file_name.directory(), "volume_results", "txt") );
        for ( size_t i( 0 ); i != file_list.size(); ++i )
        {
            FileName exp_cif_file = file_list.value( i );
            FileName min_cif_file = append_to_file_name( exp_cif_file, "_mi_ucfr" );
            CrystalStructure exp_crystal_structure;
            read_cif( exp_cif_file, exp_crystal_structure );
            double exp_volume = exp_crystal_structure.crystal_lattice().volume();
            CrystalStructure min_crystal_structure;
            read_cif( min_cif_file, min_crystal_structure );
            double min_volume = min_crystal_structure.crystal_lattice().volume();
            text_file_writer.write_line( exp_cif_file.file_name() + " " + double2string( exp_volume ) + " " + double2string( min_volume ) );
        }
    MACRO_END_GAME

    try // Analyse similarities for Daniela
    {
        // Input: multiple lists of the same structures. They are compared and similarity matrices are constructed
        // All similarities are recorded.
        // Input: a FileList of FileList objects
        // The files must have a directory specified
        FileName file_list_of_file_lists_file_name;
        if ( (false) )
        {
            if ( argc != 2 )
                throw std::runtime_error( "Please give the name of a FileList of FileList files." );
            file_list_of_file_lists_file_name.set_full_name( argv[ 1 ] );
        }
        else
            file_list_of_file_lists_file_name.set_full_name( "C:\\Data\\SpaceGroupRevisions\\FileList.txt" );
        FileList file_list_of_file_lists( file_list_of_file_lists_file_name );
        if ( file_list_of_file_lists.empty() )
            throw std::runtime_error( std::string( "No files in file list " ) + file_list_of_file_lists_file_name.full_name() );
        TextFileWriter text_file_writer( FileName( file_list_of_file_lists.base_directory(), "Analysis_mi_vs_mi", "txt" ) );
        TextFileWriter text_file_writer_2( FileName( file_list_of_file_lists.base_directory(), "Analysis_mi_vs_exp", "txt" ) );
        std::vector< FileList > file_lists;
        file_lists.reserve( file_list_of_file_lists.size() );
        std::vector< std::string > directories;
        directories.reserve( file_list_of_file_lists.size() );
        for ( size_t i( 0 ); i != file_list_of_file_lists.size(); ++i )
        {
            file_lists.push_back( FileList( file_list_of_file_lists.value( i ) ) );
            file_lists[i].set_prepend_file_name_with_basedirectory( true );
            if ( file_lists[i].size() != file_lists[0].size() )
                throw std::runtime_error( "The FileLists have different sizes" );
            directories.push_back( file_list_of_file_lists.value( i ).directory() );
        }
        Splitter splitter( "_." );
        double overall_smallest_value( 1.0 );
        double overall_largest_value( 0.0 );
        std::string overall_smallest_identifier;
        std::string overall_largest_identifier;
        for ( size_t i( 0 ); i != file_lists[0].size(); ++i )
        {
            std::string identifier = file_lists[0].value( i ).file_name(); // This is now the identifier like "at2583_mi_ucfr"
            identifier = identifier.substr( 0, identifier.length() - std::string( "_mi_ucfr" ).length() );
            text_file_writer.write( identifier );
            text_file_writer_2.write( identifier );
            std::vector< CrystalStructure > crystal_structures;
            for ( size_t j( 0 ); j != file_list_of_file_lists.size(); ++j )
            {
                CrystalStructure crystal_structure;
                read_cif( file_lists[j].value( i ), crystal_structure );
                crystal_structures.push_back( crystal_structure );
            }
            // SimilarityMatrix expects numbers from 0 to 1, with 1 meaning identity.
            // Our numbers run from 0 to infinity, with 0 meaning identity...
            // The detour via a CorrelationMatrix currently does not serve any purpose any more
            CorrelationMatrix correlation_matrix( file_list_of_file_lists.size() );
            for ( size_t n( 0 ); n != file_list_of_file_lists.size(); ++n )
            {
                // Extract RMSCD vs experimental structure (bt2748_results.txt)
                double RMSCD_1;
                TextFileReader tfr( FileName( directories[n], identifier + "_results", "txt" ) );
                bool min2( false );
                std::string line;
                while ( tfr.get_next_line( line ) )
                {
                    if ( line == "Structural comparison exp-min2:" )
                        min2 = true;
                    if ( ! min2 )
                        continue;
                    if ( line.substr( 0, std::string( "RMS non-H Cartesian displacement [A]" ).length() ) == "RMS non-H Cartesian displacement [A]" )
                    {
                         std::vector< std::string > words = split( line );
                         RMSCD_1 = string2double( words[ 6 ] );
                    }
                }
                text_file_writer_2.write( " " + double2string( RMSCD_1 ) );
                for ( size_t m( n+1 ); m != file_list_of_file_lists.size(); ++m )
                {
                    double RMSCD = root_mean_square_Cartesian_displacement( crystal_structures[n], crystal_structures[m], false );
                    correlation_matrix.set_value( n, m, RMSCD );
                    text_file_writer.write( " " + double2string( RMSCD ) );
                }
            }
            text_file_writer.write_line();
            text_file_writer_2.write_line();
            double smallest_value = correlation_matrix.smallest_value();
            double largest_value = correlation_matrix.largest_value();
            if ( smallest_value < overall_smallest_value )
            {
                overall_smallest_value = smallest_value;
                overall_smallest_identifier = identifier;
            }
            if ( largest_value > overall_largest_value )
            {
                overall_largest_value = largest_value;
                overall_largest_identifier = identifier;
            }
        }
        text_file_writer.write_line();
        text_file_writer.write_line( overall_smallest_identifier + " " + double2string( overall_smallest_value ) );
        text_file_writer.write_line( overall_largest_identifier + " " + double2string( overall_largest_value ) );
    MACRO_END_GAME

    try // Analyse a whole T series at once: calculate average structures and ADPs from a set of frames (as cif files)
    {
        if ( argc != 2 )
            throw std::runtime_error( "Please give the name of a FileList of FileList files." );
        FileName file_list_of_file_lists_file_name( argv[ 1 ] );
        FileList file_list_of_file_lists( file_list_of_file_lists_file_name );
        if ( file_list_of_file_lists.empty() )
            throw std::runtime_error( std::string( "No files in file list " ) + file_list_of_file_lists_file_name.full_name() );
        for ( size_t i( 0 ); i != file_list_of_file_lists.size(); ++i )
        {
            FileList file_list( file_list_of_file_lists.value( i ) );
            // Guess temperature from file name, must find _nnnK_ or _nnnK.
            bool temperature_found( false );
            size_t temperature;
            Splitter splitter( "_." );
            std::vector< std::string > words = splitter.split( file_list_of_file_lists.value( i ).full_name() );
            for ( size_t j( 0 ) ; j != words.size(); ++j )
            {
                if ( words[j][ words[j].length()-1 ] == 'K' )
                {
                     try
                     {
                         temperature = string2integer( words[j].substr( 0, words[j].length()-1 ) );
                         temperature_found = true;
                     }
                     catch ( std::exception & e )
                     {
                     }
                }
            }
            file_list.set_prepend_file_name_with_basedirectory( true );
            SpaceGroup space_group; // = SpaceGroup::P21c();
//            space_group.add_inversion_at_origin(); // Space group is now P-1
            // Stores results as .cif
            AnalyseTrajectory analyse_trajectory( file_list, 6, 6, 8, space_group );
  //          analyse_trajectory.save_centres_of_mass();
        }
    MACRO_END_GAME

    try // ss-NMR
    {

// This is what we are looking for:
//                        <Atom3d ID="42" HasProperties="1" Mapping="104" Parent="2" Children="442,443,444,445" Name="Cl39" UserID="39" XYZ="-0.544406997224045,-0.119468460872243,-0.582538283348663" Connections="92" TemperatureType="Isotropic" AnisotropicTemperature="0,0,0,0,0,0,0,0,0" Components="Cl" Visible="0">
//                            <Properties NMRShielding="563.67" EFGQuadrupolarCoupling="-72.78" EFGAsymmetry="0.22" Force="-0.214441360892478,1.32464677237851,1.59841007334931"/>
//                        </Atom3d>
//
// If the user only assigned the NMR shieldings and not the EFG etc., the lines look like this:
//                        <Atom3d ID="4" HasProperties="1" Mapping="135" Parent="2" Children="137,393,394,395" Name="C1" UserID="1" XYZ="0.004362191714,0.464087446874054,0.655927752421914" Connections="69,70,100" TemperatureType="Isotropic" AnisotropicTemperature="0,0,0,0,0,0,0,0,0" Components="C">
//                            <Properties NMRShielding="8.93"/>
//                        </Atom3d>
//
// Beware of lines of this type:
//                        <Atom3d ID="566" Mapping="952" ImageOf="4"/>
//

    if ( argc != 2 )
       throw std::runtime_error( "Please give the name of a .xsd file that needs to be parsed as argument." );
    FileName input_file_name( argv[ 1 ] );
    TextFileReader text_file_reader( input_file_name );
    text_file_reader.set_skip_empty_lines( true );
    LabelsAndShieldings labels_and_shieldings;
    std::vector< std::string > words;
    while ( text_file_reader.get_next_line( words ) )
    {
        if ( words[0] == "<Atom3d" )
        {
            if ( ( words.size() == 4 ) && ( words[3].substr( words[3].length() - 2 ) == "/>" ) )
                continue;
            std::string total_line = text_file_reader.get_line();
            std::string line;
            size_t line_counter( 1 );
            while ( text_file_reader.get_next_line( line ) )
            {
                ++line_counter;
                line = strip( line );
                if ( line == "</Atom3d>" )
                    break;
                else
                    total_line += " " + line;
            }
            if ( line_counter != 3 )
                std::cout << "WARNING: number of lines for one atom not equal three." << std::endl;
            words = split( total_line );
            // Find "Name="
            bool name_found( false );
            std::string label;
            for ( size_t i( 0 ); i != words.size(); ++i )
            {
                if ( words[i].substr( 0, 5 ) == "Name=" )
                {
                    label = words[i].substr( 6, words[i].length() - 7 );
                    name_found = true;
                    break;
                }
            }
            if ( ! name_found )
                std::cout << "WARNING: name not found." << std::endl;
            // Find "NMRShielding="
            bool nmr_shielding_found( false );
            double shielding;
            for ( size_t i( 0 ); i != words.size(); ++i )
            {
                if ( words[i].substr( 0, 13 ) == "NMRShielding=" )
                {
                    if ( words[i].substr( words[i].length()-3, 3 ) == "\"/>" )
                        words[i] = words[i].substr( 0, words[i].length()-3 );
                    shielding = string2double( words[i].substr( 14, words[i].length() - 15 ) );
                    nmr_shielding_found = true;
                    break;
                }
            }
            if ( ! nmr_shielding_found )
                std::cout << "WARNING: NMR shielding not found." << std::endl;
            labels_and_shieldings.push_back( label, shielding );
        }
    }
    labels_and_shieldings.save( replace_extension( input_file_name, "out" ) );
    MACRO_END_GAME

    // This function only holds one frame at a time to allow in principle an infinite number of frames to be processed
    // without running out of memory.
    try // Average atoms of MD trajectory into single cif. Including ESDs.
    {
    if ( argc == 1 )
    {
        std::cout << "Usage:" << std::endl;
        std::cout << std::endl;
        std::cout << "AVGCIFS.exe <FileList.txt>" << std::endl;
        std::cout << std::endl;
        std::cout << "The argument <FileList.txt> is a file containing the names of the cif files." << std::endl;
        std::cout << "The file can be generated trivially by a command like:" << std::endl;
        std::cout << std::endl;
        std::cout << "dir /B *.cif > FileList.txt" << std::endl;
        std::cout << std::endl;
        std::cout << "Files with spaces in them must be enclosed in double quotes." << std::endl;
        std::cout << std::endl;
        std::cout << "If <FileList.txt> contains a path," << std::endl;
        std::cout << "then that path is taken to be the base directory for all cif files." << std::endl;
        std::cout << std::endl;
        std::cout << "Output: a file average.cif containing the average of each atom of all cif files." << std::endl;
        std::cout << "It is assumed that the number of atoms and the order of the atoms in all cif files is the same." << std::endl;
        std::cout << "The unit-cell parameters in average.cif are the averages of those in the cif files." << std::endl;
        char a;
        std::cin >> a;
        return 0;
    }
    FileName file_list_file_name( argv[ 1 ] );
    FileList file_list( file_list_file_name );
    if ( file_list.empty() )
        throw std::runtime_error( std::string( "No files in file list " ) + file_list_file_name.full_name() );
    file_list.set_prepend_file_name_with_basedirectory( true );
    std::vector< Element > elements;
    std::vector< RunningAverageAndESD< Vector3D > > positions;
    size_t natoms;
    RunningAverageAndESD< double > average_a;
    RunningAverageAndESD< double > average_b;
    RunningAverageAndESD< double > average_c;
    RunningAverageAndESD< Angle >  average_alpha;
    RunningAverageAndESD< Angle >  average_beta;
    RunningAverageAndESD< Angle >  average_gamma;
    RunningAverageAndESD< double > average_volume;
    {
    CrystalStructure crystal_structure;
    std::cout << "Now reading cif... " + file_list.value( 0 ).full_name() << std::endl;
    read_cif( file_list.value( 0 ), crystal_structure );
    crystal_structure.save_cif( FileName( file_list.value( 0 ).directory(), file_list.value( 0 ).file_name() + "_lean", "cif" ) );
    average_a.add_value( crystal_structure.crystal_lattice().a() );
    average_b.add_value( crystal_structure.crystal_lattice().b() );
    average_c.add_value( crystal_structure.crystal_lattice().c() );
    average_alpha.add_value( crystal_structure.crystal_lattice().alpha() );
    average_beta.add_value( crystal_structure.crystal_lattice().beta() );
    average_gamma.add_value( crystal_structure.crystal_lattice().gamma() );
    average_volume.add_value( crystal_structure.crystal_lattice().volume() );
    natoms = crystal_structure.natoms();
    elements.reserve( natoms );
    positions.reserve( natoms );
    for ( size_t i( 0 ); i != natoms; ++i )
    {
        elements.push_back( crystal_structure.atom( i ).element() );
        RunningAverageAndESD< Vector3D > position;
        position.add_value( crystal_structure.atom( i ).position() );
        positions.push_back( position );
    }
    }
    TextFileWriter unit_cell_esds_file( FileName( file_list.value( 0 ).directory(), "UnitCellESDs", "txt" ) );
    TextFileWriter atom_1_esds_file( FileName( file_list.value( 0 ).directory(), "Atom1ESDs", "txt" ) );
    for ( size_t i( 1 ); i != file_list.size(); ++i )
    {
        std::cout << "Now reading cif... " + file_list.value( i ).full_name() << std::endl;
        CrystalStructure crystal_structure;
        read_cif( file_list.value( i ), crystal_structure );
        if ( crystal_structure.natoms() != natoms )
            throw std::runtime_error( "The number of atoms in the cif files is not the same, the average cif could not be generated." );
        CrystalLattice crystal_lattice = crystal_structure.crystal_lattice();
        average_a.add_value( crystal_lattice.a() );
        average_b.add_value( crystal_lattice.b() );
        average_c.add_value( crystal_lattice.c() );
        average_alpha.add_value( crystal_lattice.alpha() );
        average_beta.add_value(  crystal_lattice.beta() );
        average_gamma.add_value( crystal_lattice.gamma() );
        average_volume.add_value( crystal_lattice.volume() );
        crystal_structure.save_cif( append_to_file_name( file_list.value( i ), "_lean" ) );
        unit_cell_esds_file.write_line( double2string( average_a.average() ) + " " +
                                        double2string( average_a.estimated_standard_deviation() ) + " " +
                                        double2string( average_b.average() ) + " " +
                                        double2string( average_b.estimated_standard_deviation() ) + " " +
                                        double2string( average_c.average() ) + " " +
                                        double2string( average_c.estimated_standard_deviation() ) + " " +
                                        double2string( average_alpha.average().value_in_degrees() ) + " " +
                                        double2string( average_alpha.estimated_standard_deviation().value_in_degrees() ) + " " +
                                        double2string( average_beta.average().value_in_degrees() ) + " " +
                                        double2string( average_beta.estimated_standard_deviation().value_in_degrees() ) + " " +
                                        double2string( average_gamma.average().value_in_degrees() ) + " " +
                                        double2string( average_gamma.estimated_standard_deviation().value_in_degrees() )
                                      );
        for ( size_t i( 0 ); i != natoms; ++i )
            positions[i].add_value( crystal_structure.atom( i ).position() );
        atom_1_esds_file.write_line( double2string( positions[0].average().x() ) + " " +
                                     double2string( positions[0].estimated_standard_deviation().x() ) + " " +
                                     double2string( positions[0].average().y() ) + " " +
                                     double2string( positions[0].estimated_standard_deviation().y() ) + " " +
                                     double2string( positions[0].average().z() ) + " " +
                                     double2string( positions[0].estimated_standard_deviation().z() )
                                   );
    }
    CrystalLattice crystal_lattice_average( average_a.average(),
                                            average_b.average(),
                                            average_c.average(),
                                            average_alpha.average(),
                                            average_beta.average(),
                                            average_gamma.average() );
        TextFileWriter text_file_writer( FileName( file_list.base_directory(), "average", "cif" ) );
        text_file_writer.write_line( "data_average" );
        text_file_writer.write_line( "_symmetry_space_group_name_H-M  'P 1'" );
        text_file_writer.write_line( "_symmetry_Int_Tables_number     1" );
        text_file_writer.write_line( "_symmetry_cell_setting          triclinic" );
        text_file_writer.write_line( "_cell_length_a    " + crystallographic_style( average_a.average(), average_a.estimated_standard_deviation() ) );
        text_file_writer.write_line( "_cell_length_b    " + crystallographic_style( average_b.average(), average_b.estimated_standard_deviation() ) );
        text_file_writer.write_line( "_cell_length_c    " + crystallographic_style( average_c.average(), average_c.estimated_standard_deviation() ) );
        text_file_writer.write_line( "_cell_angle_alpha " + crystallographic_style( average_alpha.average().value_in_degrees(), average_alpha.estimated_standard_deviation().value_in_degrees() ) );
        text_file_writer.write_line( "_cell_angle_beta  " + crystallographic_style( average_beta.average().value_in_degrees() , average_beta.estimated_standard_deviation().value_in_degrees()  ) );
        text_file_writer.write_line( "_cell_angle_gamma " + crystallographic_style( average_gamma.average().value_in_degrees(), average_gamma.estimated_standard_deviation().value_in_degrees() ) );
        text_file_writer.write_line( "_cell_volume      " + crystallographic_style( average_volume.average(), average_volume.estimated_standard_deviation() ) );
        text_file_writer.write_line( "loop_" );
        text_file_writer.write_line( "_symmetry_equiv_pos_site_id" );
        text_file_writer.write_line( "_symmetry_equiv_pos_as_xyz" );
        text_file_writer.write_line( "1 x,y,z"  );
        text_file_writer.write_line( "loop_" );
        text_file_writer.write_line( "_atom_site_label" ); // Needed for Materials Studio
        text_file_writer.write_line( "_atom_site_type_symbol" );
        text_file_writer.write_line( "_atom_site_fract_x" );
        text_file_writer.write_line( "_atom_site_fract_y" );
        text_file_writer.write_line( "_atom_site_fract_z" );
        size_t len( 2 );
        size_t current_size = 99;
        while ( natoms >= current_size )
        {
            ++len;
            current_size = 10 * current_size + 9;
        }
        for ( size_t i( 0 ); i != natoms; ++i )
        {
            // This is just too weird, need std::vector< DoubleWithESD > for this.
            text_file_writer.write_line( elements[i].symbol() + size_t2string( i + 1, len, '0' ) + " " +
                                         elements[i].symbol() + " " +
                                         crystallographic_style( positions[ i ].average().x(), positions[ i ].estimated_standard_deviation().x() ) + " " +
                                         crystallographic_style( positions[ i ].average().y(), positions[ i ].estimated_standard_deviation().y() ) + " " +
                                         crystallographic_style( positions[ i ].average().z(), positions[ i ].estimated_standard_deviation().z() ) );
        }
        text_file_writer.write_line();
        text_file_writer.write_line( "#END" );
    MACRO_END_GAME

    try // Collapse supercell
    {
        std::vector< SymmetryOperator > symmetry_operators;
        symmetry_operators.push_back( SymmetryOperator( Matrix3D(  1,  0,  0,
                                                                   0,  1,  0,
                                                                   0,  0,  1 ), Vector3D( 0.0, 0.0, 0.0 ) ) );
        symmetry_operators.push_back( SymmetryOperator( Matrix3D( -1,  0,  0,
                                                                   0,  1,  0,
                                                                   0,  0, -1 ), Vector3D( 0.5, 0.5, 0.5 ) ) );
        SpaceGroup space_group( symmetry_operators, "P21/n" );
        space_group.add_inversion_at_origin();
        CrystalLattice crystal_lattice( 7.0939, 9.2625, 11.657, Angle::angle_90_degrees(), Angle::from_degrees( 97.672 ), Angle::angle_90_degrees() );
        std::cout << "Step 1" << std::endl;
        CrystalStructure crystal_structure;
        read_cif( FileName( "C:\\Data\\ParacetamolMD\\Paracetamol_avg_uc_300K\\Paracetamol_fr0001.cif" ), crystal_structure );
        std::cout << "Step 2" << std::endl;
        crystal_structure.collapse_supercell( crystal_lattice, space_group );
        std::cout << "Step 3" << std::endl;
        crystal_structure.save_xyz( FileName( "C:\\Data\\collapsed.xyz" ) );
    MACRO_END_GAME

    try // Graeme's 50 ESI
    {
        FileName file_name_cif( "C:\\Data\\Graeme_50\\Graeme_50_ESI.txt" );
        TextFileReader text_file_reader( file_name_cif );
        TextFileWriter output_file( replace_extension( file_name_cif, "gcd" ) );
        std::vector< std::string > words;
        size_t iLine( 0 );
        while ( text_file_reader.get_next_line( words ) )
        {
            ++iLine;
            if ( is_odd( iLine ) )
            {
                output_file.write_line( words[0] );
            }
        }
    MACRO_END_GAME

    try // Test crystallographic_style()
    {
        {
            double value( 1.23456 );
            double esd(   0.099 );
            std::cout << crystallographic_style( value, esd ) << std::endl;
        }
        {
            double value( 1.23456 );
            double esd(   0.10 );
            std::cout << crystallographic_style( value, esd ) << std::endl;
        }
        {
            double value( 1.23456 );
            double esd(   0.011 );
            std::cout << crystallographic_style( value, esd ) << std::endl;
        }
        {
            double value( 123456 );
            double esd(        1.9 );
            std::cout << crystallographic_style( value, esd ) << std::endl;
        }
        {
            double value( 123456 );
            double esd(       50 );
            std::cout << crystallographic_style( value, esd ) << std::endl;
        }
        {
            double value( 123456 );
            double esd(       19 );
            std::cout << crystallographic_style( value, esd ) << std::endl;
        }
        // @@@ This one gives the wrong value.
        {
            double value( 123456 );
            double esd(      500 );
            std::cout << crystallographic_style( value, esd ) << std::endl;
        }
        {
            double value( 123456 );
            double esd(      190 );
            std::cout << crystallographic_style( value, esd ) << std::endl;
        }
        {
            double value( 1.23456 );
            double esd(   0.021 );
            std::cout << crystallographic_style( value, esd ) << std::endl;
        }
        {
            double value( 1.23456 );
            double esd(   0.006 );
            std::cout << crystallographic_style( value, esd ) << std::endl;
        }
    MACRO_END_GAME

    try // Test running average and ESD
    {
        RunningAverageAndESD<double> running_average_and_esd;
        std::vector< double > values;
        values.push_back( 10.1 );
        values.push_back( 10.2 );
        values.push_back(  9.9 );
        values.push_back( 10.0 );
        values.push_back(  9.8 );
        values.push_back( 10.1 );
        values.push_back(  9.9 );
        values.push_back( 10.0 );
        values.push_back( 10.0 );
        values.push_back( 10.1 );
//        values.push_back( 10.0 );
//        values.push_back( 10.0 );
//        values.push_back( 10.0 );
//        values.push_back( 10.0 );
//        values.push_back( 10.0 );
//        values.push_back( 10.0 );

        double average( 0.0 );
        for ( std::vector< double >::const_iterator it( values.begin() ); it != values.end(); ++it )
        {
            average += *it;
        }
        average /= values.size();
        double estimated_standard_deviation( 0.0 );
        for ( std::vector< double >::const_iterator it( values.begin() ); it != values.end(); ++it )
        {
            estimated_standard_deviation += square( *it - average );
        }
        estimated_standard_deviation /= values.size() - 1;
        estimated_standard_deviation = sqrt( estimated_standard_deviation );
        running_average_and_esd.add_values( values );
        RunningAverageAndESD<double> running_average_and_esd_2;
        for ( std::vector< double >::const_iterator it( values.begin() ); it != values.end(); ++it )
        {
            running_average_and_esd_2.add_value( *it );
        }
        std::cout << "average = " << running_average_and_esd.average() << std::endl;
        std::cout << "ESD = " << running_average_and_esd.estimated_standard_deviation() << std::endl;
        std::cout << "average = " << running_average_and_esd_2.average() << std::endl;
        std::cout << "ESD = " << running_average_and_esd_2.estimated_standard_deviation() << std::endl;
        std::cout << "average = " << average << std::endl;
        std::cout << "ESD = " << estimated_standard_deviation << std::endl;

        RunningAverageAndESD< Vector3D > running_average_and_esd_vector;
        running_average_and_esd_vector.add_value( Vector3D( 1.0, 1.0, 1.0 ) );

    MACRO_END_GAME

    try // Calculate normalised cross correlation for two cif files
    {
        Angle two_theta_start( 5.0, Angle::DEGREES );
        Angle two_theta_end(  50.0, Angle::DEGREES );
        Angle two_theta_step( 0.01, Angle::DEGREES );
        double FWHM( 0.1 );

        CrystalStructure crystal_structure_1;
        read_cif( FileName( "C:\\Data\\Refereeing\\GF.cif" ), crystal_structure_1 );
        PowderPatternCalculator powder_pattern_calculator_1( crystal_structure_1 );
        powder_pattern_calculator_1.set_wavelength( 1.54056 );
        powder_pattern_calculator_1.set_two_theta_start( two_theta_start );
        powder_pattern_calculator_1.set_two_theta_end( two_theta_end );
        powder_pattern_calculator_1.set_two_theta_step( two_theta_step );
        powder_pattern_calculator_1.set_FWHM( FWHM );
    //    powder_pattern_calculator_1.set_preferred_orientation( MillerIndices( 0, 0, 1 ), 0.6 );
        PowderPattern powder_pattern_1;
        powder_pattern_calculator_1.calculate( powder_pattern_1 );

        CrystalStructure crystal_structure_2;
        read_cif( FileName( "C:\\Data\\Refereeing\\GF.cif" ), crystal_structure_2 );
        PowderPatternCalculator powder_pattern_calculator_2( crystal_structure_2 );
        powder_pattern_calculator_2.set_wavelength( 1.54056 );
        powder_pattern_calculator_2.set_two_theta_start( two_theta_start );
        powder_pattern_calculator_2.set_two_theta_end( two_theta_end );
        powder_pattern_calculator_2.set_two_theta_step( two_theta_step );
        powder_pattern_calculator_2.set_FWHM( FWHM );
    //    powder_pattern_calculator_2.set_preferred_orientation( MillerIndices( 0, 0, 1 ), 0.6 );
        PowderPattern powder_pattern_2;
        powder_pattern_calculator_2.calculate( powder_pattern_2 );

        std::cout << "NWCC = " << normalised_weighted_cross_correlation( powder_pattern_1, powder_pattern_2, Angle( 3.0, Angle::DEGREES ) ) << std::endl;
        std::cout << "Rwp =  " << Rwp( powder_pattern_1, powder_pattern_2 ) << std::endl;
    MACRO_END_GAME

    try // Recalculate ESDs XRPD pattern
    {
        PowderPattern powder_pattern;
        powder_pattern.read_xye( FileName( "C:\\Data\\Files\\Beta\\Improved\\GF.xye" ) );
        powder_pattern.recalculate_estimated_standard_deviations();
        powder_pattern.save_xye( FileName( "C:\\Data\\Files\\Beta\\Improved\\GF_new_ESDs.xye" ), true );
    MACRO_END_GAME

    try // Calculate powder pattern from MD trajectory
    {
        if ( argc == 1 )
        {
            std::cout << "Usage:" << std::endl;
            std::cout << std::endl;
            std::cout << "MD2XRPD.exe <FileList.txt>" << std::endl;
            std::cout << std::endl;
            std::cout << "The argument <FileList.txt> is a file containing the names of the cif files." << std::endl;
            std::cout << "The file can be generated trivially by a command like:" << std::endl;
            std::cout << std::endl;
            std::cout << "dir /B *.cif > FileList.txt" << std::endl;
            std::cout << std::endl;
            std::cout << "Files with spaces in them must be enclosed in double quotes." << std::endl;
            std::cout << std::endl;
            std::cout << "If <FileList.txt> contains a path," << std::endl;
            std::cout << "then that path is taken to be the base directory for all cif files." << std::endl;
            char a;
            std::cin >> a;
            return 0;
        }
        FileName file_list_file_name( argv[ 1 ] );
        FileList file_list( file_list_file_name );
        if ( file_list.empty() )
            throw std::runtime_error( std::string( "No files in file list " ) + file_list_file_name.full_name() );
        Angle two_theta_start( 0.0, Angle::DEGREES );
        Angle two_theta_end(  60.0, Angle::DEGREES );
        Angle two_theta_step( 0.01, Angle::DEGREES );
        double FWHM( 0.1 );
        PowderPattern powder_pattern_sum( two_theta_start, two_theta_end, two_theta_step );
        for ( size_t i( 0 ); i != file_list.size(); ++i )
        {
            CrystalStructure crystal_structure;
            std::cout << "Now reading cif... " + file_list.value( i ).full_name() << std::endl;
            read_cif( file_list.value( i ), crystal_structure );
            std::cout << "Now calculating powder pattern... " + size_t2string( i, 4, '0' ) << std::endl;
            PowderPatternCalculator powder_pattern_calculator( crystal_structure );
            powder_pattern_calculator.set_wavelength( 1.54056 );
            powder_pattern_calculator.set_two_theta_start( two_theta_start );
            powder_pattern_calculator.set_two_theta_end( two_theta_end );
            powder_pattern_calculator.set_two_theta_step( two_theta_step );
            powder_pattern_calculator.set_FWHM( FWHM );
            PowderPattern powder_pattern;
            powder_pattern_calculator.calculate( powder_pattern );
            powder_pattern_sum += powder_pattern;
            powder_pattern.save_xye( FileName( file_list.base_directory(), "MD_fr" + size_t2string( i, 4, '0' ), "xye" ), true );
        }
        powder_pattern_sum.normalise_highest_peak();
        powder_pattern_sum.recalculate_estimated_standard_deviations();
        powder_pattern_sum.save_xye( FileName( file_list.base_directory(), "MD_sum_" + size_t2string( 0, 4, '0' )+"_"+size_t2string( file_list.size(), 4, '0' ), "xye" ), true );
    MACRO_END_GAME

    try // Calculate normalised cross correlation for two powder patterns
    {
        // GF_exp_no_background_5.0-39.98.xye
        // GF_Final_5.0-39.98.xye
        // GF_MD_5.0-39.98.xye
        // GF_mi_ucf_5.0-39.98.xye
        FileName file_name_1( "C:\\Data\\Tatiana\\GF.xye" );
        PowderPattern powder_pattern_1;
        PowderPattern powder_pattern_2;
        powder_pattern_1.read_xye( file_name_1 );
        powder_pattern_2.read_xye( FileName( "C:\\Data\\Tatiana\\GF_mi_ucf_5.0-39.98.xye" ) );
        std::cout << "NWCC = " << normalised_weighted_cross_correlation( powder_pattern_1, powder_pattern_2 ) << std::endl;
        std::cout << "Rwp =  " << Rwp( powder_pattern_1, powder_pattern_2 ) << std::endl;
        powder_pattern_2.read_xye( FileName( "C:\\Data\\Tatiana\\GF_MD_5.0-39.98.xye" ) );
        std::cout << "NWCC = " << normalised_weighted_cross_correlation( powder_pattern_1, powder_pattern_2 ) << std::endl;
        std::cout << "Rwp =  " << Rwp( powder_pattern_1, powder_pattern_2 ) << std::endl;
        powder_pattern_2.read_xye( FileName( "C:\\Data\\Tatiana\\GF_Final_5.0-39.98.xye" ) );
        std::cout << "NWCC = " << normalised_weighted_cross_correlation( powder_pattern_1, powder_pattern_2 ) << std::endl;
        std::cout << "Rwp =  " << Rwp( powder_pattern_1, powder_pattern_2 ) << std::endl;
    MACRO_END_GAME

    try // HMBENZ06
    {
        FileName file_name( "C:\\Data\\ActaCryst_powder\\RevisedSpaceGroups\\Old\\HMBENZ06_CSD.cif" );
        CrystalStructure crystal_structure;
        std::cout << "Now reading cif... " + file_name.full_name() << std::endl;
        read_cif( file_name, crystal_structure );
        crystal_structure.position_all_atoms_within_unit_cell();
        std::vector< SymmetryOperator > symmetry_operators;
        symmetry_operators.push_back( SymmetryOperator( Matrix3D( 1, 0, 0,
                                                                  0, 1, 0,
                                                                  0, 0, 1 ), Vector3D( 0, 0, 0 ) ) );
        symmetry_operators.push_back( SymmetryOperator( Matrix3D( 0, 0, 1,
                                                                  1, 0, 0,
                                                                  0, 1, 0 ), Vector3D( 0, 0, 0 ) ) );
        symmetry_operators.push_back( SymmetryOperator( Matrix3D( 0, 1, 0,
                                                                  0, 0, 1,
                                                                  1, 0, 0 ), Vector3D( 0, 0, 0 ) ) );
        std::cout << "Number of atoms in the asymmetric unit = " << crystal_structure.natoms() << std::endl;
        SpaceGroup space_group( symmetry_operators, "R-3" );
        space_group.add_inversion_at_origin();
        size_t old_natoms = crystal_structure.natoms();
        TextFileWriter text_file_writer( FileName( "C:\\Data\\ActaCryst_powder\\RevisedSpaceGroups\\HMBENZ06_CSD.out" ) );
        crystal_structure.set_space_group( space_group );
        for ( size_t i( 0 ); i != old_natoms; ++i )
        {
            Atom atom = crystal_structure.atom( i );
            for ( size_t j( 1 ); j != space_group.nsymmetry_operators(); ++j )
            {
                Atom new_atom( atom.element(), adjust_for_translations( space_group.symmetry_operator(j) * atom.position() ), atom.label() );
                crystal_structure.add_atom( new_atom );
            }
        }
        std::cout << "Number of atoms in the unit cell = " << crystal_structure.natoms() << std::endl;

        std::vector< bool > done( crystal_structure.natoms(), false );
        for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
        {
            if ( done[ i ] )
                continue;
            done[ i ] = true;
            std::vector< Atom > atoms;
            atoms.push_back( crystal_structure.atom( i ) );
            // Find all atoms that are close
            for ( size_t j( i+1 ); j != crystal_structure.natoms(); ++j )
            {
                if ( done[ j ] )
                    continue;
                if ( crystal_structure.crystal_lattice().shortest_distance( crystal_structure.atom( i ).position(), crystal_structure.atom( j ).position() ) < 0.5 )
                {
                    atoms.push_back( crystal_structure.atom( j ) );
                    done[ j ] = true;
                }
            }
                // Average the coordinates
        std::cout << atoms.size() << std::endl;
                Vector3D sum;
                for ( size_t k( 0 ); k != atoms.size(); ++k )
                {
                    sum += atoms[ k ].position();
                }
                sum = (1.0 / atoms.size()) * sum;
                text_file_writer.write_line( sum.to_string() );
        }
    MACRO_END_GAME

    try // Test Poisson generators.
    {
    TextFileWriter text_file_writer_exact( FileName( "C:\\Data\\for_testing\\Poisson_exact.txt" ) );
    TextFileWriter text_file_writer_Gaussian_0( FileName( "C:\\Data\\for_testing\\Poisson_Gaussian_0.txt" ) );
    TextFileWriter text_file_writer_Gaussian_1( FileName( "C:\\Data\\for_testing\\Poisson_Gaussian_1.txt" ) );
    TextFileWriter text_file_writer_Gaussian_2( FileName( "C:\\Data\\for_testing\\Poisson_Gaussian_2.txt" ) );

    double mean( 50.0 );
    size_t nsamples( 1000000 );

    std::cout << "Poisson_distribution_exact" << std::endl;
    for ( size_t i( 0 ); i != nsamples; ++i )
    {
        double value = Poisson_distribution_exact( mean );
        text_file_writer_exact.write_line( double2string( value ) );
    }
    std::cout << "Poisson_distribution_Gaussian_0" << std::endl;
    for ( size_t i( 0 ); i != nsamples; ++i )
    {
        double value = Poisson_distribution_Gaussian_0( mean );
        text_file_writer_Gaussian_0.write_line( double2string( value ) );
    }
    std::cout << "Poisson_distribution_Gaussian_1" << std::endl;
    for ( size_t i( 0 ); i != nsamples; ++i )
    {
        double value = Poisson_distribution_Gaussian_1( mean );
        if ( value > 10000 )
            std::cout << "1 value > 10000 = " << value << std::endl;
        text_file_writer_Gaussian_1.write_line( double2string( value ) );
    }
    std::cout << "Poisson_distribution_Gaussian_2" << std::endl;
    for ( size_t i( 0 ); i != nsamples; ++i )
    {
        double value = Poisson_distribution_Gaussian_2( mean );
        if ( value > 10000 )
            std::cout << "2 value > 10000 = " << value << std::endl;
        text_file_writer_Gaussian_2.write_line( double2string( value ) );
    }
    MACRO_END_GAME

    try // Add Poisson noise to powder pattern.
    {
        Angle two_theta_start( 0.0, Angle::DEGREES );
        Angle two_theta_end(  80.0, Angle::DEGREES );
        Angle two_theta_step( 0.01, Angle::DEGREES );
        double FWHM( 0.1 );
        CrystalStructure crystal_structure;
        read_cif( FileName( "C:\\Data\\for_testing\\HXACAN27_CSD_P1.cif" ), crystal_structure );
        PowderPatternCalculator powder_pattern_calculator( crystal_structure );
        powder_pattern_calculator.set_wavelength( 1.54056 );
        powder_pattern_calculator.set_two_theta_start( two_theta_start );
        powder_pattern_calculator.set_two_theta_end( two_theta_end );
        powder_pattern_calculator.set_two_theta_step( two_theta_step );
        powder_pattern_calculator.set_FWHM( FWHM );
        PowderPattern powder_pattern;
        powder_pattern_calculator.calculate( powder_pattern );
        powder_pattern.add_constant_background( 200.0 );
        powder_pattern.save_xye( FileName( "C:\\Data\\for_testing\\HXACAN27_CSD_P1.xye" ), true );
        powder_pattern.add_Poisson_noise();
        powder_pattern.save_xye( FileName( "C:\\Data\\for_testing\\HXACAN27_Poisson.xye" ), true );
    MACRO_END_GAME

    try // Test shortest distance in lattice with acute angle
    {
        CrystalLattice crystal_lattice( 1.0, 10.0, 1.0, Angle::angle_90_degrees(),
                                                        Angle::angle_90_degrees(),
                                                        Angle::from_degrees( 10.0 ) );
        for ( int i( -3 ); i != 4; ++i )
        {
            for ( int j( -3 ); j != 4; ++j )
            {
                Vector3D new_difference_vector( 0.9 + i, 0.9 + j, 0 );
                double distance = crystal_lattice.fractional_to_orthogonal( new_difference_vector ).length();
                std::cout << "distance = " << distance << ", i = " << i << ", j = " << j << std::endl;
            }
        }
        double shortest_distance = crystal_lattice.shortest_distance( Vector3D(), Vector3D( 0.9, 0.9, 0 ) );
        std::cout << "shortest_distance = " << shortest_distance << std::endl;
    MACRO_END_GAME

    try // apply_space_group_symmetry()
    {
        CrystalStructure crystal_structure;
        read_cif( FileName( "C:\\Data\\for_testing\\HXACAN27.cif" ), crystal_structure );
        crystal_structure.apply_space_group_symmetry();
        crystal_structure.save_cif( FileName( "C:\\Data\\for_testing\\HXACAN27_2.cif" ) );
    MACRO_END_GAME

    try // QAMXUY
    {
        TextFileReader text_file_reader_1( FileName( "C:\\Data\\ActaCryst_powder\\Errors\\QAMXUY\\QAMXUY_2.hkl" ) );
        ReflectionList reflection_list_1;
        std::vector< std::string > words;
        while ( text_file_reader_1.get_next_line( words ) )
        {
            if ( words.size() != 10 )
                throw std::runtime_error( "Wrong number of words" );
            // Miller indices, F^2, d-spacing, multiplicity
            reflection_list_1.push_back( MillerIndices( string2integer( words[0] ), string2integer( words[1] ), string2integer( words[2] ) ), string2double( words[5] ), string2double( words[8] ), 0 );
        }
        TextFileReader text_file_reader_2( FileName( "C:\\Data\\ActaCryst_powder\\Errors\\QAMXUY\\QAMXUY.hkl" ) );
        ReflectionList reflection_list_2;
        while ( text_file_reader_2.get_next_line( words ) )
        {
            if ( words.size() != 6 )
                throw std::runtime_error( "Wrong number of words" );
            // Miller indices, F^2, d-spacing, multiplicity
            reflection_list_2.push_back( MillerIndices( string2integer( words[0] ), string2integer( words[1] ), string2integer( words[2] ) ), string2double( words[4] ), string2double( words[3] ), string2integer( words[5] ) );
        }
        ReflectionList reflection_list_3;
        for ( size_t i( 0 ); i != reflection_list_1.size(); ++i )
        {
            MillerIndices miller_indices = reflection_list_1.miller_indices( i );
            size_t index = reflection_list_2.index( miller_indices );
            if ( index == reflection_list_2.size() )
                continue;
            // Miller indices, F^2, d-spacing, multiplicity
            reflection_list_3.push_back( miller_indices, reflection_list_1.F_squared( i ), reflection_list_1.d_spacing( i ), reflection_list_2.multiplicity( index ) );
        }
        Angle two_theta_start( 0.0, Angle::DEGREES );
        Angle two_theta_end(  80.0, Angle::DEGREES );
        Angle two_theta_step( 0.01, Angle::DEGREES );
        double FWHM( 0.1 );
        CrystalStructure crystal_structure;
        PowderPatternCalculator powder_pattern_calculator( crystal_structure );
        powder_pattern_calculator.set_wavelength( 1.54056 );
        powder_pattern_calculator.set_two_theta_start( two_theta_start );
        powder_pattern_calculator.set_two_theta_end( two_theta_end );
        powder_pattern_calculator.set_two_theta_step( two_theta_step );
        powder_pattern_calculator.set_FWHM( FWHM );
        PowderPattern powder_pattern;
        powder_pattern_calculator.calculate( reflection_list_3, powder_pattern );
        powder_pattern.save_xye( FileName( "C:\\Data\\ActaCryst_powder\\Errors\\QAMXUY\\QAMXUY.xye" ), true );
    MACRO_END_GAME

    try // Similarity transformation YIRVOL
    {

//        Matrix3D rotation( -1.0, -1.0, -1.0,
//                            2.0,  2.0,  1.0,
//                           -2.0, -1.0,  0.0 );
//        Vector3D r( 0.46306, 0.89770, 0.62438 );
//        r = rotation * r;
//        r = r + Vector3D( 0.25, 0.25, 0.75 );

        Matrix3D rotation(  1.0,  1.0, -1.0,
                            0.0,  0.0,  1.0,
                            0.0,  1.0,  0.0 );
        Vector3D r( 0.46306, 0.89770, 0.62438 );
        r = rotation * r;
        r = r + Vector3D( 0.25, 0.25, 0.75 );

        std::cout << r << std::endl;

        Vector3D translation( 0.0, 0.0, 0.0 );
        SymmetryOperator symmetry_operator( rotation, translation );
    MACRO_END_GAME

    try // Similarity transformation
    {
        // This is a bit silly because this is the unit-cell redefinition, not a rotation
        Matrix3D rotation(  1.0,  0.0,  0.0,
                            0.0,  1.0,  0.0,
                           -1.0, -1.0,  1.0 );
        Vector3D translation( 0.0, 0.0, 0.0 );
        SymmetryOperator symmetry_operator( rotation, translation );

        std::vector< SymmetryOperator > symmetry_operators;
        symmetry_operators.push_back( SymmetryOperator( "x,y,z" ) );
        symmetry_operators.push_back( SymmetryOperator( "y,x,-x-y-z" ) );
        symmetry_operators.push_back( SymmetryOperator( "-z+1/4,x+y+z+1/4,-x+1/4" ) );
        symmetry_operators.push_back( SymmetryOperator( "x+y+z+1/4,-z+1/4,-y+1/4" ) );
        SpaceGroup space_group( symmetry_operators, "Fdd2" );

        symmetry_operator.invert();
        space_group.apply_similarity_transformation( symmetry_operator );
        std::ofstream output_file;
        output_file.open( "C:\\Data\\symmetry_operators.txt" );
        if ( ! output_file )
           throw std::runtime_error( std::string( "Could not open file " ) );
        output_file << space_group;
    MACRO_END_GAME

    try // Coordinate rotation CASTEP mailing list
    {
        CrystalLattice crystal_lattice( 12.568, 12.568, 16.00, Angle::angle_90_degrees(), Angle::angle_90_degrees(), Angle::angle_120_degrees() );
        Matrix3D transformation_matrix = crystal_lattice.fractional_to_orthogonal_matrix();
        transformation_matrix.show();
        Matrix3D rotation_matrix = rotation_about_z( Angle::from_degrees( -60.0 ) );
        Matrix3D new_matrix = rotation_matrix * transformation_matrix;
        new_matrix.show();
        Vector3D p1( -1.570, 4.540, 5.050 );
        Vector3D p2(  1.570, 4.540, 5.050 );
        Vector3D p3(  0.000, 3.630, 0.950 );
        Matrix3D rotation_matrix_2 = rotation_about_z( Angle::from_degrees( 60.0 ) );
        Matrix3D matrix_2 = crystal_lattice.orthogonal_to_fractional_matrix() * rotation_matrix_2;
        p1 = matrix_2 * p1;
        p2 = matrix_2 * p2;
        p3 = matrix_2 * p3;
        std::cout << "p1 = " << p1 << std::endl;
        std::cout << "p2 = " << p2 << std::endl;
        std::cout << "p3 = " << p3 << std::endl;
        Vector3D r1( -0.3334794893476222, 0.0836386236569772, 0.3156250000000000 );
        Vector3D r2( -0.0836386236569798, 0.3334794893476196, 0.3156250000000000 );
        Vector3D r3( -0.1667555555555555, 0.1667555555555555, 0.0593750000000000 );
        r1 = new_matrix * r1;
        r2 = new_matrix * r2;
        r3 = new_matrix * r3;
        std::cout << "r1 = " << r1 << std::endl;
        std::cout << "r2 = " << r2 << std::endl;
        std::cout << "r3 = " << r3 << std::endl;
    MACRO_END_GAME

    try // WUBDOM
    {
    CrystalLattice crystal_lattice(  9.484,
                                    16.235,
                                     2.0*6.907,
                                    Angle::angle_90_degrees(),
                                    Angle::from_degrees( 118.25 ),
                                    Angle::angle_90_degrees() );
    std::vector< Atom > atoms;
    read_xyz( FileName( "C:\\Data\\ActaCryst_powder\\WUBDOM_2.xyz" ), atoms );
    TextFileWriter text_file_writer( FileName( "C:\\Data\\ActaCryst_powder\\WUBDOM_2_output.txt" ) );
    for ( size_t i( 0 ); i != atoms.size(); ++i )
        text_file_writer.write_line( atoms[ i ].element().symbol() + " " + crystal_lattice.orthogonal_to_fractional( atoms[ i ].position() ).to_string() );
    MACRO_END_GAME

    try // Convert MAQDAJ powder pattern to .xye
    {
        TextFileReader text_file_reader( FileName( "C:\\Data\\ActaCryst_powder\\Errors\\MAQDAJ.txt" ) );
        std::vector< std::string > words;
        PowderPattern powder_pattern;
        powder_pattern.set_wavelength( 1.15008 );
        size_t i( 0 );
        while ( text_file_reader.get_next_line( words ) )
        {
            powder_pattern.push_back( Angle::from_degrees( 8.0 + i*0.02 ), string2double( words[0] ) );
            ++i;
        }
        powder_pattern.save_xye( FileName( "C:\\Data\\ActaCryst_powder\\Errors\\MAQDAJ.xye" ), true );
    MACRO_END_GAME

    try // WIMWOE
    {
        CrystalLattice crystal_lattice( 6.62, 6.62, 9.31, Angle::angle_90_degrees(), Angle::angle_90_degrees(), Angle::angle_120_degrees() );
        std::vector< Vector3D > positions;
        positions.push_back( cylindrical2Cartesian( 1.24, Angle::from_degrees(    0.0 ), 0.00 ) );
        positions.push_back( cylindrical2Cartesian( 0.27, Angle::from_degrees(   12.0 ), 1.16 ) );
        positions.push_back( cylindrical2Cartesian( 1.04, Angle::from_degrees(  -75.3 ), 1.95 ) );
        positions.push_back( cylindrical2Cartesian( 1.14, Angle::from_degrees(  114.6 ), 1.34 ) );
        positions.push_back( cylindrical2Cartesian( 1.24, Angle::from_degrees( -120.0 ), 3.10 ) );
        positions.push_back( cylindrical2Cartesian( 2.64, Angle::from_degrees( -105.2 ), 3.40 ) );
        positions.push_back( cylindrical2Cartesian( 3.19, Angle::from_degrees(  -79.2 ), 2.86 ) );
        positions.push_back( cylindrical2Cartesian( 2.45, Angle::from_degrees(  -67.3 ), 1.66 ) );
        TextFileWriter text_file_writer( FileName( "C:\\Data\\ActaCryst_powder\\WIMWOE\\output.txt" ) );
        for ( size_t i( 0 ); i != positions.size(); ++i )
        {
            Vector3D result = crystal_lattice.orthogonal_to_fractional( positions[ i ] );
            std::cout << result << std::endl;
            text_file_writer.write_line( double2string( result.x() ) + " " + double2string( result.y() ) + " " + double2string( result.z() ) );
        }
    MACRO_END_GAME

    try // Skip Bo simulator
    {
        SkipBoGame skip_bo_game( 4, 20 );
        std::vector< SkipBoPlayer > players;
        players.push_back( SkipBoPlayer() );
        players.push_back( SkipBoPlayer() );
        players.push_back( SkipBoPlayer() );
        players.push_back( SkipBoPlayer() );
        CyclicInteger player( players.size() );
        --player;
        while ( ! players[player.next_value()].finished() )
        {
            ++player;
            players[player.next_value()].play( skip_bo_game );
        }
    MACRO_END_GAME

    try // Compare two C6 files. Only compare the elements they have in common.
    {
        TextFileReader text_file_reader_1( FileName( "C:\\Data\\C6ReferenceValues_old.txt" ) );
        text_file_reader_1.set_skip_empty_lines( true );
        text_file_reader_1.set_allow_single_quotes( true );
        TextFileReader text_file_reader_2( FileName( "C:\\Data\\C6ReferenceValues_new.txt" ) );
        text_file_reader_2.set_skip_empty_lines( true );
        text_file_reader_2.set_allow_single_quotes( true );
        std::vector< std::string > words;
        std::vector< C6Record > C6_records_old;
        std::vector< C6Record > C6_records_new;
        std::set< size_t > old_elements;
        std::set< size_t > new_elements;
        while ( text_file_reader_1.get_next_line( words ) )
        {
            C6Record C6_record;
            C6_record.element_1 = string2integer( words[0] );
            C6_record.element_2 = string2integer( words[1] );
            old_elements.insert( C6_record.element_1 );
            old_elements.insert( C6_record.element_2 );
            C6_record.CN_1 = string2double( words[2] );
            C6_record.CN_2 = string2double( words[3] );
            C6_record.C6 = string2double( words[4] );
            C6_records_old.push_back( C6_record );
        }
        while ( text_file_reader_2.get_next_line( words ) )
        {
            C6Record C6_record;
            C6_record.element_1 = string2integer( words[0] );
            C6_record.element_2 = string2integer( words[1] );
            new_elements.insert( C6_record.element_1 );
            new_elements.insert( C6_record.element_2 );
            C6_record.CN_1 = string2double( words[2] );
            C6_record.CN_2 = string2double( words[3] );
            C6_record.C6 = string2double( words[4] );
            C6_records_new.push_back( C6_record );
        }
        std::cout << old_elements.size() << std::endl;
        std::cout << new_elements.size() << std::endl;
    MACRO_END_GAME

    try // TRIZIN05, crystal structures of s-triazine
    {
        Angle phi( 15.2, Angle::DEGREES );
        Angle NCN( 125.2, Angle::DEGREES );
        Angle CNC = Angle::from_degrees( 240.0 ) - NCN;
        Angle hNCN = NCN / 2.0;
        Angle hCNC = CNC / 2.0;
        double c = 7.030;
        double CN( 1.338 );
        double CH( 0.96 );
        double y = 0.0;
        Vector3D O( 0.0, y, c/4.0 );
        Angle sixty( 60.0, Angle::DEGREES );
        double tan_60 = sixty.tangent();

        double C1_x = hCNC.sine() * CN * phi.cosine();
        double C1_y = ( hCNC.sine() / tan_60 ) * CN;
        double C1_z = phi.sine() * hCNC.sine() * CN;
        Vector3D C1( C1_x, C1_y, C1_z );
        double N1_x = hNCN.sine() * CN * phi.cosine();
        double N1_y = -( hNCN.sine() / tan_60 ) * CN;
        double N1_z = phi.sine() * hNCN.sine() * CN;
        Vector3D N1( N1_x, N1_y, N1_z );
        double C2_x = 0.0;
        double C2_y = -( hNCN.sine() / tan_60 ) * CN - hNCN.cosine() * CN;
        double C2_z = 0.0;
        Vector3D C2( C2_x, C2_y, C2_z );
        double N2_x = 0.0;
        double N2_y = ( hCNC.sine() / tan_60 ) * CN + hCNC.cosine() * CN;
        double N2_z = 0.0;
        Vector3D N2( N2_x, N2_y, N2_z );
        Vector3D H1 = C1 + C1 * ( CH / C1.length() );
        double H2_x = 0.0;
        double H2_y = C2_y - CH;
        double H2_z = 0.0;
        Vector3D H2( H2_x, H2_y, H2_z );
        C1 += O;
        N1 += O;
        C2 += O;
        N2 += O;
        H1 += O;
        H2 += O;

        CrystalLattice crystal_lattice( 6.719,
                                        9.528,
                                        c,
                                        Angle::from_degrees( 90.0 ),
                                        Angle::from_degrees( 125.28 ),
                                        Angle::from_degrees( 90.0 ) );
        Matrix3D orthogonal_to_fractional = crystal_lattice.for_CASTEP();
        orthogonal_to_fractional.transpose();
        orthogonal_to_fractional.invert();

        std::vector< SymmetryOperator > symmetry_operators;
        symmetry_operators.push_back( SymmetryOperator( "x,y,z" ) );
        symmetry_operators.push_back( SymmetryOperator( "1/2+x,1/2+y,z" ) );
        symmetry_operators.push_back( SymmetryOperator( "1/2-x,1/2+y,1/2-z" ) );
        symmetry_operators.push_back( SymmetryOperator( "-x,y,1/2-z" ) );
        symmetry_operators.push_back( SymmetryOperator( "-x,-y,-z" ) );
        symmetry_operators.push_back( SymmetryOperator( "1/2-x,1/2-y,-z" ) );
        symmetry_operators.push_back( SymmetryOperator( "1/2+x,1/2-y,1/2+z" ) );
        symmetry_operators.push_back( SymmetryOperator( "x,-y,1/2+z" ) );
        SpaceGroup space_group( symmetry_operators, "C2/c");
        CrystalStructure crystal_structure;
        crystal_structure.set_crystal_lattice( crystal_lattice );
        crystal_structure.set_space_group( space_group );
        crystal_structure.add_atom( Atom( Element( "C" ), orthogonal_to_fractional * C1, "C1" ) );
        crystal_structure.add_atom( Atom( Element( "N" ), orthogonal_to_fractional * N1, "N1" ) );
        crystal_structure.add_atom( Atom( Element( "C" ), orthogonal_to_fractional * C2, "C2" ) );
        crystal_structure.add_atom( Atom( Element( "N" ), orthogonal_to_fractional * N2, "N2" ) );
        crystal_structure.add_atom( Atom( Element( "H" ), orthogonal_to_fractional * H1, "H1" ) );
        crystal_structure.add_atom( Atom( Element( "H" ), orthogonal_to_fractional * H2, "H2" ) );
        crystal_structure.save_cif( FileName( "C:\\Data\\CCDC\\TRIZIN05.cif" ) );
    MACRO_END_GAME

    try // TRIZIN04, crystal structures of s-triazine
    {
        Angle phi( 11.1, Angle::DEGREES );
        Angle NCN( 125.2, Angle::DEGREES );
        Angle CNC = Angle::from_degrees( 240.0 ) - NCN;
        Angle hNCN = NCN / 2.0;
        Angle hCNC = CNC / 2.0;
        double c = 7.093;
        double CN( 1.338 );
        double CH( 0.96 );
        double y = -0.004;
        Vector3D O( 0.0, y, c/4.0 );
        Angle sixty( 60.0, Angle::DEGREES );
        double tan_60 = sixty.tangent();

        double C1_x = hCNC.sine() * CN * phi.cosine();
        double C1_y = ( hCNC.sine() / tan_60 ) * CN;
        double C1_z = phi.sine() * hCNC.sine() * CN;
        Vector3D C1( C1_x, C1_y, C1_z );
        double N1_x = hNCN.sine() * CN * phi.cosine();
        double N1_y = -( hNCN.sine() / tan_60 ) * CN;
        double N1_z = phi.sine() * hNCN.sine() * CN;
        Vector3D N1( N1_x, N1_y, N1_z );
        double C2_x = 0.0;
        double C2_y = -( hNCN.sine() / tan_60 ) * CN - hNCN.cosine() * CN;
        double C2_z = 0.0;
        Vector3D C2( C2_x, C2_y, C2_z );
        double N2_x = 0.0;
        double N2_y = ( hCNC.sine() / tan_60 ) * CN + hCNC.cosine() * CN;
        double N2_z = 0.0;
        Vector3D N2( N2_x, N2_y, N2_z );
        Vector3D H1 = C1 + C1 * ( CH / C1.length() );
        double H2_x = 0.0;
        double H2_y = C2_y - CH;
        double H2_z = 0.0;
        Vector3D H2( H2_x, H2_y, H2_z );
        C1 += O;
        N1 += O;
        C2 += O;
        N2 += O;
        H1 += O;
        H2 += O;

        CrystalLattice crystal_lattice( 6.884,
                                        9.569,
                                        c,
                                        Angle::from_degrees( 90.0 ),
                                        Angle::from_degrees( 126.61 ),
                                        Angle::from_degrees( 90.0 ) );
        Matrix3D orthogonal_to_fractional = crystal_lattice.for_CASTEP();
        orthogonal_to_fractional.transpose();
        orthogonal_to_fractional.invert();

        std::vector< SymmetryOperator > symmetry_operators;
        symmetry_operators.push_back( SymmetryOperator( "x,y,z" ) );
        symmetry_operators.push_back( SymmetryOperator( "1/2+x,1/2+y,z" ) );
        symmetry_operators.push_back( SymmetryOperator( "1/2-x,1/2+y,1/2-z" ) );
        symmetry_operators.push_back( SymmetryOperator( "-x,y,1/2-z" ) );
        symmetry_operators.push_back( SymmetryOperator( "-x,-y,-z" ) );
        symmetry_operators.push_back( SymmetryOperator( "1/2-x,1/2-y,-z" ) );
        symmetry_operators.push_back( SymmetryOperator( "1/2+x,1/2-y,1/2+z" ) );
        symmetry_operators.push_back( SymmetryOperator( "x,-y,1/2+z" ) );
        SpaceGroup space_group( symmetry_operators, "C2/c");
        CrystalStructure crystal_structure;
        crystal_structure.set_crystal_lattice( crystal_lattice );
        crystal_structure.set_space_group( space_group );
        crystal_structure.add_atom( Atom( Element( "C" ), orthogonal_to_fractional * C1, "C1" ) );
        crystal_structure.add_atom( Atom( Element( "N" ), orthogonal_to_fractional * N1, "N1" ) );
        crystal_structure.add_atom( Atom( Element( "C" ), orthogonal_to_fractional * C2, "C2" ) );
        crystal_structure.add_atom( Atom( Element( "N" ), orthogonal_to_fractional * N2, "N2" ) );
        crystal_structure.add_atom( Atom( Element( "H" ), orthogonal_to_fractional * H1, "H1" ) );
        crystal_structure.add_atom( Atom( Element( "H" ), orthogonal_to_fractional * H2, "H2" ) );
        crystal_structure.save_cif( FileName( "C:\\Data\\CCDC\\TRIZIN04.cif" ) );
    MACRO_END_GAME

    try // Disorder 02.
    {
        CrystalStructure crystal_structure;
        size_t n( 2 );
        crystal_structure.set_crystal_lattice( CrystalLattice( 10.0 * n,
                                                               10.0 * n,
                                                               10.0 * n,
                                                               Angle::angle_90_degrees(),
                                                               Angle::angle_90_degrees(),
                                                               Angle::angle_90_degrees() ) );
        bool ordered( true );
        RandomNumberGenerator_double rng;
        std::vector< std::vector< bool > > initial_configurations; // Stores the configuration at the "start" of the column, which is 50:50
        for ( size_t iX( 0 ); iX != n; ++iX )
        {
            std::vector< bool > temp_initial_configurations;
            for ( size_t iY( 0 ); iY != n; ++iY )
            {
                if ( ordered )
                    temp_initial_configurations.push_back( true );
                else
                    temp_initial_configurations.push_back( rng.next_number() < 0.5 );
            }
            initial_configurations.push_back( temp_initial_configurations );
        }
        // Cartesian coordinates in cell 0,0,0
        Vector3D r1_1( 5.0, 2.5, 2.5 );
        Vector3D r2_1( 5.0, 2.5, 7.5 );
        Vector3D r3_1( 5.0, 7.5, 2.5 );
        Vector3D r4_1( 5.0, 7.5, 7.5 );

        Vector3D r1_2( 2.5, 5.0, 2.5 );
        Vector3D r2_2( 2.5, 5.0, 7.5 );
        Vector3D r3_2( 7.5, 5.0, 2.5 );
        Vector3D r4_2( 7.5, 5.0, 7.5 );
        for ( size_t iX( 0 ); iX != n; ++iX )
        {
            for ( size_t iY( 0 ); iY != n; ++iY )
            {
                for ( size_t iZ( 0 ); iZ != n; ++iZ )
                {
                    bool current_configuration = initial_configurations[ iX ][ iY ];
                    if ( is_odd( iZ ) )
                        current_configuration = ( ! current_configuration );
                    if ( current_configuration )
                    {
                        {
                        Atom atom( Element( "C" ), crystal_structure.crystal_lattice().orthogonal_to_fractional( r1_1 + Vector3D( iX*10.0, iY*10.0, iZ*10.0 ) ), "C1" );
                        crystal_structure.add_atom( atom );
                        }
                        {
                        Atom atom( Element( "C" ), crystal_structure.crystal_lattice().orthogonal_to_fractional( r2_1 + Vector3D( iX*10.0, iY*10.0, iZ*10.0 ) ), "C1" );
                        crystal_structure.add_atom( atom );
                        }
                        {
                        Atom atom( Element( "C" ), crystal_structure.crystal_lattice().orthogonal_to_fractional( r3_1 + Vector3D( iX*10.0, iY*10.0, iZ*10.0 ) ), "C1" );
                        crystal_structure.add_atom( atom );
                        }
                        {
                        Atom atom( Element( "C" ), crystal_structure.crystal_lattice().orthogonal_to_fractional( r4_1 + Vector3D( iX*10.0, iY*10.0, iZ*10.0 ) ), "C1" );
                        crystal_structure.add_atom( atom );
                        }
                    }
                    else
                    {
                        {
                        Atom atom( Element( "C" ), crystal_structure.crystal_lattice().orthogonal_to_fractional( r1_2 + Vector3D( iX*10.0, iY*10.0, iZ*10.0 ) ), "C1" );
                        crystal_structure.add_atom( atom );
                        }
                        {
                        Atom atom( Element( "C" ), crystal_structure.crystal_lattice().orthogonal_to_fractional( r2_2 + Vector3D( iX*10.0, iY*10.0, iZ*10.0 ) ), "C1" );
                        crystal_structure.add_atom( atom );
                        }
                        {
                        Atom atom( Element( "C" ), crystal_structure.crystal_lattice().orthogonal_to_fractional( r3_2 + Vector3D( iX*10.0, iY*10.0, iZ*10.0 ) ), "C1" );
                        crystal_structure.add_atom( atom );
                        }
                        {
                        Atom atom( Element( "C" ), crystal_structure.crystal_lattice().orthogonal_to_fractional( r4_2 + Vector3D( iX*10.0, iY*10.0, iZ*10.0 ) ), "C1" );
                        crystal_structure.add_atom( atom );
                        }
                    }
                }
            }
        }
        crystal_structure.save_cif( FileName( "C:\\Users\\jacco\\Documents\\disorder_04.cif" ) );

        if ( false )
        {
            crystal_structure.apply_space_group_symmetry();
            std::cout << "Now calculating powder pattern... " << std::endl;
            PowderPatternCalculator powder_pattern_calculator( crystal_structure );
            Angle two_theta_start( 1.0, Angle::DEGREES );
            Angle two_theta_end(  40.0, Angle::DEGREES );
            Angle two_theta_step( 0.015, Angle::DEGREES );
            double wavelength( 1.54056 );
            powder_pattern_calculator.set_wavelength( wavelength );
            powder_pattern_calculator.set_two_theta_start( two_theta_start );
            powder_pattern_calculator.set_two_theta_end( two_theta_end );
            powder_pattern_calculator.set_two_theta_step( two_theta_step );
            powder_pattern_calculator.set_FWHM( 0.20 );
            PowderPattern powder_pattern;
            powder_pattern_calculator.calculate( powder_pattern );
            powder_pattern.save_xye( FileName( "C:\\Users\\jacco\\Documents\\disorder_02.xye" ), true );
        }
    MACRO_END_GAME

    try // Disorder 01.
    {
        CrystalStructure crystal_structure;
        crystal_structure.set_crystal_lattice( CrystalLattice( 10.0,
                                                               10.0,
                                                               10.0,
                                                               Angle::angle_90_degrees(),
                                                               Angle::angle_90_degrees(),
                                                               Angle::angle_90_degrees() ) );
        {
        Atom atom( Element( "C" ), Vector3D( 0.5, 0.25, 0.25 ), "C1" );
        atom.set_occupancy( 0.5 );
        crystal_structure.add_atom( atom );
        }
        {
        Atom atom( Element( "C" ), Vector3D( 0.5, 0.25, 0.75 ), "C2" );
        atom.set_occupancy( 0.5 );
        crystal_structure.add_atom( atom );
        }
        {
        Atom atom( Element( "C" ), Vector3D( 0.5, 0.75, 0.25 ), "C3" );
        atom.set_occupancy( 0.5 );
        crystal_structure.add_atom( atom );
        }
        {
        Atom atom( Element( "C" ), Vector3D( 0.5, 0.75, 0.75 ), "C4" );
        atom.set_occupancy( 0.5 );
        crystal_structure.add_atom( atom );
        }

        {
        Atom atom( Element( "C" ), Vector3D( 0.25, 0.5, 0.25 ), "C5" );
        atom.set_occupancy( 0.5 );
        crystal_structure.add_atom( atom );
        }
        {
        Atom atom( Element( "C" ), Vector3D( 0.25, 0.5, 0.75 ), "C6" );
        atom.set_occupancy( 0.5 );
        crystal_structure.add_atom( atom );
        }
        {
        Atom atom( Element( "C" ), Vector3D( 0.75, 0.5, 0.25 ), "C7" );
        atom.set_occupancy( 0.5 );
        crystal_structure.add_atom( atom );
        }
        {
        Atom atom( Element( "C" ), Vector3D( 0.75, 0.5, 0.75 ), "C8" );
        atom.set_occupancy( 0.5 );
        crystal_structure.add_atom( atom );
        }
        crystal_structure.save_cif( FileName( "C:\\Users\\jacco\\Documents\\disorder_01.cif" ) );
        crystal_structure.apply_space_group_symmetry();

        std::cout << "Now calculating powder pattern... " << std::endl;
        PowderPatternCalculator powder_pattern_calculator( crystal_structure );
        Angle two_theta_start( 1.0, Angle::DEGREES );
        Angle two_theta_end(  40.0, Angle::DEGREES );
        Angle two_theta_step( 0.015, Angle::DEGREES );
        double wavelength( 1.54056 );
        powder_pattern_calculator.set_wavelength( wavelength );
        powder_pattern_calculator.set_two_theta_start( two_theta_start );
        powder_pattern_calculator.set_two_theta_end( two_theta_end );
        powder_pattern_calculator.set_two_theta_step( two_theta_step );
        powder_pattern_calculator.set_FWHM( 0.20 );
        PowderPattern powder_pattern;
        powder_pattern_calculator.calculate( powder_pattern );
        powder_pattern.save_xye( FileName( "C:\\Users\\jacco\\Documents\\disorder_01.xye" ), true );
    MACRO_END_GAME

    try // Pythagoras pie problem.
    {
        double total_left( 1.0 );
        double largest_share( 0.0 );
        size_t guest_with_largest_share( 0 );
        for ( size_t i( 1 ); i != 101; ++i )
        {
            double this_percentage( i / 100.0 );
            std::cout << "this_percentage = " << this_percentage << std::endl;
            double this_share = this_percentage * total_left;
            std::cout << "this_share = " << this_share << std::endl;
            total_left -= this_share;
            std::cout << "total_left = " << total_left << std::endl;
            if ( this_share > largest_share )
            {
                largest_share = this_share;
                guest_with_largest_share = i;
            }
        }
        std::cout << "largest_share = " << largest_share << std::endl;
        std::cout << "guest_with_largest_share = " << guest_with_largest_share << std::endl;
        std::cout << "total_left = " << total_left << std::endl;
    MACRO_END_GAME

    try // Compile all headers stand-alone
    {
        FileName file_list_file_name( "C:\\Cpp\\FileList.txt" );
        FileList file_list( file_list_file_name );
        file_list.set_prepend_file_name_with_basedirectory( false );
        if ( file_list.empty() )
            throw std::runtime_error( std::string( "No files in file list " ) + file_list_file_name.full_name() );
        for ( size_t i( 0 ); i != file_list.size(); ++i )
        {
            TextFileWriter text_file_writer( FileName( "C:\\Compile_Headers", file_list.value( i ).file_name() + "_header", "cpp" ) );
            text_file_writer.write_line( "#include \"" + file_list.value( i ).file_name() + ".h\"" );
        }
    MACRO_END_GAME

}
