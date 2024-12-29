/* *********************************************
Copyright (c) 2013-2025, Cornelis Jan (Jacco) van de Streek
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
#include "Centring.h"
#include "ChebyshevBackground.h"
#include "CheckFoundItem.h"
#include "ChemicalFormula.h"
#include "CollectionOfPoints.h"
#include "Complex.h"
#include "Constraints.h"
#include "CopyTextFile.h"
#include "CorrelationMatrix.h"
#include "CrystallographicCalculations.h"
#include "CrystalStructure.h"
#include "CrystalStructuresDatabase.h"
#include "CyclicInteger.h"
#include "DoublesAsTable.h"
#include "DoubleWithESD.h"
#include "DrunkardsWalk.h"
#include "Eigenvalue.h"
#include "EndGame.h"
#include "FileList.h"
#include "FileName.h"
#include "FingerCoxJephcoat.h"
#include "FingerCoxJephcoat_functions.h"
#include "Finish_inp.h"
#include "Fraction.h"
#include "GCDFile.h"
#include "GeneratePowderCIF.h"
#include "Histogram.h"
#include "InpWriter.h"
#include "LabelsAndShieldings.h"
#include "LinearRegression.h"
#include "Mapping.h"
#include "MathsFunctions.h"
#include "MatrixFraction3D.h"
#include "MC_alkanes.h"
#include "ModelBuilding.h"
#include "OrientationalOrderParameters.h"
#include "Plane.h"
#include "PowderPattern.h"
#include "PowderPatternCalculator.h"
#include "RandomNumberGenerator.h"
#include "ReadCell.h"
#include "ReadCif.h"
#include "ReadCifOrCell.h"
#include "ReadXSD.h"
#include "ReadXYZ.h"
#include "RealisticXRPDSimulator.h"
#include "RealisticXRPDSimulatorSettings.h"
#include "Refcode.h"
#include "RefcodeList.h"
#include "ReflectionList.h"
#include "RunningAverageAndESD.h"
#include "RunTests.h"
#include "SimilarityAnalysis.h"
#include "SkipBo.h"
#include "Sort.h"
#include "SphericalHarmonics.h"
#include "StringFunctions.h"
#include "Sudoku.h"
#include "SudokuSolver.h"
#include "SymmetryOperator.h"
#include "Tally.h"
#include "TextFile.h"
#include "TextFileReader.h"
#include "TextFileReader_2.h"
#include "TextFileWriter.h"
#include "TLS.h"
#include "TLSWriter.h"
#include "TOPAS.h"
#include "Utilities.h"
#include "Vector3D.h"
#include "Vector3DCalculations.h"
#include "VoidsFinder.h"
#include "WriteCASTEPFile.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#define MACRO_ONE_FILELISTNAME_AS_ARGUMENT \
        if ( argc != 2 ) \
            throw std::runtime_error( "Please give the name of a FileList.txt file." ); \
        FileName file_list_file_name( argv[ 1 ] ); \
        FileList file_list( file_list_file_name ); \
        if ( file_list.empty() ) \
            throw std::runtime_error( std::string( "No files in file list " ) + file_list_file_name.full_name() );

#define MACRO_ONE_FILELISTNAME_OR_LIST_OF_FILES_AS_ARGUMENT \
        if ( argc == 1 ) \
            throw std::runtime_error( "Please give the name of a FileList.txt file." ); \
        FileList file_list; \
        if ( argc == 2 ) \
        { \
            FileName file_list_file_name( argv[ 1 ] ); \
            file_list = FileList( file_list_file_name ); \
            if ( file_list.empty() ) \
                throw std::runtime_error( std::string( "No files in file list " ) + file_list_file_name.full_name() ); \
        } \
        else \
        { \
            std::vector< FileName > files; \
            for ( int i( 1 ); i != argc; ++i ) \
                files.push_back( FileName( argv[ i ] ) ); \
            file_list = FileList( files ); \
        }

#define MACRO_LIST_OF_FILES_AS_ARGUMENT \
        if ( argc == 1 ) \
            throw std::runtime_error( "Please give the names of one or more files as argument." ); \
        FileList file_list; \
        std::vector< FileName > files; \
        for ( int i( 1 ); i != argc; ++i ) \
            files.push_back( FileName( argv[ i ] ) ); \
        file_list = FileList( files );

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

// Test if .h files compile stand alone.
// FractionalCoordinate / OrthonormalCoordinate.
// Smart pointers. Because I currently have no smart pointers, it is almost impossible to hold objects by pointer so everything is always copied.
// Add class "CheckedItemReadFromFile" or "CheckedItemAssignedValue" Store reference (&) to an existing variable
// + bool, upon destruction writes error message if variable has not been assigned a value.
// Maybe a "format checker and beautifier", e.g. an atom label must be "Aa11", a refcode must be "AAAAAA11".
// A "Distance" and "Length" class, which for < and nearly_equal() uses the norm2 and calculates the norm on demand (so cached).
// void check_if_quotes_correct( const std::string & input );
// Need a "add_two_strings_with_quotes()".
// Alternatively, a QuotedString class.
// Add Marcus' cell deformation measure.
// Number of molecules for hemihydrate must be 1 API plus 1/2 H2O.
// Could now add CrystalStructure::CrystalStructure( const FileName & ).
// FileList::save().
// A "DoubleChecked" class: check that adding two numbers changes the result (i.e. check if one is too small to be added to the other),
//    a bool to keep track whether the value has been initialised. Function to convert everything smaller than e.g. 1E-6 to 0.0.
//    add specialised function for adding a std::vector of numbers. Is there anything special that can be done by alternating positive and negative numbers?
//    Dimension (as in the unit, kg, meter, second) as std::string.
//    Internal bool for "very big" (infinity), although "not initialised" may fulfil that purpose.
// Topological equivalence of torsion angles.
// Two kinds of crystal structure: asymmetric unit only and all molecules in one unit cell.
// Add appending + checks for existence to TextFileWriter.
// Remove trailing whitespace from text file.
// The "finish .inp after DASH 6" could reshuffle some keywords (move 2theta step to just after 2theta begin and 2theta end, for example) and beautify some others.
// Replace Rietveld by Loopstra-Rietveld.
// Maybe a file, or class, CrystallographicCalculator.
//    The space group and the lattice are always needed when e.g. calculating if two reflections are equivalent
//    (currently member functions of PowderPatternCalculator).
// Could implement a "VeryBigInteger" class by storing a std::vector< size_t > which stores e.g. the number 8242 as [ 8, 2, 4, 2 ]. MAXINT = 2,147,483,647
// After transformation, check that the structure has not changed by calculating the XRPD patterns.
// Simulate a single crystal diffraction pattern.
// Find all O-H and N-H, remove all those H, find all hydrogen bonds between all O, N and Cl-.
// Equality operator for chemical formulae that can ignore order and D vs H.
// std::vector< double > occurs so often that it would make sense if it had its own class.
//    Then also a classs that holds two such std::vector< double > objects, of the same length.
//    Then also a classs that holds three such std::vector< double > objects, of the same length.
// Ideally, we should change show() to to_string() eveywhere. This has three advantages. 1. We do not have to include the expensive
//    iostream header everywhere. 2. We cannot write test suite code for show() but we can test to_string(). 3. It is more flexible. E.g. show()
//    forces a std::endl upon the user. There is only one disadvantage, show() can distribute an object over multiple lines, this would be
//    unexpected for to_string().
//    Could create a show( class ) for every class in a separate file to avoid having to include <iostream> or <string> in classes like Matrix3D and Angle.
// nearly_zero( Complex ) and nearly_zero( Vector3D ) are *NOT* implemented as nearly_zero( norm2() ) but as
//    nearly_zero( real() ) && nearly_zero( imaginary() ) and nearly_zero( x() ) && nearly_zero( y() ) && nearly_zero( z() ).
//    This is to ensure that nearly_equal( Vector3D( a, b, c ), Vector3D( 0.0, 0.0, 0.0 ) ) and nearly_zero( Vector3D( a, b, c ) )
//    are consistent, and nearly_equal( Vector3D( a, b, c ), Vector3D( x, y, z ) ) cannot be simplified to
//    only testing if the norm()- or norm2()-values are nearly equal. It also has as a nice consequence that Vector3D and Vector2D behave the same.
//    It also avoids taking the square root and it makes nearly_zero( Vector2D ) and nearly_zero( Vector3D ) behave the same.
// Replace quaternions by rotors?
// abs() is C and only works for int "int abs( int )". std::abs() is C++ and in the C98 standard would always return at least float "double std::abs( int )", since C11, "int std::abs( int )" works.
// We probably need an include file with constants, such as TOLERANCE and CONSTANT_PI
// std:vector< double > e.g. are_all_nearly_equal( const double tolerance ); calculates average and checks if each value is less than tolerance from average.
// We currently never check if a determinant for a transformation is negative. Should allow this but should give warning. I think this is wrong:
//    we check in CrystalLattice.transform(), and that is always used, and there a warning is written out.
//    but we currently throw in add_centring_to_space_group_after_transformation().
// The whole Wavelength class must be reprogrammed to store an enum for synchrotron or one of the five lab sources and with or without monochromator.
// In ListOfDoubles: boundary checking for access.
// In ListOfDoubles: implement +=.
// Eulerian_angles() should use the EulerianAngles class.
// Add an RMSCD that uses the original atom labels to match up atoms (currently we either rely on the order or do a full matching from scratch).
// Why do we ever write anything to a file? Why not always use std::vector< std::string > and add a function that writes a std::vector< std::string > to a file?
// I probably need a class SettingsFile.

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

// Turn BUFBEM11_01 and BUFBEM11_02 etc. into BUFBEM11
std::string disordered_identifier_to_refcode( const std::string & input )
{
    if ( input.length() == 9 )
        return input.substr( 0, 6 );
    if ( input.length() == 11 )
        return input.substr( 0, 8 );
    return input;
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

    try // Generate a powder diffraction pattern from SHELX .hkl file.
    {
        // We need a CrystalLattice and a space group.
        std::string directory( "Experimental" );
        FileName input_file_name( directory, "GF", "cif" );
        CrystalStructure crystal_structure;
        std::cout << "Now reading cif... " + input_file_name.full_name() << std::endl;
        read_cif( input_file_name, crystal_structure );
        crystal_structure.apply_space_group_symmetry();
        PointGroup Laue_class = crystal_structure.space_group().Laue_class();
        ReflectionList reflection_list;
        reflection_list.read_hkl( FileName( "GP.hkl" ) );
        std::cout << ".hkl file has been read." << std::endl;
        ReflectionList reflection_list_final;
        std::vector< bool > done( reflection_list.size(), false );
        for ( size_t i( 0 ); i != reflection_list.size(); ++i )
        {
            if ( done[i] )
                continue;
            MillerIndices miller_indices = reflection_list.miller_indices( i );

// @@ This really should be a stand-alone function.

// Would it not be much easier to add member functions to ReflectionsList that
// first convert all MillerIndices to the representative of the equivalent MillerIndices
// and then to merge all those with the same MillerIndices?

            // Calculate equivalent reflections.
            std::set< MillerIndices > equivalent_reflections;
            equivalent_reflections.insert( miller_indices );
            for ( size_t j( 0 ); j != Laue_class.nsymmetry_operators(); ++j )
                equivalent_reflections.insert( miller_indices * Laue_class.symmetry_operator( j ) );
            double F_squared = 0.0;
            size_t number_of_equivalent_reflections_that_have_actually_been_measured( 0 );
            // It is possible that the same reflection has been measured multiple times, so we cannot simply find "a" reflection based on
            // a set of hkl indices, as the list may contain the same set of hkl indices multiple times. Therefore, the only way is to go through the
            // entire list from start to finish every single time.
            for ( size_t j( i+1 ); j != reflection_list.size(); ++j )
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
        reflection_list_final.save( FileName( directory, "GF_cal", "hkl" ) );
        PowderPatternCalculator powder_pattern_calculator( crystal_structure );
        powder_pattern_calculator.set_FWHM( 0.05 );
        powder_pattern_calculator.set_two_theta_end( Angle::from_degrees( 40.0 ) );
        powder_pattern_calculator.set_two_theta_step( Angle::from_degrees( 0.01 ) );
        PowderPattern powder_pattern;
        powder_pattern_calculator.calculate( reflection_list_final, powder_pattern );
        powder_pattern.save_xye( FileName( "GP_new.xye" ), true );
    MACRO_END_GAME

    // Transformation of the crystal structure (unit cell + atomic coordinates including ADPs + space group)
    // followed by a transformation of the atomic coordinates including ADPs.
    try
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        Vector3D com = crystal_structure.centre_of_mass( false );
        std::cout << "Centre of mass = " << std::endl;
        com.show();
        Matrix3D tranformation_matrix(  1.0,  0.0,  0.0,
                                        0.0,  1.0,  0.0,
                                        0.0,  0.0,  1.0 );
        if ( ! tranformation_matrix.is_nearly_the_identity() )
        {
            crystal_structure.transform( tranformation_matrix );
            std::cout << "Inverse transformation matrix:" << std::endl;
            std::cout << inverse( tranformation_matrix ) << std::endl;
        }
        if ( (false) )
        {
            SymmetryOperator symmetry_operator( "x,y+1/4,z" );
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
        bool there_is_a_rotation = ( ! rotation.is_nearly_the_identity() );
        bool there_is_a_translation = ( ! shift.nearly_zero() );
        if ( there_is_a_rotation || there_is_a_translation )
        {
            for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
            {
                Atom new_atom( crystal_structure.atom( i ) );
                if ( there_is_a_rotation )
                {
                    new_atom.set_position( ( rotation * crystal_structure.atom( i ).position() ) + shift );
                    if ( new_atom.ADPs_type() == Atom::ANISOTROPIC )
                        new_atom.set_anisotropic_displacement_parameters( rotate_adps( new_atom.anisotropic_displacement_parameters(), rotation, crystal_structure.crystal_lattice() ) );
                }
                else if ( there_is_a_translation )
                    new_atom.set_position( ( crystal_structure.atom( i ).position() ) + shift );
                crystal_structure.set_atom( i, new_atom );
            }
        }
        // In Mercury, if the space-group name and the set of symmetry operators do not match up,
        // the space-group name takes precedence, so we have to erase it to ensure that the
        // symmetry operators are used instead.
        if ( (false) )
        {
            SpaceGroup space_group = crystal_structure.space_group();
            space_group.set_name( "" );
            // If the determinant is 1/2, 1/3, 1/4 etc. then we are probably
            // reducing a centred unit cell to primitive, so the number of symmetry
            // operators is also reduced, namely by a factor of 2, 3, 4 etc.
            // the determinant is either 1 or <= 1/2, so I use 3/4 as a value to
            // circumvent rounding errors 
            if ( tranformation_matrix.determinant() < 0.75 )
                space_group.remove_duplicate_symmetry_operators();
            if ( tranformation_matrix.determinant() > 1.5 )
                add_centring_to_space_group_after_transformation( tranformation_matrix, space_group );
            crystal_structure.set_space_group( space_group );
        }
        crystal_structure.save_cif( replace_extension( append_to_file_name( input_file_name, "_transformed" ), "cif" ) );
    MACRO_END_GAME

    try // Simulate a realistic experimental powder diffraction pattern.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        crystal_structure.apply_space_group_symmetry();
        bool save_all_for_debugging( true );
        RealisticXRPDSimulatorSettings settings;
        settings.set_wavelength( Wavelength() );
        settings.set_two_theta_start( Angle( 2.0, Angle::DEGREES ) );
        settings.set_two_theta_end( Angle( 40.0, Angle::DEGREES ) );
        settings.set_two_theta_step( Angle( 0.015, Angle::DEGREES ) );
        settings.set_FWHM( 0.1 );
        settings.set_zero_point_error( Angle( 0.02, Angle::DEGREES ) );
        if ( true )
            settings.set_preferred_orientation( select_realistic_preferred_orientation_direction( crystal_structure.crystal_lattice() ), 0.9 );
        if ( true )
            settings.set_finger_cox_jephcoat( 10.0 / 400.0, 10.0 / 400.0 ); // Finger-Cox-Jephcoat
        settings.set_include_background( true );
        settings.set_include_noise( true );
        settings.set_include_noise_for_zero_background( 20 );
        settings.set_Bragg_total_signal_normalisation( 10000.0 );
        settings.set_background_total_signal_normalisation( 0.2*10000.0 );
        settings.set_highest_peak( 10000.0 );
        RealisticXRPDSimulator realistic_XRPD_simulator( crystal_structure, settings );
        PowderPattern powder_pattern = realistic_XRPD_simulator.calculate();
        if ( save_all_for_debugging )
            realistic_XRPD_simulator.Bragg_diffraction().save_xye( replace_extension( append_to_file_name( input_file_name, "_Bragg" ), "xye" ), true );
        if ( save_all_for_debugging )
            realistic_XRPD_simulator.background().save_xye( replace_extension( append_to_file_name( input_file_name, "_background" ), "xye" ), true );
        if ( save_all_for_debugging )
            realistic_XRPD_simulator.noise().save_xye( replace_extension( append_to_file_name( input_file_name, "_noise" ), "xye" ), true );
        // A second one for the NaCl.
        if ( false )
        {
            CrystalStructure crystal_structure_NaCl = NaCl();
            crystal_structure_NaCl.apply_space_group_symmetry();
            RealisticXRPDSimulatorSettings settings_NaCl( settings );
            settings_NaCl.set_FWHM( 0.1 );
            settings_NaCl.unset_preferred_orientation();
            settings_NaCl.set_include_background( false );
            settings_NaCl.set_Bragg_total_signal_normalisation( 0.02 * settings.Bragg_total_signal_normalisation() );
            RealisticXRPDSimulator realistic_XRPD_simulator_NaCl( crystal_structure_NaCl, settings_NaCl );
            PowderPattern powder_pattern_NaCl = realistic_XRPD_simulator_NaCl.calculate();
            powder_pattern_NaCl.scale( realistic_XRPD_simulator.scale_factor() / realistic_XRPD_simulator_NaCl.scale_factor() );
            powder_pattern_NaCl.make_counts_integer();
            powder_pattern_NaCl.recalculate_estimated_standard_deviations();
            if ( save_all_for_debugging )
                powder_pattern_NaCl.save_xye( replace_extension( append_to_file_name( input_file_name, "_NaCl" ), "xye" ), true );
            powder_pattern += powder_pattern_NaCl;
        }
        powder_pattern.save_xye( replace_extension( input_file_name, "xye" ), true );
        settings.save( replace_extension( input_file_name, "txt" ) );
        if ( false )
        {
            PowderPattern estimated_background = calculate_Brueckner_background( powder_pattern,
                                                                                 50, // niterations
                                                                                 round_to_int( 50.0 * ( Angle::from_degrees( 0.015 ) / powder_pattern.average_two_theta_step() ) ), // window
                                                                                 true, // apply_smoothing
                                                                                 5 ); // smoothing_window
            if ( save_all_for_debugging )
                estimated_background.save_xye( replace_extension( append_to_file_name( input_file_name, "_Brueckner_BKGR" ), "xye" ), true );
            powder_pattern -= estimated_background;
            powder_pattern.save_xye( replace_extension( append_to_file_name( input_file_name, "_Brueckner_BKGR_subtracted" ), "xye" ), true );
        }
    MACRO_END_GAME

    try // Take all disordered atoms that have been modelled as large ADPs and change them into a split-atom model.
    {
        if ( argc < 3 )
            throw std::runtime_error( "Please give the name of a .cif file, a threshold (ca. 0.1 / 0.125) and a factor (ca. 0.75)." );
        FileName input_file_name( argv[ 1 ] );
        double threshold = 0.1;
        if ( argc > 2 )
            threshold = string2double( argv[ 2 ] );
        std::cout << "Threshold = " << threshold << " (good values are 0.1 to 0.125)." << std::endl;
        double factor = 0.75;
        if ( argc == 4 )
            factor = string2double( argv[ 3 ] );
        std::cout << "Factor = " << factor << " (a good value is 0.75)." << std::endl;
        CrystalStructure crystal_structure;
        std::cout << "Now reading cif... " + input_file_name.full_name() << std::endl;
        read_cif( input_file_name, crystal_structure );
        for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
        {
            Atom atom = crystal_structure.atom( i );
            if ( atom.ADPs_type() != Atom::ANISOTROPIC )
                continue;
            SymmetricMatrix3D Ucart = atom.anisotropic_displacement_parameters().U_cart();
            std::vector< double > eigenvalues;
            std::vector< NormalisedVector3D > eigenvectors;
            calculate_eigenvalues( Ucart, eigenvalues, eigenvectors );
            if ( eigenvalues[2] < threshold )
                continue;
            // It is clear that some random scaling factor must be involved, because the size of the ADP depends on the probability level.
            // Empirically, 1.58 corresponds to 50%, which is the default in Mercury. It is probably something like pi/sqrt(4), but I have
            // not been able to find the exact value.
            // It may be 1.5382, which is the C in Table 6.1 in CCDC/develop_related/chap6.pdf .
            // But we do not want to be on the outer edge of the ADP, we want the two new atoms to both be within the ADP. Then a value of 0.75 is much better.
            Vector3D new_position_1 = crystal_structure.crystal_lattice().fractional_to_orthogonal( atom.position() ) + factor * sqrt( eigenvalues[2] ) * eigenvectors[2];
            new_position_1 = crystal_structure.crystal_lattice().orthogonal_to_fractional( new_position_1 );
            Atom new_atom_1( atom.element(), new_position_1, atom.label() + "a" );
            new_atom_1.set_occupancy( 0.5 );
            // Uiso = U_cart().trace() / 3.0;
            // U_cart().trace() = eigenvalues[2] + eigenvalues[1] + eigenvalues[0]
            // By splitting an atom over two positions, we say that only half of the largest principal axis is to be assigned to this atom, the other half should be assigned to the other atom.
            // So the U_cart().trace() for this atom should be ( 0.5 * eigenvalues[2] ) + eigenvalues[1] + eigenvalues[0]
            // double new_Uiso = ( ( 0.5 * eigenvalues[2] ) + eigenvalues[1] + eigenvalues[0] ) / 3.0;
            // I did not like these Uisos, they were too big. The volume of an ellipsoid is k * a * b * c where a, b, and c are the principal axes. For a sphere, it is k * r^3.
            // So just calculate the volume, and calculate the radius that would give a sphere with half the volume.
            double new_Uiso = std::pow( 0.5 * eigenvalues[2] * eigenvalues[1] * eigenvalues[0], 1.0/3.0 );
            new_atom_1.set_Uiso( new_Uiso );
            crystal_structure.add_atom( new_atom_1 );
            Vector3D new_position_2 = crystal_structure.crystal_lattice().fractional_to_orthogonal( atom.position() ) - factor * sqrt( eigenvalues[2] ) * eigenvectors[2];
            new_position_2 = crystal_structure.crystal_lattice().orthogonal_to_fractional( new_position_2 );
            Atom new_atom_2( atom.element(), new_position_2, atom.label() + "b" );
            new_atom_2.set_occupancy( 0.5 );
            new_atom_2.set_Uiso( new_Uiso );
            crystal_structure.add_atom( new_atom_2 );
        }
        crystal_structure.save_cif( append_to_file_name( input_file_name, "_split" ) );
    MACRO_END_GAME

    try // Supercell with and without original space group.
    {
        if ( argc != 5 )
            throw std::runtime_error( "Usage: u v w <name>.cif." );
        size_t u = string2integer( argv[ 1 ] );
        size_t v = string2integer( argv[ 2 ] );
        size_t w = string2integer( argv[ 3 ] );
        FileName input_file_name( argv[ 4 ] );
        CrystalStructure crystal_structure;
        read_cif( input_file_name, crystal_structure );
        CrystalStructure crystal_structure_2( crystal_structure );
        crystal_structure_2.supercell( u, v, w );
        crystal_structure_2.save_cif( append_to_file_name( input_file_name, "_" + size_t2string( u ) + "_" + size_t2string( v ) + "_" + size_t2string( w ) ) );
        Matrix3D transformation_matrix(    u,  0.0,  0.0,
                                        0.0,    v,  0.0,
                                        0.0,  0.0,    w );
        crystal_structure.transform( transformation_matrix );
        // In Mercury, if the space-group name and the set of symmetry operators do not match up,
        // the space-group name takes precedence, so we have to erase it to ensure that the
        // symmetry operators are used instead.
        SpaceGroup space_group = crystal_structure.space_group();
        space_group.set_name( "" );
        add_centring_to_space_group_after_transformation( transformation_matrix, space_group );
        crystal_structure.set_space_group( space_group );
        crystal_structure.save_cif( append_to_file_name( input_file_name, "_" + size_t2string( u ) + "_" + size_t2string( v ) + "_" + size_t2string( w ) + "_SpGr" ) );
    MACRO_END_GAME
        
    try // CrystalStructure::convert_to_P1().
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        crystal_structure.convert_to_P1();
        crystal_structure.make_atom_labels_unique();
        crystal_structure.save_cif( append_to_file_name( input_file_name, "_P1" ) );
    MACRO_END_GAME

    try // Write .inp from .cif + .xye file.
    {
        std::vector< std::string > extensions;
        extensions.push_back( "cif" );
        extensions.push_back( "xye" );
        std::vector< FileName > input_file_names = sort_file_names_by_extension( argc, argv, extensions );
        inp_writer_distance_restraints( input_file_names[0], input_file_names[1] );
    MACRO_END_GAME

    try // Write .inp from .cif + two _restraints.txt files + .xye file.
    {
        std::vector< std::string > extensions;
        extensions.push_back( "cif" );
        extensions.push_back( "xye" );
        std::vector< FileName > input_file_names = sort_file_names_by_extension( argc, argv, extensions );
        inp_writer( input_file_names[0], input_file_names[1] );
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
        generate_powder_cif.generate_R_input_file( GeneratePowderCIF::ALWAYS_ZOOM );
    MACRO_END_GAME

    try // Rotate group.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        if ( true )
        {
        std::string atom_1_label( "C31" );
        std::string atom_2_label( "C32" );
        Vector3D C1 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_1_label ) ).position() );
        Vector3D C2 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_2_label ) ).position() );
        NormalisedVector3D n = normalised_vector( C1 - C2 );
        std::vector< std::string > labels_to_be_rotated;
        labels_to_be_rotated.push_back( "C33a" );
        Angle additional_angle = Angle::from_degrees( -12.0 );
        for ( size_t i( 0 ); i != labels_to_be_rotated.size(); ++i )
        {
            std::cout << labels_to_be_rotated[i] << std::endl;
            Atom new_atom( crystal_structure.atom( crystal_structure.find_label( labels_to_be_rotated[i] ) ) );
            Vector3D tVector = crystal_structure.crystal_lattice().fractional_to_orthogonal( new_atom.position() );
            tVector = rotate_point_about_axis( tVector, C1, n, additional_angle );
            new_atom.set_position( crystal_structure.crystal_lattice().orthogonal_to_fractional( tVector ) );
        //    new_atom.set_label( labels_to_be_rotated[i] + "a" );
        //    crystal_structure.add_atom( new_atom ); // Copy, then rotate

       //     new_atom = Atom( crystal_structure.atom( crystal_structure.find_label( labels_to_be_rotated[i] ) ) );
       //     tVector = crystal_structure.crystal_lattice().fractional_to_orthogonal( new_atom.position() );
       //     tVector = rotate_point_about_axis( tVector, C1, n, -additional_angle );
       //     new_atom.set_position( crystal_structure.crystal_lattice().orthogonal_to_fractional( tVector ) );
       //     new_atom.set_label( labels_to_be_rotated[i] + "b" );
            crystal_structure.set_atom( crystal_structure.find_label( labels_to_be_rotated[i] ), new_atom );
        }
        }

        if ( false )
        {
        std::string atom_1_label( "N9" );
        std::string atom_2_label( "C213" );
        Vector3D C1 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_1_label ) ).position() );
        Vector3D C2 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_2_label ) ).position() );
        NormalisedVector3D n = normalised_vector( C1 - C2 );
        std::vector< std::string > labels_to_be_rotated;
        labels_to_be_rotated.push_back( "O29" );
        labels_to_be_rotated.push_back( "N13" );
        labels_to_be_rotated.push_back( "H165" );
        labels_to_be_rotated.push_back( "H161" );
        Angle additional_angle = Angle::from_degrees( 180.0 );
        for ( size_t i( 0 ); i != labels_to_be_rotated.size(); ++i )
        {
            std::cout << labels_to_be_rotated[i] << std::endl;
            Atom new_atom( crystal_structure.atom( crystal_structure.find_label( labels_to_be_rotated[i] ) ) );
            Vector3D tVector = crystal_structure.crystal_lattice().fractional_to_orthogonal( new_atom.position() );
            tVector = rotate_point_about_axis( tVector, C1, n, additional_angle );
            new_atom.set_position( crystal_structure.crystal_lattice().orthogonal_to_fractional( tVector ) );
        //    new_atom.set_label( labels_to_be_rotated[i] + "a" );
        //    crystal_structure.add_atom( new_atom ); // Copy, then rotate

       //     new_atom = Atom( crystal_structure.atom( crystal_structure.find_label( labels_to_be_rotated[i] ) ) );
       //     tVector = crystal_structure.crystal_lattice().fractional_to_orthogonal( new_atom.position() );
       //     tVector = rotate_point_about_axis( tVector, C1, n, -additional_angle );
       //     new_atom.set_position( crystal_structure.crystal_lattice().orthogonal_to_fractional( tVector ) );
       //     new_atom.set_label( labels_to_be_rotated[i] + "b" );
            crystal_structure.set_atom( crystal_structure.find_label( labels_to_be_rotated[i] ), new_atom );
        }
        }

        if ( false )
        {
        std::string atom_1_label( "C117" );
        std::string atom_2_label( "C121" );
        Vector3D C1 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_1_label ) ).position() );
        Vector3D C2 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_2_label ) ).position() );
        NormalisedVector3D n = normalised_vector( C1 - C2 );
        std::vector< std::string > labels_to_be_rotated;
        labels_to_be_rotated.push_back( "O17" );
        labels_to_be_rotated.push_back( "O21" );
        labels_to_be_rotated.push_back( "H105" );
        Angle additional_angle = Angle::from_degrees( 180.0 );
        for ( size_t i( 0 ); i != labels_to_be_rotated.size(); ++i )
        {
            std::cout << labels_to_be_rotated[i] << std::endl;
            Atom new_atom( crystal_structure.atom( crystal_structure.find_label( labels_to_be_rotated[i] ) ) );
            Vector3D tVector = crystal_structure.crystal_lattice().fractional_to_orthogonal( new_atom.position() );
            tVector = rotate_point_about_axis( tVector, C1, n, additional_angle );
            new_atom.set_position( crystal_structure.crystal_lattice().orthogonal_to_fractional( tVector ) );
        //    new_atom.set_label( labels_to_be_rotated[i] + "a" );
        //    crystal_structure.add_atom( new_atom ); // Copy, then rotate

       //     new_atom = Atom( crystal_structure.atom( crystal_structure.find_label( labels_to_be_rotated[i] ) ) );
       //     tVector = crystal_structure.crystal_lattice().fractional_to_orthogonal( new_atom.position() );
       //     tVector = rotate_point_about_axis( tVector, C1, n, -additional_angle );
       //     new_atom.set_position( crystal_structure.crystal_lattice().orthogonal_to_fractional( tVector ) );
       //     new_atom.set_label( labels_to_be_rotated[i] + "b" );
            crystal_structure.set_atom( crystal_structure.find_label( labels_to_be_rotated[i] ), new_atom );
        }
        }

        if ( false )
        {
        std::string atom_1_label( "C61" );
        std::string atom_2_label( "C65" );
        Vector3D C1 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_1_label ) ).position() );
        Vector3D C2 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( crystal_structure.find_label( atom_2_label ) ).position() );
        NormalisedVector3D n = normalised_vector( C1 - C2 );
        std::vector< std::string > labels_to_be_rotated;
        labels_to_be_rotated.push_back( "O5" );
        labels_to_be_rotated.push_back( "O9" );
        labels_to_be_rotated.push_back( "H49" );
        Angle additional_angle = Angle::from_degrees( 180.0 );
        for ( size_t i( 0 ); i != labels_to_be_rotated.size(); ++i )
        {
            std::cout << labels_to_be_rotated[i] << std::endl;
            Atom new_atom( crystal_structure.atom( crystal_structure.find_label( labels_to_be_rotated[i] ) ) );
            Vector3D tVector = crystal_structure.crystal_lattice().fractional_to_orthogonal( new_atom.position() );
            tVector = rotate_point_about_axis( tVector, C1, n, additional_angle );
            new_atom.set_position( crystal_structure.crystal_lattice().orthogonal_to_fractional( tVector ) );
        //    new_atom.set_label( labels_to_be_rotated[i] + "a" );
        //    crystal_structure.add_atom( new_atom ); // Copy, then rotate

       //     new_atom = Atom( crystal_structure.atom( crystal_structure.find_label( labels_to_be_rotated[i] ) ) );
       //     tVector = crystal_structure.crystal_lattice().fractional_to_orthogonal( new_atom.position() );
       //     tVector = rotate_point_about_axis( tVector, C1, n, -additional_angle );
       //     new_atom.set_position( crystal_structure.crystal_lattice().orthogonal_to_fractional( tVector ) );
       //     new_atom.set_label( labels_to_be_rotated[i] + "b" );
            crystal_structure.set_atom( crystal_structure.find_label( labels_to_be_rotated[i] ), new_atom );
        }
        }

        crystal_structure.save_cif( append_to_file_name( input_file_name, "_rotated" ) );
    MACRO_END_GAME

    try // Write TOPAS input file for a temperature series.
    {
        std::vector< std::string > extensions;
        extensions.push_back( "cif" );
        extensions.push_back( "xye" );
        double temperature_start    = 120.0;
        double temperature_end      = 340.0;
        double temperature_interval = 20.0;
        size_t ntemperatures = 1 + ( ( temperature_end - temperature_start ) / temperature_interval );
        if ( std::abs( temperature_end - ( temperature_start + ntemperatures * temperature_interval ) ) > 0.5 )
            std::cout << "Warning: temperature_start, temperature_interval and temperature_end are inconsistent." << std::endl;
        std::vector< FileName > input_file_names = sort_file_names_by_extension( argc, argv, extensions );
        CrystalStructure crystal_structure;
        read_cif( input_file_names[0], crystal_structure );
        PowderPattern powder_pattern;
        powder_pattern.read_xye( input_file_names[1] );
        TextFileWriter text_file_writer( replace_extension( input_file_names[1], "inp" ) );
        write_preamble( text_file_writer );
        size_t nbackground_terms = 20;
        for ( size_t i( 0 ); i != nbackground_terms; ++i )
            text_file_writer.write_line( "prm JvdS_BKGR_" + size_t2string( i+1, 2 ) + " 0.0" );
        text_file_writer.write_line( "prm JvdS_zero_point_error 0.02" );
        text_file_writer.write_line( "prm JvdS_filament_length 10.0" );
        text_file_writer.write_line( "prm JvdS_sample_length 2.0" );
        text_file_writer.write_line( "prm JvdS_receiving_slit_length 4.0" );
        text_file_writer.write_line( "prm JvdS_CS_G 500.0" );
        text_file_writer.write_line( "prm JvdS_CS_L 500.0" );
        text_file_writer.write_line( "prm JvdS_Strain_G 1.0" );
        text_file_writer.write_line( "prm JvdS_Strain_L 1.0" );
        write_unit_cell_variables( text_file_writer, crystal_structure.crystal_lattice() );
        text_file_writer.write_line( "prm JvdS_Zero_Error_0 0.03" );
        text_file_writer.write_line( "prm JvdS_Zero_Error_1_0 0.00000" );
        text_file_writer.write_line( "prm JvdS_Zero_Error_1 = JvdS_Zero_Error_1_0 / 100000.0;" );
        text_file_writer.write_line( "prm JvdS_bnonh_0 3.0" );
        text_file_writer.write_line( "prm JvdS_bnonh_1_0 0.10" );
        text_file_writer.write_line( "prm JvdS_bnonh_1 = JvdS_bnonh_1_0 / 1000.0;" );
        text_file_writer.write_line( "prm JvdS_scale 1" );
        text_file_writer.write_line( "prm JvdS_PO 1.0" );
        for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
        {
            text_file_writer.write_line( "prm " + crystal_structure.atom( i ).label() + "_x " + double2string_pad_plus( crystal_structure.atom( i ).position().x(), 5, ' ' ) );
            text_file_writer.write_line( "prm " + crystal_structure.atom( i ).label() + "_y " + double2string_pad_plus( crystal_structure.atom( i ).position().y(), 5, ' ' ) );
            text_file_writer.write_line( "prm " + crystal_structure.atom( i ).label() + "_z " + double2string_pad_plus( crystal_structure.atom( i ).position().z(), 5, ' ' ) );
        }
        for ( size_t i( 0 ); i != ntemperatures; ++i )
        {
            double temperature = temperature_start + i * temperature_interval;
            text_file_writer.write_line( "xdd " + input_file_names[1].name() + "_" + size_t2string( i+1, 2 ) + ".xye xye_format" );
            text_file_writer.write_line( "  local !JvdS_temperature " + double2string( temperature ) );
            text_file_writer.write_line( "  bkg" );
            for ( size_t j( 0 ); j != nbackground_terms; ++j )
                text_file_writer.write_line( "  =JvdS_BKGR_" + size_t2string( j+1, 2 ) + ";" );
            text_file_writer.write_line( "  start_X       " + double2string( powder_pattern.two_theta_start().value_in_degrees() ) );
            text_file_writer.write_line( "  finish_X      " + double2string( powder_pattern.two_theta_end().value_in_degrees() ) );
            text_file_writer.write_line( "  x_calculation_step " + double2string( powder_pattern.average_two_theta_step().value_in_degrees() ) );
            text_file_writer.write_line( "  Zero_Error(@ , JvdS_Zero_Error_0 + JvdS_Zero_Error_1 * JvdS_temperature; : 0.0" );
            text_file_writer.write_line( "  LP_Factor( 26.5 )" );
            text_file_writer.write_line( "  axial_conv" );
            text_file_writer.write_line( "    filament_length =JvdS_filament_length;" );
            text_file_writer.write_line( "    sample_length =JvdS_sample_length;" );
            text_file_writer.write_line( "    receiving_slit_length =JvdS_receiving_slit_length;" );
            text_file_writer.write_line( "    axial_n_beta 50" );
            text_file_writer.write_line( "  lam" );
            text_file_writer.write_line( "    ymin_on_ymax 0.001" );
            text_file_writer.write_line( "    la 1 lo 1.540560" );
            text_file_writer.write_line( "  str" );
            text_file_writer.write_line( "    r_bragg 0.0" );
            text_file_writer.write_line( "    CS_G(@ , =JvdS_CS_G; )" );
            text_file_writer.write_line( "    CS_L(@ , =JvdS_CS_L; )" );
            text_file_writer.write_line( "    Strain_G(@ , =JvdS_Strain_G; )" );
            text_file_writer.write_line( "    Strain_L(@ , =JvdS_Strain_L; )" );
//    text_file_writer.write_line( "    prm  sh_scale_l" + aal + " 0.01" );
//    text_file_writer.write_line( "    spherical_harmonics_hkl sh_l" + aal );
//    text_file_writer.write_line( "      sh_order 6" );
//    text_file_writer.write_line( "    lor_fwhm = Abs( sh_scale_l" + aal + " * sh_l" + aal + " );" );
//    text_file_writer.write_line( "    prm  sh_scale_g" + aal + "  0.01" );
//    text_file_writer.write_line( "    spherical_harmonics_hkl sh_g" + aal );
//    text_file_writer.write_line( "      sh_order 6" );
//    text_file_writer.write_line( "    gauss_fwhm = Abs( sh_scale_g" + aal + " * sh_g" + aal + " );" );
            write_unit_cell( text_file_writer, crystal_structure.crystal_lattice(), true, true );
            text_file_writer.write_line( "    MVW( 0.0, 0.0, 0.0 )");
            text_file_writer.write_line( "    space_group \"" + remove( remove( crystal_structure.space_group().name(), '_' ), ' ') + "\"");
            text_file_writer.write_line( "    scale = JvdS_scale/100000;" );
//            text_file_writer.write_line( "'    PO_Spherical_Harmonics( sh, 6 )" );
            text_file_writer.write_line( "    PO(@ , =JvdS_PO;, , 0 0 1 )" );
            text_file_writer.write_line( "    local bnonh = JvdS_bnonh_0 + JvdS_bnonh_1 * JvdS_temperature; : 0.0" );
            text_file_writer.write_line( "    local bh = 1.2 * bnonh; : 0.0" );
            for ( size_t j( 0 ); j != crystal_structure.natoms(); ++j )
            {
                text_file_writer.write( "    site " + crystal_structure.atom( j ).label() + " x =" + crystal_structure.atom( j ).label() + "_x" + ";" +
                                                                                            " y =" + crystal_structure.atom( j ).label() + "_y" + ";" +
                                                                                            " z =" + crystal_structure.atom( j ).label() + "_z" + ";" +
                                             " occ " + pad( crystal_structure.atom( j ).element().symbol(), 2, ' ' ) + " 1 beq = " );
                if ( crystal_structure.atom( j ).element().is_H_or_D() )
                    text_file_writer.write_line( "bh;" );
                else
                    text_file_writer.write_line( "bnonh;" );
            }
            if ( i == ntemperatures-1 )
                text_file_writer.write_line( "    Out_CIF_STR( " + input_file_names[0].name() + "_RR.cif" + " )" );
        }
//        text_file_writer.write_line( "    prm !flatten_width 0" );
//        text_file_writer.write_line( "    prm !flatten_weight    100000" );
//        text_file_writer.write_line( "    'Flatten( C19 N11 N35 C47 H54 C50 H58 C45 H49 C34 C46, , 0.0, flatten_width, flatten_weight )" );
//        text_file_writer.write_line( "  xdd_out " + FileName( input_cif_file_name.directory(), input_cif_file_name.name() + "_profile", "txt" ).full_name() + " load out_record out_fmt out_eqn" );
//        text_file_writer.write_line( "  {" );
//        text_file_writer.write_line( "      \" %11.5f \" = X;" );
//        text_file_writer.write_line( "      \" %11.5f \" = Yobs;" );
//        text_file_writer.write_line( "      \" %11.5f \" = Ycalc;" );
//        text_file_writer.write_line( "      \" %11.5f\\n\" = SigmaYobs;" );
//        text_file_writer.write_line( "  }" );
//        text_file_writer.write_line( "  phase_out " + FileName( input_cif_file_name.directory(), input_cif_file_name.name() + "_tickmarks", "txt" ).full_name() + " load out_record out_fmt out_eqn" );
//        text_file_writer.write_line( "  {" );
//        text_file_writer.write_line( "      \" %11.5f -200\\n\" = 2.0 * Rad * Th;" );
//        text_file_writer.write_line( "  }" );
    MACRO_END_GAME

    try // Split powder pattern into multiple patterns.
    {
        if ( argc != 4 )
            throw std::runtime_error( "Please give the name of a .xye file and the number of patterns it has to be split into and a .cif file." );
        FileName input_file_name( argv[ 1 ] );
        PowderPattern powder_pattern;
        powder_pattern.read_xye( input_file_name );
        size_t n = string2integer( argv[ 2 ] );
        FileName input_cif_file_name( argv[ 3 ] );
        CrystalStructure crystal_structure;
        read_cif( input_cif_file_name, crystal_structure );
        std::vector< PowderPattern > powder_patterns = split( powder_pattern, n, false );
        for ( size_t i( 0 ); i != n; ++i )
            powder_patterns[i].save_xye( append_to_file_name( input_file_name, "_" + size_t2string( i+1, 2 ) ), false );
            
        TextFileWriter text_file_writer( replace_extension( input_cif_file_name, "inp" ) );
        write_preamble( text_file_writer );

        size_t nbackground_terms = 20;
        for ( size_t i( 0 ); i != nbackground_terms; ++i )
            text_file_writer.write_line( "prm JvdS_BKGR_" + size_t2string( i+1, 2 ) + " 0.0" );

        text_file_writer.write_line( "prm JvdS_zero_point_error 0.0" );
        text_file_writer.write_line( "prm JvdS_filament_length 10.0" );
        text_file_writer.write_line( "prm JvdS_sample_length 2.0" );
        text_file_writer.write_line( "prm JvdS_receiving_slit_length 4.0" );
        text_file_writer.write_line( "prm JvdS_CS_G 500.0" );
        text_file_writer.write_line( "prm JvdS_CS_L 500.0" );
        text_file_writer.write_line( "prm JvdS_Strain_G 1.0" );
        text_file_writer.write_line( "prm JvdS_Strain_L 1.0" );
        write_unit_cell_variables( text_file_writer, crystal_structure.crystal_lattice() );
        text_file_writer.write_line( "prm bnonh 3.0" );
        text_file_writer.write_line( "prm bh = 1.2 * bnonh;" );
        text_file_writer.write_line( "prm JvdS_scale 1" );
        text_file_writer.write_line( "prm JvdS_PO 1.0" );
        for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
        {
            text_file_writer.write_line( "prm " + crystal_structure.atom( i ).label() + "x " + double2string_pad_plus( crystal_structure.atom( i ).position().x(), 5, ' ' ) );
            text_file_writer.write_line( "prm " + crystal_structure.atom( i ).label() + "y " + double2string_pad_plus( crystal_structure.atom( i ).position().y(), 5, ' ' ) );
            text_file_writer.write_line( "prm " + crystal_structure.atom( i ).label() + "z " + double2string_pad_plus( crystal_structure.atom( i ).position().z(), 5, ' ' ) );
        }
        size_t padding_length = 4;
        for ( size_t i( 0 ); i != n; ++i )
        {
            text_file_writer.write_line( "xdd " + append_to_file_name( input_file_name, "_" + size_t2string( i+1, padding_length ) ).file_name() + " xye_format" );
            text_file_writer.write_line( "  bkg" );
            for ( size_t j( 0 ); j != nbackground_terms; ++j )
                text_file_writer.write_line( "  =JvdS_BKGR_" + size_t2string( j+1, 2 ) + ";" );
            text_file_writer.write_line( "  start_X       " + double2string( powder_pattern.two_theta_start().value_in_degrees() ) );
            text_file_writer.write_line( "  finish_X      " + double2string( powder_pattern.two_theta_end().value_in_degrees() ) );
            text_file_writer.write_line( "  x_calculation_step " + double2string( powder_pattern.average_two_theta_step().value_in_degrees() ) );
            text_file_writer.write_line( "  Zero_Error(@ , =JvdS_zero_point_error; )" );
            text_file_writer.write_line( "  LP_Factor( 26.5 )" );
            text_file_writer.write_line( "  axial_conv" );
            text_file_writer.write_line( "    filament_length =JvdS_filament_length;" );
            text_file_writer.write_line( "    sample_length =JvdS_sample_length;" );
            text_file_writer.write_line( "    receiving_slit_length =JvdS_receiving_slit_length;" );
            text_file_writer.write_line( "    axial_n_beta 50" );
            text_file_writer.write_line( "  lam" );
            text_file_writer.write_line( "    ymin_on_ymax 0.001" );
            text_file_writer.write_line( "    la 1 lo 1.540560" );
            text_file_writer.write_line( "  str" );
            text_file_writer.write_line( "    r_bragg 0.0" );
            text_file_writer.write_line( "    CS_G(@ , =JvdS_CS_G; )" );
            text_file_writer.write_line( "    CS_L(@ , =JvdS_CS_L; )" );
            text_file_writer.write_line( "    Strain_G(@ , =JvdS_Strain_G; )" );
            text_file_writer.write_line( "    Strain_L(@ , =JvdS_Strain_L; )" );
//    text_file_writer.write_line( "    prm  sh_scale_l" + aal + " 0.01" );
//    text_file_writer.write_line( "    spherical_harmonics_hkl sh_l" + aal );
//    text_file_writer.write_line( "      sh_order 6" );
//    text_file_writer.write_line( "    lor_fwhm = Abs( sh_scale_l" + aal + " * sh_l" + aal + " );" );
//    text_file_writer.write_line( "    prm  sh_scale_g" + aal + "  0.01" );
//    text_file_writer.write_line( "    spherical_harmonics_hkl sh_g" + aal );
//    text_file_writer.write_line( "      sh_order 6" );
//    text_file_writer.write_line( "    gauss_fwhm = Abs( sh_scale_g" + aal + " * sh_g" + aal + " );" );
            write_unit_cell( text_file_writer, crystal_structure.crystal_lattice() );
            text_file_writer.write_line( "    MVW( 0.0, 0.0, 0.0 )");
            text_file_writer.write_line( "    space_group \"" + remove( remove( crystal_structure.space_group().name(), '_' ), ' ') + "\"");
            text_file_writer.write_line( "    scale =JvdS_scale/100000;" );
//            text_file_writer.write_line( "'    PO_Spherical_Harmonics( sh, 6 )" );
            text_file_writer.write_line( "    PO(@ , =JvdS_PO;, , 0 0 1 )" );
            for ( size_t j( 0 ); j != crystal_structure.natoms(); ++j )
            {
                text_file_writer.write( "    site " + crystal_structure.atom( j ).label() + " x =" + crystal_structure.atom( j ).label() + "x" + ";" +
                                                                                            " y =" + crystal_structure.atom( j ).label() + "y" + ";" +
                                                                                            " z =" + crystal_structure.atom( j ).label() + "z" + ";" +
                                             " occ " + pad( crystal_structure.atom( j ).element().symbol(), 2, ' ' ) + " 1 beq = " );
                if ( crystal_structure.atom( j ).element().is_H_or_D() )
                    text_file_writer.write_line( "bh;" );
                else
                    text_file_writer.write_line( "bnonh;" );
            }
            if ( i == n-1 )
                text_file_writer.write_line( "    Out_CIF_STR( " + FileName( input_cif_file_name.directory(), input_cif_file_name.name() + "_RR", "cif" ).full_name() + " )" );
//    text_file_writer.write_line( "    prm !flatten_width 0" );
//    text_file_writer.write_line( "    prm !flatten_weight    100000" );
//    text_file_writer.write_line( "    'Flatten( C19 N11 N35 C47 H54 C50 H58 C45 H49 C34 C46, , 0.0, flatten_width, flatten_weight )" );
//    text_file_writer.write_line( "    Out_CIF_STR( " + FileName( input_cif_file_name.directory(), input_cif_file_name.name() + "_RR", "cif" ).full_name() + " )" );
//    text_file_writer.write_line( "  xdd_out " + FileName( input_cif_file_name.directory(), input_cif_file_name.name() + "_profile", "txt" ).full_name() + " load out_record out_fmt out_eqn" );
//    text_file_writer.write_line( "  {" );
//    text_file_writer.write_line( "      \" %11.5f \" = X;" );
//    text_file_writer.write_line( "      \" %11.5f \" = Yobs;" );
//    text_file_writer.write_line( "      \" %11.5f \" = Ycalc;" );
//    text_file_writer.write_line( "      \" %11.5f\\n\" = SigmaYobs;" );
//    text_file_writer.write_line( "  }" );
//    text_file_writer.write_line( "  phase_out " + FileName( input_cif_file_name.directory(), input_cif_file_name.name() + "_tickmarks", "txt" ).full_name() + " load out_record out_fmt out_eqn" );
//    text_file_writer.write_line( "  {" );
//    text_file_writer.write_line( "      \" %11.5f -200\\n\" = 2.0 * Rad * Th;" );
//    text_file_writer.write_line( "  }" );
        }
    MACRO_END_GAME

    try // // Write TOPAS input file for a phase transition captured by XRPD.
    {
// C:\Users\jacco\Documents\PowderPatterns\PhaseTransition\structure_000001.cif C:\Users\jacco\Documents\PowderPatterns\PhaseTransition\structure_000007.cif C:\Users\jacco\Documents\PowderPatterns\PhaseTransition\structure_000001_0.xye
        std::vector< std::string > extensions;
        extensions.push_back( "cif" );
        extensions.push_back( "cif" );
        extensions.push_back( "xye" );
        std::vector< FileName > input_file_names = sort_file_names_by_extension( argc, argv, extensions );
        CrystalStructure crystal_structure_1;
        read_cif( input_file_names[0], crystal_structure_1 );
        CrystalStructure crystal_structure_2;
        read_cif( input_file_names[1], crystal_structure_2 );
        PowderPattern powder_pattern;
        powder_pattern.read_xye( input_file_names[2] );
        TextFileWriter text_file_writer( replace_extension( input_file_names[2], "inp" ) );
        write_preamble( text_file_writer );
        size_t nbackground_terms = 20;
        for ( size_t i( 0 ); i != nbackground_terms; ++i )
            text_file_writer.write_line( "prm JvdS_BKGR_" + size_t2string( i+1, 2 ) + " 0.0" );
        text_file_writer.write_line( "prm JvdS_zero_point_error 0.02" );
        text_file_writer.write_line( "prm JvdS_filament_length 10.0" );
        text_file_writer.write_line( "prm JvdS_sample_length 2.0" );
        text_file_writer.write_line( "prm JvdS_receiving_slit_length 4.0" );
        text_file_writer.write_line( "prm JvdS_CS_G_1 500.0" );
        text_file_writer.write_line( "prm JvdS_CS_L_1 500.0" );
        text_file_writer.write_line( "prm JvdS_Strain_G_1 1.0" );
        text_file_writer.write_line( "prm JvdS_Strain_L_1 1.0" );
        write_unit_cell_variables( text_file_writer, crystal_structure_1.crystal_lattice() );
        text_file_writer.write_line( "prm bnonh_1 3.0" );
        text_file_writer.write_line( "prm bh_1 = 1.2 * bnonh_1;" );
        text_file_writer.write_line( "prm JvdS_scale_1 1" );
        text_file_writer.write_line( "prm JvdS_PO_1 1.0" );
        for ( size_t i( 0 ); i != crystal_structure_1.natoms(); ++i )
        {
            text_file_writer.write_line( "prm " + crystal_structure_1.atom( i ).label() + "_1_x " + double2string_pad_plus( crystal_structure_1.atom( i ).position().x(), 5, ' ' ) );
            text_file_writer.write_line( "prm " + crystal_structure_1.atom( i ).label() + "_1_y " + double2string_pad_plus( crystal_structure_1.atom( i ).position().y(), 5, ' ' ) );
            text_file_writer.write_line( "prm " + crystal_structure_1.atom( i ).label() + "_1_z " + double2string_pad_plus( crystal_structure_1.atom( i ).position().z(), 5, ' ' ) );
        }
        text_file_writer.write_line( "prm JvdS_CS_G_2 500.0" );
        text_file_writer.write_line( "prm JvdS_CS_L_2 500.0" );
        text_file_writer.write_line( "prm JvdS_Strain_G_2 1.0" );
        text_file_writer.write_line( "prm JvdS_Strain_L_2 1.0" );
        write_unit_cell_variables( text_file_writer, crystal_structure_2.crystal_lattice() );
        text_file_writer.write_line( "prm bnonh_2 3.0" );
        text_file_writer.write_line( "prm bh_2 = 1.2 * bnonh_2;" );
        text_file_writer.write_line( "prm JvdS_scale_2 1" );
        text_file_writer.write_line( "prm JvdS_PO_2 1.0" );
        for ( size_t i( 0 ); i != crystal_structure_2.natoms(); ++i )
        {
            text_file_writer.write_line( "prm " + crystal_structure_2.atom( i ).label() + "_2_x " + double2string_pad_plus( crystal_structure_2.atom( i ).position().x(), 5, ' ' ) );
            text_file_writer.write_line( "prm " + crystal_structure_2.atom( i ).label() + "_2_y " + double2string_pad_plus( crystal_structure_2.atom( i ).position().y(), 5, ' ' ) );
            text_file_writer.write_line( "prm " + crystal_structure_2.atom( i ).label() + "_2_z " + double2string_pad_plus( crystal_structure_2.atom( i ).position().z(), 5, ' ' ) );
        }
        text_file_writer.write_line( "prm inflection_point 0.0" );
        text_file_writer.write_line( "prm slope 1.0" );
        size_t n = 21;
        for ( size_t i( 0 ); i != n; ++i )
        {
            text_file_writer.write_line( "xdd " + input_file_names[2].name() + "_" + size_t2string( i+1, 2 ) + ".xye xye_format" );
            text_file_writer.write_line( "  bkg" );
            for ( size_t j( 0 ); j != nbackground_terms; ++j )
                text_file_writer.write_line( "  =JvdS_BKGR_" + size_t2string( j+1, 2 ) + ";" );
            text_file_writer.write_line( "  start_X       " + double2string( powder_pattern.two_theta_start().value_in_degrees() ) );
            text_file_writer.write_line( "  finish_X      " + double2string( powder_pattern.two_theta_end().value_in_degrees() ) );
            text_file_writer.write_line( "  x_calculation_step " + double2string( powder_pattern.average_two_theta_step().value_in_degrees() ) );
            text_file_writer.write_line( "  Zero_Error(@ , =JvdS_zero_point_error; )" );
            text_file_writer.write_line( "  LP_Factor( 26.5 )" );
            text_file_writer.write_line( "  axial_conv" );
            text_file_writer.write_line( "    filament_length =JvdS_filament_length;" );
            text_file_writer.write_line( "    sample_length =JvdS_sample_length;" );
            text_file_writer.write_line( "    receiving_slit_length =JvdS_receiving_slit_length;" );
            text_file_writer.write_line( "    axial_n_beta 50" );
            text_file_writer.write_line( "  lam" );
            text_file_writer.write_line( "    ymin_on_ymax 0.001" );
            text_file_writer.write_line( "    la 1 lo 1.540560" );
            text_file_writer.write_line( "  str" );
            text_file_writer.write_line( "    r_bragg 0.0" );
            text_file_writer.write_line( "    CS_G(@ , =JvdS_CS_G_1; )" );
            text_file_writer.write_line( "    CS_L(@ , =JvdS_CS_L_1; )" );
            text_file_writer.write_line( "    Strain_G(@ , =JvdS_Strain_G_1; )" );
            text_file_writer.write_line( "    Strain_L(@ , =JvdS_Strain_L_1; )" );
//    text_file_writer.write_line( "    prm  sh_scale_l" + aal + " 0.01" );
//    text_file_writer.write_line( "    spherical_harmonics_hkl sh_l" + aal );
//    text_file_writer.write_line( "      sh_order 6" );
//    text_file_writer.write_line( "    lor_fwhm = Abs( sh_scale_l" + aal + " * sh_l" + aal + " );" );
//    text_file_writer.write_line( "    prm  sh_scale_g" + aal + "  0.01" );
//    text_file_writer.write_line( "    spherical_harmonics_hkl sh_g" + aal );
//    text_file_writer.write_line( "      sh_order 6" );
//    text_file_writer.write_line( "    gauss_fwhm = Abs( sh_scale_g" + aal + " * sh_g" + aal + " );" );
            write_unit_cell( text_file_writer, crystal_structure_1.crystal_lattice(), true );
            text_file_writer.write_line( "    MVW( 0.0, 0.0, 0.0 )");
            text_file_writer.write_line( "    space_group \"" + remove( remove( crystal_structure_1.space_group().name(), '_' ), ' ') + "\"");
            text_file_writer.write_line( "    scale =(JvdS_scale_1 / ( 1.0 + Exp( -1.0 * slope * ( " + size_t2string( i ) + " - inflection_point ) ) ) )/100000;" );
//            text_file_writer.write_line( "'    PO_Spherical_Harmonics( sh, 6 )" );
            text_file_writer.write_line( "    PO(@ , =JvdS_PO_1;, , 0 0 1 )" );
            for ( size_t j( 0 ); j != crystal_structure_1.natoms(); ++j )
            {
                text_file_writer.write( "    site " + crystal_structure_1.atom( j ).label() + "_1 x =" + crystal_structure_1.atom( j ).label() + "_1_x" + ";" +
                                                                                                " y =" + crystal_structure_1.atom( j ).label() + "_1_y" + ";" +
                                                                                                " z =" + crystal_structure_1.atom( j ).label() + "_1_z" + ";" +
                                             " occ " + pad( crystal_structure_1.atom( j ).element().symbol(), 2, ' ' ) + " 1 beq = " );
                if ( crystal_structure_1.atom( j ).element().is_H_or_D() )
                    text_file_writer.write_line( "bh_1;" );
                else
                    text_file_writer.write_line( "bnonh_1;" );
            }
            if ( i == n-1 )
                text_file_writer.write_line( "    Out_CIF_STR( " + input_file_names[0].name() + "_RR.cif" + " )" );
            text_file_writer.write_line( "  str" );
            text_file_writer.write_line( "    r_bragg 0.0" );
            text_file_writer.write_line( "    CS_G(@ , =JvdS_CS_G_2; )" );
            text_file_writer.write_line( "    CS_L(@ , =JvdS_CS_L_2; )" );
            text_file_writer.write_line( "    Strain_G(@ , =JvdS_Strain_G_2; )" );
            text_file_writer.write_line( "    Strain_L(@ , =JvdS_Strain_L_2; )" );
//    text_file_writer.write_line( "    prm  sh_scale_l" + aal + " 0.01" );
//    text_file_writer.write_line( "    spherical_harmonics_hkl sh_l" + aal );
//    text_file_writer.write_line( "      sh_order 6" );
//    text_file_writer.write_line( "    lor_fwhm = Abs( sh_scale_l" + aal + " * sh_l" + aal + " );" );
//    text_file_writer.write_line( "    prm  sh_scale_g" + aal + "  0.01" );
//    text_file_writer.write_line( "    spherical_harmonics_hkl sh_g" + aal );
//    text_file_writer.write_line( "      sh_order 6" );
//    text_file_writer.write_line( "    gauss_fwhm = Abs( sh_scale_g" + aal + " * sh_g" + aal + " );" );
            write_unit_cell( text_file_writer, crystal_structure_2.crystal_lattice(), true );
            text_file_writer.write_line( "    MVW( 0.0, 0.0, 0.0 )");
            text_file_writer.write_line( "    space_group \"" + remove( remove( crystal_structure_2.space_group().name(), '_' ), ' ') + "\"");
            text_file_writer.write_line( "    scale = ( JvdS_scale_2 * ( 1.0 - ( 1.0 / ( 1.0 + Exp( -1.0 * slope * ( " + size_t2string( i ) + " - inflection_point ) ) ) ) ) )/100000;" );
//            text_file_writer.write_line( "'    PO_Spherical_Harmonics( sh, 6 )" );
            text_file_writer.write_line( "    PO(@ , =JvdS_PO_2;, , 0 0 1 )" );
            for ( size_t j( 0 ); j != crystal_structure_2.natoms(); ++j )
            {
                text_file_writer.write( "    site " + crystal_structure_2.atom( j ).label() + "_2 x =" + crystal_structure_2.atom( j ).label() + "_2_x" + ";" +
                                                                                                " y =" + crystal_structure_2.atom( j ).label() + "_2_y" + ";" +
                                                                                                " z =" + crystal_structure_2.atom( j ).label() + "_2_z" + ";" +
                                             " occ " + pad( crystal_structure_2.atom( j ).element().symbol(), 2, ' ' ) + " 1 beq = " );
                if ( crystal_structure_2.atom( j ).element().is_H_or_D() )
                    text_file_writer.write_line( "bh_2;" );
                else
                    text_file_writer.write_line( "bnonh_2;" );
            }
            if ( i == n-1 )
                text_file_writer.write_line( "    Out_CIF_STR( " + input_file_names[1].name() + "_RR.cif" + " )" );
        }
//        text_file_writer.write_line( "    prm !flatten_width 0" );
//        text_file_writer.write_line( "    prm !flatten_weight    100000" );
//        text_file_writer.write_line( "    'Flatten( C19 N11 N35 C47 H54 C50 H58 C45 H49 C34 C46, , 0.0, flatten_width, flatten_weight )" );
//        text_file_writer.write_line( "  xdd_out " + FileName( input_cif_file_name.directory(), input_cif_file_name.name() + "_profile", "txt" ).full_name() + " load out_record out_fmt out_eqn" );
//        text_file_writer.write_line( "  {" );
//        text_file_writer.write_line( "      \" %11.5f \" = X;" );
//        text_file_writer.write_line( "      \" %11.5f \" = Yobs;" );
//        text_file_writer.write_line( "      \" %11.5f \" = Ycalc;" );
//        text_file_writer.write_line( "      \" %11.5f\\n\" = SigmaYobs;" );
//        text_file_writer.write_line( "  }" );
//        text_file_writer.write_line( "  phase_out " + FileName( input_cif_file_name.directory(), input_cif_file_name.name() + "_tickmarks", "txt" ).full_name() + " load out_record out_fmt out_eqn" );
//        text_file_writer.write_line( "  {" );
//        text_file_writer.write_line( "      \" %11.5f -200\\n\" = 2.0 * Rad * Th;" );
//        text_file_writer.write_line( "  }" );
    MACRO_END_GAME

    try // Simulate a phase transition captured by XRPD.
    {
// C:\Users\jacco\Documents\PowderPatterns/PhaseTransition/structure_000001.cif C:\Users\jacco\Documents\PowderPatterns/PhaseTransition/structure_000007.cif
        std::vector< std::string > extensions;
        extensions.push_back( "cif" );
        extensions.push_back( "cif" );
        std::vector< FileName > input_file_names = sort_file_names_by_extension( argc, argv, extensions );
        CrystalStructure crystal_structure_1;
        read_cif( input_file_names[0], crystal_structure_1 );
        CrystalStructure crystal_structure_2;
        read_cif( input_file_names[1], crystal_structure_2 );
        crystal_structure_1.apply_space_group_symmetry();
        crystal_structure_2.apply_space_group_symmetry();

        bool save_all_for_debugging( true );

        RealisticXRPDSimulatorSettings settings_1;
        settings_1.set_wavelength( Wavelength() );
        settings_1.set_two_theta_start( Angle( 2.0, Angle::DEGREES ) );
        settings_1.set_two_theta_end( Angle( 40.0, Angle::DEGREES ) );
        settings_1.set_two_theta_step( Angle( 0.015, Angle::DEGREES ) );
        settings_1.set_FWHM( 0.2 );
        settings_1.set_zero_point_error( Angle( 0.02, Angle::DEGREES ) );
        bool include_PO = true;
        if ( include_PO )
            settings_1.set_preferred_orientation( select_realistic_preferred_orientation_direction( crystal_structure_1.crystal_lattice() ), 0.9 );
        if ( true )
            settings_1.set_finger_cox_jephcoat( 10.0 / 400.0, 10.0 / 400.0 ); // Finger-Cox-Jephcoat
        settings_1.set_include_background( true );
        settings_1.set_include_noise( true );
        settings_1.set_Bragg_total_signal_normalisation( 10000.0 );
        settings_1.set_background_total_signal_normalisation( 0.2*10000.0 );
        settings_1.set_highest_peak( 10000.0 );
//        if ( save_all_for_debugging )
//            realistic_XRPD_simulator.Bragg_diffraction().save_xye( replace_extension( append_to_file_name( input_file_names[0], "_Bragg" ), "xye" ), true );
//        if ( save_all_for_debugging )
//            realistic_XRPD_simulator.background().save_xye( replace_extension( append_to_file_name( input_file_names[0], "_background" ), "xye" ), true );
//        if ( save_all_for_debugging )
//            realistic_XRPD_simulator.noise().save_xye( replace_extension( append_to_file_name( input_file_names[0], "_noise" ), "xye" ), true );
        RealisticXRPDSimulatorSettings settings_2( settings_1 );
        if ( include_PO )
            settings_2.set_preferred_orientation( select_realistic_preferred_orientation_direction( crystal_structure_2.crystal_lattice() ), 0.9 );
        RealisticXRPDSimulator realistic_XRPD_simulator_1( crystal_structure_1, settings_1 );
        RealisticXRPDSimulator realistic_XRPD_simulator_2( crystal_structure_2, settings_2 );
        // We could let the amount of background increase to simulate that some of the crystalline phase
        // turns amorphous.
        NormalisedLogisticFunction sigmoidal_function( 0.0, 1.0 ); // Inflection point, slope.
        // Fairly good points to choose: -10, -9, -8, ..., 8, 9, 10
        for ( int t( 0 ); t != 1; ++t )
        {
            double scale_factor_structure_1 = sigmoidal_function( static_cast<double>( t ) );
            double scale_factor_structure_2 = 1.0 - scale_factor_structure_1;
            PowderPattern powder_pattern_1 = realistic_XRPD_simulator_1.calculate();
            powder_pattern_1.scale( scale_factor_structure_1 );
            PowderPattern powder_pattern_2 = realistic_XRPD_simulator_2.calculate();
            powder_pattern_2.scale( scale_factor_structure_2 );
            PowderPattern result = powder_pattern_1;
            result += powder_pattern_2;
            result.make_counts_integer();
            result.recalculate_estimated_standard_deviations();
            result.save_xye( replace_extension( append_to_file_name( input_file_names[0], "_" + int2string(t) ), "xye" ), true );
        }
    MACRO_END_GAME

    try // Split powder pattern over n powder patterns then only save one.
    {
        MACRO_ONE_XYEFILENAME_AS_ARGUMENT
        size_t n( 2 );
        std::vector< PowderPattern > powder_patterns = split( powder_pattern, n );
        powder_patterns[0].save_xye( append_to_file_name( input_file_name, "_1o2" ), true );
    MACRO_END_GAME

    try // Add class.
    {
        add_class( "RealisticXRPDSimulator" );
    MACRO_END_GAME

    try // Convert .cif to .xye.
    {
        if ( argc != 2 )
            throw std::runtime_error( "Please give the name of a .cif file." );
        FileName input_file_name( argv[ 1 ] );
        PowderPattern powder_pattern;
        powder_pattern.read_cif( input_file_name );
        powder_pattern.set_wavelength( Wavelength::determine_from_wavelength( 0.81906 ) );
        powder_pattern.save_xye( replace_extension( input_file_name, "xye" ), true );
    MACRO_END_GAME

    try // Insert one hydrogen atom between two atoms.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        if ( (true) )
        {
        std::string origin_atom_label( "O4" );
        std::string neighbour_atom_label( "O3" );
        Vector3D origin_atom_frac = crystal_structure.atom( crystal_structure.find_label( origin_atom_label ) ).position();
        Vector3D neighbour_atom_frac = crystal_structure.atom( crystal_structure.find_label( neighbour_atom_label ) ).position();
        // Find shortest distance between the two taking periodicity and space-group symmetry into account.
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
        // Find shortest distance between the two taking periodicity and space-group symmetry into account.
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
        // Find shortest distance between the two taking periodicity and space-group symmetry into account.
        double distance;
        Vector3D difference_frac;
        crystal_structure.second_shortest_distance( origin_atom_frac, neighbour_atom_frac, distance, difference_frac );
            std::cout << "Warning: distance = " << distance << std::endl;
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
        // Find shortest distance between the two taking periodicity and space-group symmetry into account.
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

    try // Calculate amount of PO.
    {
        double min_PO( 100000.0 );
        double max_PO( 0.0 );
        double r( 1.14 );
        size_t nsteps( 50 );
        for ( size_t i( 0 ); i != nsteps; ++i )
        {
            Angle alpha = ( static_cast<double>(i) / static_cast<double>(nsteps) ) * Angle::angle_180_degrees();
            double PO = std::pow( square(r) * square(alpha.cosine()) + square(alpha.sine())/r, -3.0/2.0 );
            if ( PO < min_PO )
                min_PO = PO;
            if ( PO > max_PO )
                max_PO = PO;
        }
        std::cout << min_PO/max_PO<< std::endl;
    MACRO_END_GAME

    try // Apply zero-point error correction to a powder pattern.
    {
        if ( argc != 3 )
        {
            std::cout << "This corrects a zero-point error by subtracting it, so the zero-point error you provide should probably be a positive value." << std::endl;
            throw std::runtime_error( "Please give the name of a powder diffraction pattern and a zero-point error." );
        }
        FileName input_file_name( argv[ 1 ] );
        PowderPattern powder_pattern;
        powder_pattern.read_xye( input_file_name );
        Angle zero_point_error = Angle::from_degrees( string2double( argv[ 2 ] ) );
        powder_pattern.correct_zero_point_error( zero_point_error );
        powder_pattern.save_xye( append_to_file_name( input_file_name, "_zp" ), true );
    MACRO_END_GAME

    try // Split powder pattern over n powder patterns.
    {
        MACRO_ONE_XYEFILENAME_AS_ARGUMENT
        size_t n( 3 );
        std::vector< PowderPattern > powder_patterns = split( powder_pattern, n );
        for ( size_t i( 0 ); i != n; ++i )
            powder_patterns[i].save_xye( append_to_file_name( input_file_name, "_" + size_t2string( i+1 ) ), true );
    MACRO_END_GAME

    try // Rebin powder pattern.
    {
        if ( argc != 3 )
            throw std::runtime_error( "Please give the name of a .xye file and the new bin size." );
        FileName input_file_name( argv[ 1 ] );
        PowderPattern powder_pattern;
        powder_pattern.read_xye( input_file_name );
        size_t bin_size = string2integer( argv[ 2 ] );
        powder_pattern.rebin( bin_size );
        powder_pattern.save_xye( append_to_file_name( input_file_name, "_rebin_" + size_t2string( bin_size ) ), true );
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
    crystal_structure.save_cif( FileName( file_list.value( 0 ).directory(), file_list.value( 0 ).name() + "_lean", "cif" ) );
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

    try // Calculate powder pattern from MD trajectory.
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

    try // Invert a matrix.
    {
        Matrix3D matrix( -1.0, -1.0,  0.0,
                          1.0, -1.0,  0.0,
                          1.0,  1.0,  1.0 );
        std::cout << "Determinant of the input matrix = " << matrix.determinant() << std::endl;
        matrix.invert();
        matrix.show();
        std::cout << "Determinant of the inverse matrix = " << matrix.determinant() << std::endl;
    MACRO_END_GAME

    try // Add two powder diffraction patterns.
    {
        std::vector< std::string > extensions;
        extensions.push_back( "xye" );
        extensions.push_back( "xye" );
        std::vector< FileName > input_file_names = sort_file_names_by_extension( argc, argv, extensions );
        std::vector< PowderPattern > powder_patterns;
        std::vector< double > noscp2ts;
        for ( size_t i( 0 ); i != input_file_names.size(); ++i )
        {
            PowderPattern powder_pattern;
            powder_pattern.read_xye( input_file_names[i] );
            powder_patterns.push_back( powder_pattern );
            noscp2ts.push_back( 1.0 );
        }
        PowderPattern result = add_powder_patterns( powder_patterns, noscp2ts );
        result.save_xye( FileName( "PowderPatternsAdded.xye" ), true );
    MACRO_END_GAME

    try // Sort FileList.txt.
    {
        if ( argc != 2 )
        {
            std::cout << "This sorts a FileList.txt if it contains names like:" << std::endl;
            std::cout << "structure_9.cif" << std::endl;
            std::cout << "structure_99.cif" << std::endl;
            std::cout << "because structure_100.cif would come before structure_99.cif." << std::endl;
            throw std::runtime_error( "Please give the name of a FileList.txt file." );
        }
        FileName file_list_file_name( argv[ 1 ] );
        FileList file_list( file_list_file_name );
        if ( file_list.empty() )
            throw std::runtime_error( std::string( "No files in file list " ) + file_list_file_name.full_name() );
        std::vector< size_t > sequential_numbers;
        sequential_numbers.reserve( file_list.size() );
        // Literally cut out the bit between "_" and ".".
        for ( size_t i( 0 ); i != file_list.size(); ++i )
        {
            std::string file_name = file_list.value( i ).file_name();
            size_t start_pos = file_name.find( '_' );
            if ( start_pos == std::string::npos )
                throw std::runtime_error( "Something was wrong with start_pos." );
            ++start_pos;
            size_t end_pos = file_name.rfind( '.' );
            if ( end_pos == std::string::npos )
                throw std::runtime_error( "Something was wrong with end_pos."  );
            std::string sequential_number_string = file_name.substr( start_pos, end_pos - start_pos );
            sequential_numbers.push_back( string2integer( sequential_number_string ) );
        }
        Mapping mapping = sort( sequential_numbers );
        FileList result;
        result.reserve( file_list.size() );
        for ( size_t i( 0 ); i != mapping.size(); ++i )
            result.push_back( file_list.value( mapping[i] ) );
        result.save( append_to_file_name( file_list_file_name, "_sorted" ) );
    MACRO_END_GAME

    try // Calculate ADPs from a set of frames (as cif files).
    {
        MACRO_ONE_FILELISTNAME_AS_ARGUMENT
        SpaceGroup space_group; // = SpaceGroup::P21c();
//        SpaceGroup space_group = SpaceGroup::P21c();
//        space_group.add_inversion_at_origin(); // Space group is now P-1
        // Stores results as .cif
        AnalyseTrajectory analyse_trajectory( file_list, 6, 12, 6, space_group );
  //      analyse_trajectory.save_centres_of_mass();
    MACRO_END_GAME

    try // Sudoku.
    {
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
        if ( false ) // AI Escargot 6313 guesses, stack pointer = 1
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

    try // Collapse supercell.
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
        read_cif( FileName( "Paracetamol_fr0001.cif" ), crystal_structure );
        std::cout << "Step 2" << std::endl;
        crystal_structure.collapse_supercell( crystal_lattice, space_group );
        std::cout << "Step 3" << std::endl;
        crystal_structure.save_xyz( FileName( "collapsed.xyz" ) );
    MACRO_END_GAME

    try // Primitive to centred.
    {
        std::vector< Centring > centrings;
        centrings.push_back( Centring( "A" ) );
        centrings.push_back( Centring( "B" ) );
        centrings.push_back( Centring( "C" ) );
        centrings.push_back( Centring( "I" ) );
        centrings.push_back( Centring( "F" ) );
        centrings.push_back( Centring( "R" ) );
        for ( size_t i( 0 ); i != centrings.size(); ++i )
        {
            Matrix3D centred2primitive = centrings[i].to_primitive();
            centred2primitive.invert();
            CrystalStructure crystal_structure;
            SpaceGroup space_group;
            crystal_structure.set_space_group( space_group );
            // Apply transformation to primitive.
            crystal_structure.transform( centred2primitive );
            space_group = crystal_structure.space_group();
            space_group.set_name( "" );
            add_centring_to_space_group_after_transformation( centred2primitive, space_group );
            crystal_structure.set_space_group( space_group );
            space_group.show();
        }
    MACRO_END_GAME

    try // Centred to primitive.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        if ( crystal_structure.space_group().centring().is_primitive() )
        {
            std::cout << "The structure is not centred, nothing to do." << std::endl;
            return 0;
        }
        Centring original_centring = crystal_structure.space_group().centring();
        crystal_structure.reduce_to_primitive();
        crystal_structure.save_cif( replace_extension( append_to_file_name( input_file_name, "_reduced_no90" ), "cif" ) );
        Matrix3D best_transformation_matrix = crystal_structure.crystal_lattice().choose_angles_close_to_90();
        crystal_structure.transform( best_transformation_matrix );
        // Write out crystal structure with new symmetry operators and in P1.
        crystal_structure.save_cif( replace_extension( append_to_file_name( input_file_name, "_reduced" ), "cif" ) );
        bool has_inversion_at_origin = crystal_structure.space_group().has_inversion_at_origin();
        crystal_structure.convert_to_P1();
        crystal_structure.save_cif( replace_extension( append_to_file_name( input_file_name, "_reduced_P1" ), "cif" ) );
        if ( has_inversion_at_origin )
        {
            SpaceGroup space_group;
            space_group.add_inversion_at_origin();
            space_group.set_name( "P-1" );
            crystal_structure.set_space_group( space_group );
            crystal_structure.reduce_to_asymmetric_unit( 0.01 );
            crystal_structure.save_cif( replace_extension( append_to_file_name( input_file_name, "_reduced_P-1" ), "cif" ) );
        }
        // Print inverse.
        std::cout << "Best transformation = " << std::endl;
        std::cout << best_transformation_matrix << std::endl;
        Matrix3D combined_transformation_matrix = best_transformation_matrix * original_centring.to_primitive();
        std::cout << "Combined transformation = " << std::endl;
        std::cout << combined_transformation_matrix << std::endl;
        std::cout << "Inverse = " << std::endl;
        std::cout << inverse( combined_transformation_matrix ) << std::endl;
    MACRO_END_GAME

    try // Centred to primitive.
    {
        MACRO_ONE_FILELISTNAME_OR_LIST_OF_FILES_AS_ARGUMENT
        for ( size_t i( 0 ); i != file_list.size(); ++i )
        {
            CrystalStructure crystal_structure;
            std::cout << "Now reading cif... " + file_list.value( i ).full_name() << std::endl;
            read_cif( file_list.value( i ), crystal_structure );
            if ( crystal_structure.space_group().centring().is_primitive() )
                continue;
            Centring original_centring = crystal_structure.space_group().centring();
            crystal_structure.reduce_to_primitive();
            crystal_structure.save_cif( replace_extension( append_to_file_name( file_list.value( i ), "_reduced_no90" ), "cif" ) );
            Matrix3D best_transformation_matrix = crystal_structure.crystal_lattice().choose_angles_close_to_90();
            crystal_structure.transform( best_transformation_matrix );
            // Write out crystal structure with new symmetry operators and in P1.
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
            Matrix3D combined_transformation_matrix = best_transformation_matrix * original_centring.to_primitive();
            std::cout << "Combined transformation = " << std::endl;
            std::cout << combined_transformation_matrix << std::endl;
            std::cout << "Inverse = " << std::endl;
            std::cout << inverse( combined_transformation_matrix ) << std::endl;
        }
    MACRO_END_GAME

    try // Calculate voids / volume as per Wood / Parsons paper.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        double void_volume_2 = void_volume( crystal_structure );
        std::cout << "Unit-cell volume = " << double2string( crystal_structure.crystal_lattice().volume() ) << std::endl;
        std::cout << "Void volume = " << double2string( void_volume_2 ) << std::endl;
        std::cout << "Network volume = " << double2string( crystal_structure.crystal_lattice().volume() - void_volume_2 ) << std::endl;
    MACRO_END_GAME

    try // Find voids.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        crystal_structure.apply_space_group_symmetry();
        double probe_radius = 1.75;
        double volume = find_voids( crystal_structure, probe_radius );
        std::cout << double2string( volume ) + " " + double2string( volume / crystal_structure.space_group().nsymmetry_operators() ) << std::endl;
    MACRO_END_GAME

    try // Find voids for a list of files.
    {
        MACRO_ONE_FILELISTNAME_OR_LIST_OF_FILES_AS_ARGUMENT
        std::cout << "WARNING: the molecular volume is estimated assuming that the smallest molecular volume corresponds to Z'=1." << std::endl;
        std::cout << "WARNING: if the smallest molecular volume corresponds to Z'>1 or Z'<1 then the results will be wrong." << std::endl;
        TextFileWriter text_file_writer( FileName( "Voids.txt" ) );
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
        Mapping sorted_map = sort( voids_volumes_per_Z );
        size_t iStart;
        for ( iStart = 0; iStart != nfiles; ++iStart )
        {
            if ( voids_volumes_per_Z[ sorted_map[ iStart ] ] > 20.0 )
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
                text_file_writer.write_line( identifiers[ sorted_map[ i ] ] + " " + double2string( voids_volumes_per_Z[ sorted_map[ i ] ] ) );
            text_file_writer.write_line();
            if ( (nfiles - iStart) == 1 )
            {
                text_file_writer.write( "Rank " );
                text_file_writer.write( size_t2string( sorted_map[ iStart ] + 1 ) );
                text_file_writer.write( " contains voids amounting to " );
                text_file_writer.write( double2string_2( voids_volumes_per_Z[ sorted_map[ iStart ] ], 0 ) );
                text_file_writer.write( " \\\\AA$^{3}$/Z." );
            }
            else
            {
//        Ranks 12, 22 5, 17, 1, 9 and 10 contain voids amounting to 20, 21, 21, 24, 28, 40 and 45 A3/Z, respectively.
                text_file_writer.write( "Ranks " );
                for ( size_t i( iStart ); i != nfiles; ++i )
                {
                    if ( i == nfiles - 1 )
                        text_file_writer.write( " and "  );
                    else if ( i != iStart )
                        text_file_writer.write( ", "  );
                    text_file_writer.write( size_t2string( sorted_map[ i ] + 1 ) );
                }
                text_file_writer.write( " contain voids amounting to " );
                for ( size_t i( iStart ); i != nfiles; ++i )
                {
                    if ( i == nfiles - 1 )
                        text_file_writer.write( " and "  );
                    else if ( i != iStart )
                        text_file_writer.write( ", "  );
                    text_file_writer.write( double2string_2( voids_volumes_per_Z[ sorted_map[ i ] ], 0 ) );
                }
                text_file_writer.write( " \\\\AA$^{3}$/Z, respectively." );
            }
            text_file_writer.write( " Of interest are voids that are greater than about 20 A3/Z: 21.5 A3/Z suffices to store a water molecule (at least in terms of volume), a chloride ion is about 25 A3/Z." );
            text_file_writer.write_line( " Voids between 15 and 20 A3/Z are quite common, but voids over 25 A3/Z are rare." );
        }
    MACRO_END_GAME

    try // Cyclic permute unit-cell axes.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        Matrix3D tranformation_matrix(  0.0,  1.0,  0.0,
                                        0.0,  0.0,  1.0,
                                        1.0,  0.0,  0.0 );
        crystal_structure.transform( tranformation_matrix );
        crystal_structure.save_cif( append_to_file_name( input_file_name, "_abc2bca" ) );
    MACRO_END_GAME

    try // Swap a- and c-axes.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        Matrix3D tranformation_matrix(  0.0,  0.0,  1.0,
                                        0.0, -1.0,  0.0,
                                        1.0,  0.0,  0.0 );
        crystal_structure.transform( tranformation_matrix );
        crystal_structure.save_cif( append_to_file_name( input_file_name, "_ac_swapped" ) );
    MACRO_END_GAME

    try // Generate and save all orthorhombic permutations.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        std::vector< Matrix3D > tranformation_matrices = orthorhombic_unit_cell_axes_permutations();
        std::vector< std::string > permutation_labels = orthorhombic_unit_cell_axes_permutation_labels();
        for ( size_t i( 1 ); i != tranformation_matrices.size(); ++i )
        {
            CrystalStructure crystal_structure_2( crystal_structure );
            crystal_structure_2.transform( tranformation_matrices[i] );
            crystal_structure_2.save_cif( append_to_file_name( input_file_name, "_" + permutation_labels[i] ) );
        }
    MACRO_END_GAME

    try // Copy hydrogen atoms from minimisation to exp. structure.
    {
        if ( argc != 3 )
            throw std::runtime_error( "Please give the name of two .cif files, _exp.cif and _mi.cif." );
// We first have to find out if there was any drift.
        FileName input_file_name_exp( argv[ 1 ] );
        CrystalStructure crystal_structure_exp;
        read_cif_or_cell( input_file_name_exp, crystal_structure_exp );
        FileName input_file_name_mi( argv[ 2 ] );
        CrystalStructure crystal_structure_mi;
        read_cif_or_cell( input_file_name_mi, crystal_structure_mi );
        std::vector< Vector3D > drift_vectors;
        for ( size_t i( 0 ); i != crystal_structure_mi.natoms(); ++i )
        {
            if ( crystal_structure_mi.atom( i ).element().is_H_or_D() )
                continue;
            Vector3D exp_frac = crystal_structure_exp.atom( i ).position();
            Vector3D mi_frac = crystal_structure_mi.atom( i ).position();
            drift_vectors.push_back( mi_frac - exp_frac );
        }
        Vector3D average_drift_vector = average( drift_vectors );
        std::cout << "Average drift vector = " << average_drift_vector << std::endl;
        for ( size_t i( 0 ); i != drift_vectors.size(); ++i)
        {
            if ( ! nearly_equal( drift_vectors[i], average_drift_vector, 0.0001 ) )
            {
                std::cout << "Drift vector differs by more than 0.0001" << std::endl;
            }
        }
        TextFile exp_file( input_file_name_exp );
        std::string line;
        std::vector< std::string > words;
        for ( size_t i( 0 ); i != crystal_structure_mi.natoms(); ++i )
        {
            Atom atom = crystal_structure_mi.atom( i );
            if ( ! crystal_structure_mi.atom( i ).element().is_H_or_D() )
                continue;
            size_t iLine = exp_file.find( atom.label() + " H " );
            if ( iLine == std::string::npos )
                throw std::runtime_error( "String not found: |" + atom.label() + " H |");
            line = exp_file.line( iLine );
            words = split( line );
            if ( words.size() != 8 )
                continue;
            line = pad( words[0], 4 ) + " " + pad( words[1], 2 ) + " " + words[2] + " " +
                   pad( double2string_pad_plus( atom.position().x() - average_drift_vector.x(), 5, ' ' ), 11 ) + " " +
                   pad( double2string_pad_plus( atom.position().y() - average_drift_vector.y(), 5, ' ' ), 11 ) + " " +
                   pad( double2string_pad_plus( atom.position().z() - average_drift_vector.z(), 5, ' ' ), 12 ) + " " + words[6] + " " + words[7];
            exp_file.set_line( iLine, line );
        }
        exp_file.save( append_to_file_name( input_file_name_exp, "_Hcal" ) );
    MACRO_END_GAME

    try // Fix non-hydrogen atoms. Initialise H atoms positions from energy-minimised structure.
    {
        if ( argc != 3 )
            throw std::runtime_error( "Please give the name of two .cif files, _exp.cif and _mi.cif." );
        FileName input_file_name_exp( argv[ 1 ] );
        CrystalStructure crystal_structure_exp;
        read_cif_or_cell( input_file_name_exp, crystal_structure_exp );
        FileName input_file_name_mi( argv[ 2 ] );
        CrystalStructure crystal_structure_mi;
        read_cif_or_cell( input_file_name_mi, crystal_structure_mi );
        for ( size_t i( 0 ); i != crystal_structure_mi.natoms(); ++i )
        {
            Atom new_atom = crystal_structure_exp.atom( i );
            if ( crystal_structure_mi.atom( i ).element().is_H_or_D() )
            {
                Vector3D H_mi_cart = crystal_structure_mi.crystal_lattice().fractional_to_orthogonal( crystal_structure_mi.atom( i ).position() );
                size_t j = crystal_structure_mi.nearest_atom( i );
                Vector3D X_mi_cart = crystal_structure_mi.crystal_lattice().fractional_to_orthogonal( crystal_structure_mi.atom( j ).position() );
                Vector3D X_exp_cart = crystal_structure_exp.crystal_lattice().fractional_to_orthogonal( crystal_structure_exp.atom( j ).position() );
                new_atom.set_position( crystal_structure_exp.crystal_lattice().orthogonal_to_fractional( X_exp_cart - X_mi_cart + H_mi_cart ) );
            }
            else
                new_atom.set_label( new_atom.label() + "_fx" );
            crystal_structure_exp.set_atom( i, new_atom );
        }
        crystal_structure_exp.save_cif( replace_extension( append_to_file_name( input_file_name_exp, "_nHfix" ) , "cif" ) );
    MACRO_END_GAME

    try // Fix non-hydrogen atoms.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
        {
            if ( crystal_structure.atom( i ).element().is_H_or_D() )
                continue;
            Atom new_atom( crystal_structure.atom( i ) );
            new_atom.set_label( new_atom.label() + "_fx" );
            crystal_structure.set_atom( i, new_atom );
        }
        crystal_structure.save_cif( replace_extension( append_to_file_name( input_file_name, "_nHfix" ) , "cif" ) );
    MACRO_END_GAME

    try // Determine shift in GRACE when non-H atoms fixed.
    {
        if ( argc != 3 )
            throw std::runtime_error( "Please give the name of two .cif files, as input and from GRACE." );
        FileName input_file_name_input( argv[ 1 ] );
        CrystalStructure crystal_structure_input;
        read_cif_or_cell( input_file_name_input, crystal_structure_input );
        FileName input_file_name_GRACE( argv[ 2 ] );
        CrystalStructure crystal_structure_GRACE;
        read_cif_or_cell( input_file_name_GRACE, crystal_structure_GRACE );
        for ( size_t i( 0 ); i != crystal_structure_GRACE.natoms(); ++i )
        {
            if ( crystal_structure_input.atom( i ).element().is_H_or_D() )
                continue;
            Vector3D input_cart = crystal_structure_input.crystal_lattice().fractional_to_orthogonal( crystal_structure_input.atom( i ).position() );
            Vector3D GRACE_cart = crystal_structure_GRACE.crystal_lattice().fractional_to_orthogonal( crystal_structure_GRACE.atom( i ).position() );
            std::cout << GRACE_cart - input_cart << std::endl;
        }
    MACRO_END_GAME

    try // Analyse a whole T series at once: calculate average structures and ADPs from a set of frames (as cif files).
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

    try // Convert Biso to Uiso, including ESDs.
    {
        TextFileReader_2 input_file( FileName( "C:\\Users\\jacco\\Documents\\Research\\Lipitor\\Lipitor_Final_Biso.cif" ) );
        double conversion_factor = 1.0 / ( 8.0 * CONSTANT_PI * CONSTANT_PI );
        std::vector< std::string > words;
        for ( size_t i( 0 ); i != input_file.size(); ++i )
        {
            words = split( input_file.line( i ) );
            if ( words.size() != 1 )
                continue;
            DoubleWithESD dwe_1( words[0] );
            DoubleWithESD dwe_2( conversion_factor * dwe_1.value(), conversion_factor * dwe_1.estimated_standard_deviation() );
            std::cout << dwe_2.crystallographic_style() << std::endl;
        }
    MACRO_END_GAME

    try // Fit exponential.
    {
        TextFileReader text_file_reader( FileName( "./Chlorpropamide_Generation_table.txt" ) );
        std::vector< double > all_energies;
        std::vector< double > Z1_energies;
        std::vector< double > Z2_energies;
        std::string line;
        std::vector< std::string > words;
        double lowest_energy;
        size_t natoms = 30;
        // First three lines are headers.
        if ( ! text_file_reader.get_next_line( line ) )
            throw std::runtime_error( "File is empty." );
        if ( ! text_file_reader.get_next_line( line ) )
            throw std::runtime_error( "File is empty." );
        if ( ! text_file_reader.get_next_line( line ) )
            throw std::runtime_error( "File is empty." );
        if ( text_file_reader.get_next_line( line ) )
        {
            words = split( line, '|' );
            strip( words );
            lowest_energy = string2double( words[1] ) * natoms;
            size_t Z_prime = string2integer( words[4] ) / natoms;
            double energy = 0.0;
            all_energies.push_back( energy );
            if ( Z_prime == 1 )
                Z1_energies.push_back( energy );
            else if ( Z_prime == 2 )
                Z2_energies.push_back( energy );
        }
        else
            throw std::runtime_error( "File is empty." );
        while ( text_file_reader.get_next_line( line ) )
        {
            words = split( line, '|' );
            strip( words );
            size_t Z_prime = string2integer( words[4] ) / natoms;
            double energy = string2double( words[1] ) * natoms - lowest_energy;
            all_energies.push_back( energy );
            if ( Z_prime == 1 )
                Z1_energies.push_back( energy );
            else if ( Z_prime == 2 )
                Z2_energies.push_back( energy );
        }
        { // Scoping brackets
        size_t nbins = 26;
        Histogram histogram( 0.0, all_energies[ all_energies.size()-1 ], nbins );
        histogram.add_data( all_energies );
        histogram.plot();
        histogram.show();
        std::vector< double > x_values;
        std::vector< double > y_values;
        // The last bin is not completely filled, so we have to discard it, and because
        // everything is zero-based, the last bin is nbins-1, so we need nbins-2.
        histogram.values( 0, nbins-2, x_values, y_values );
        // x values, y values, returns a*exp(b*x)
        ExponentialFunction exponential_function = fit_exponential( x_values, y_values );
        double delta = histogram.middle_of_bin( 1 ) - histogram.middle_of_bin( 0 );
        std::cout << "a = " << exponential_function.a() << std::endl;
        std::cout << "b = " << exponential_function.b() << std::endl;
        }
        { // Scoping brackets
        size_t nbins = 26;
        Histogram histogram( 0.0, all_energies[ all_energies.size()-1 ], nbins );
        histogram.add_data( Z1_energies );
        histogram.plot();
        histogram.show();
        std::vector< double > x_values;
        std::vector< double > y_values;
        // The last bin is not completely filled, so we have to discard it, and because
        // everything is zero-based, the last bin is nbins-1, so we need nbins-2.
        histogram.values( 0, nbins-2, x_values, y_values );
        // x values, y values, returns a*exp(b*x)
        ExponentialFunction exponential_function = fit_exponential( x_values, y_values );
        double delta = histogram.middle_of_bin( 1 ) - histogram.middle_of_bin( 0 );
        std::cout << "a = " << exponential_function.a() << std::endl;
        std::cout << "b = " << exponential_function.b() << std::endl;
        }
        { // Scoping brackets
        size_t nbins = 26;
        Histogram histogram( 0.0, all_energies[ all_energies.size()-1 ], nbins );
        histogram.add_data( Z2_energies );
        histogram.plot();
        histogram.show();
        std::vector< double > x_values;
        std::vector< double > y_values;
        // The last bin is not completely filled, so we have to discard it, and because
        // everything is zero-based, the last bin is nbins-1, so we need nbins-2.
        histogram.values( 0, nbins-2, x_values, y_values );
        // x values, y values, returns a*exp(b*x)
        ExponentialFunction exponential_function = fit_exponential( x_values, y_values );
        double delta = histogram.middle_of_bin( 1 ) - histogram.middle_of_bin( 0 );
        std::cout << "a = " << exponential_function.a() << std::endl;
        std::cout << "b = " << exponential_function.b() << std::endl;
        }
    MACRO_END_GAME

    try // Find unit-cell transformation.
    {
        if ( argc != 3 )
            throw std::runtime_error( "Please give the name of two .cif files, the second is the target." );
        FileName input_file_name( argv[ 1 ] );
        CrystalStructure crystal_structure;
        read_cif_or_cell( input_file_name, crystal_structure );
        CrystalLattice old_crystal_lattice = crystal_structure.crystal_lattice();
        FileName file_name_2( argv[ 2 ] );
        CrystalStructure crystal_structure_2;
        read_cif_or_cell( file_name_2, crystal_structure_2 );
        CrystalLattice target_crystal_lattice = crystal_structure_2.crystal_lattice();
        double determinant = target_crystal_lattice.volume() / old_crystal_lattice.volume();
        Fraction fraction = Farey( determinant, 8 );
        determinant = fraction.to_double();
        std::cout << "Determinant = " << determinant << std::endl;
        if ( determinant + 0.000001 < 1.0 )
            throw std::runtime_error( "The determinant is < 1, that is not possible with the current algorithm." );
        double length_tolerance_percent( 10.0 );
        Angle angle_tolerance = Angle::from_degrees( 10.0 );
        double best_FoM( 1000000.0 );
        Matrix3D best_transformation_matrix;
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
                    Matrix3D transformation_matrix( i1, i2, i3, j1, j2, j3, k1, k2, k3 );
                    if ( ! nearly_equal( transformation_matrix.determinant(), determinant ) )
                        continue;
                    // Make a copy.
                    CrystalLattice new_lattice( old_crystal_lattice );
                    new_lattice.transform( transformation_matrix );
                    if ( ! nearly_equal( new_lattice, target_crystal_lattice, length_tolerance_percent, angle_tolerance ) )
                        continue;
                    {
                        double FoM = ((target_crystal_lattice.a_vector()+target_crystal_lattice.b_vector()+target_crystal_lattice.c_vector()) - (new_lattice.a_vector()+new_lattice.b_vector()+new_lattice.c_vector())).length();
                        transformation_matrix.show();
                        std::cout << "Inverse =" << std::endl;
                        inverse( transformation_matrix ).show();
                        new_lattice.print();
                        std::cout << "FoM = " << FoM << std::endl;
                        std::cout << std::endl;
                        if ( FoM < best_FoM )
                        {
                            best_FoM = FoM;
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
        }
        crystal_structure.transform( best_transformation_matrix );
        SpaceGroup space_group = crystal_structure.space_group();
        space_group.set_name( "" );
        crystal_structure.set_space_group( space_group );
        crystal_structure.save_cif( replace_extension( append_to_file_name( input_file_name, "_transformed" ) , "cif" ) );
    MACRO_END_GAME

    try // Find transformation using similarity as a target.
    {
        if ( argc < 2 )
            throw std::runtime_error( "Please give the name of two .cif files, the first is the target." );
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
                // Make a copy.
                CrystalStructure new_crystal_structure( crystal_structure );
                // Transform.
                for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
                {
                    Atom new_atom( crystal_structure.atom( i ) );
                    Vector3D current_position = crystal_structure.space_group().symmetry_operator( k ) * ( crystal_structure.atom( i ).position() + shifts[iShifts] );
                    new_atom.set_position( current_position );
                    new_crystal_structure.set_atom( i, new_atom );
                }
                new_crystal_structure.apply_space_group_symmetry();
                PowderPatternCalculator powder_pattern_calculator( new_crystal_structure );
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

    try // Find unit-cell transformation with a space-group setting as the target.
    {
        if ( argc < 2 )
            throw std::runtime_error( "Space-group setting search. Please give the name of two .cif files, the second is the target." );
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
                    Matrix3D transformation_matrix( i1, i2, i3, j1, j2, j3, k1, k2, k3 );
                    if ( ! nearly_equal( transformation_matrix.determinant(), 1.0 ) )
                        continue;
                    // Make a copy.
                    SpaceGroup new_space_group( old_space_group );
                    // Transform.
                    Matrix3D transformation_matrix_inverse_transpose( transformation_matrix );
                    transformation_matrix_inverse_transpose.invert();
                    transformation_matrix_inverse_transpose.transpose();

                    new_space_group.apply_similarity_transformation( transformation_matrix_inverse_transpose );

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

    try // Find unit-cell with all angles greater or smaller than 90 degrees (necessary to e.g. convert P1 to standard setting).
    // Also allows unit-cell axes to be ordered by length.
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
                    Matrix3D transformation_matrix( i1, i2, i3, j1, j2, j3, k1, k2, k3 );
                    if ( ! nearly_equal( transformation_matrix.determinant(), 1.0 ) )
                        continue;
                    // Make a copy.
                    CrystalLattice new_crystal_lattice( old_crystal_lattice );
                    // Transform.
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

    try // Calculate normalised cross correlation for two powder patterns.
    {
        if ( argc < 3 )
            throw std::runtime_error( "Please give the names of two .xye files triangle_width." );
        FileName file_name_1( argv[ 1 ] );
        FileName file_name_2( argv[ 2 ] );
        PowderPattern powder_pattern_1;
        PowderPattern powder_pattern_2;
        powder_pattern_1.read_xye( file_name_1 );
        powder_pattern_2.read_xye( file_name_2 );
        Angle triangle_width = Angle( 3.0, Angle::DEGREES );
        if ( argc > 3 )
            triangle_width = Angle( string2double( argv[ 3 ] ), Angle::DEGREES );
        std::cout << "NWCC = " << normalised_weighted_cross_correlation( powder_pattern_1, powder_pattern_2, triangle_width ) << std::endl;
        std::cout << "Rwp =  " << Rwp( powder_pattern_1, powder_pattern_2 ) << std::endl;
    MACRO_END_GAME

    try // Change one bond length.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        Atom origin = crystal_structure.atom( crystal_structure.find_label( "C43A" ) );
        Atom to_be_moved = crystal_structure.atom( crystal_structure.find_label( "C44" ) );
        double smallest_distance2 = crystal_structure.crystal_lattice().shortest_distance2( to_be_moved.position(), origin.position() );
        double smallest_distance = sqrt( smallest_distance2 );
        if ( smallest_distance < 0.001 )
            throw std::runtime_error( "Points too close together." );
        if ( smallest_distance > 3.0 )
            throw std::runtime_error( "Atoms not bound." );
        // We need to be able to set the length of the bond, so we must work in Cartesian coordinates.
        Vector3D difference_frac;
        double distance;
        crystal_structure.shortest_distance( origin.position(), to_be_moved.position(), distance, difference_frac );
        if ( ! nearly_equal( distance, smallest_distance ) )
            throw std::runtime_error( "Distances differ." );
        Vector3D difference_cart = crystal_structure.crystal_lattice().fractional_to_orthogonal( difference_frac );
        double target_bond_length( 1.47 );
        difference_cart *= target_bond_length / distance;
        difference_frac = crystal_structure.crystal_lattice().orthogonal_to_fractional( difference_cart );
        Vector3D new_atom_frac = origin.position() + difference_frac;
        Atom new_atom = to_be_moved;
        new_atom.set_position( new_atom_frac );
        crystal_structure.set_atom( crystal_structure.find_label( "C44" ), new_atom );
        crystal_structure.save_cif( append_to_file_name( input_file_name, "_BondChanged" ) );
    MACRO_END_GAME

    try // RMSCD with / without matching.
    {
        if ( argc < 3 )
            throw std::runtime_error( "Please give the names of two .cif files shift_steps add_inversion." );
        FileName file_name_1( argv[ 1 ] );
        CrystalStructure crystal_structure_1;
        read_cif_or_cell( file_name_1, crystal_structure_1 );
        if ( to_upper( file_name_1.extension() ) == "CELL" )
            crystal_structure_1.save_cif( replace_extension( file_name_1, "cif" ) );
        else
            crystal_structure_1.reduce_to_asymmetric_unit();
        FileName file_name_2( argv[ 2 ] );
        CrystalStructure crystal_structure_2;
        read_cif_or_cell( file_name_2, crystal_structure_2 );
        if ( to_upper( file_name_2.extension() ) == "CELL" )
            crystal_structure_2.save_cif( replace_extension( file_name_2, "cif" ) );
        else
            crystal_structure_1.reduce_to_asymmetric_unit();
        size_t shift_steps( 2 );
        if ( argc > 3 )
            shift_steps = string2integer( argv[ 3 ] );
        bool add_inversion( false );
        if ( argc > 4 )
            add_inversion = string2bool( argv[ 4 ] );
        double result;
        try
        {
        result = root_mean_square_Cartesian_displacement( crystal_structure_1, crystal_structure_2, false );
        std::cout << "Without matching RMSCD = " << result << " A" << std::endl;
        }
        catch ( std::exception & e )
        {
        }
        try
        {
        result = RMSCD_with_matching( crystal_structure_1, crystal_structure_2, shift_steps, add_inversion, false );
        std::cout << "With matching RMSCD = " << result << " A" << std::endl;
        }
        catch ( std::exception & e )
        {
        }
    MACRO_END_GAME

    try // Interpolate value.
    {
        if ( argc != 6 )
            throw std::runtime_error( "Usage: x1 y1 x2 y2 x." );
        double x1 = string2double( argv[ 1 ] );
        double y1 = string2double( argv[ 2 ] );
        double x2 = string2double( argv[ 3 ] );
        double y2 = string2double( argv[ 4 ] );
        double x  = string2double( argv[ 5 ] );
        std::vector< double > x_s;
        x_s.push_back( x1 );
        x_s.push_back( x2 );
        std::vector< double > y_s;
        y_s.push_back( y1 );
        y_s.push_back( y2 );
        LinearFunction linear_function = linear_regression( x_s, y_s );
        std::cout << double2string( linear_function( x ), 12 ) << std::endl;
    MACRO_END_GAME

    // Detect pseudo-inversion symmetry.
    // Crude algorithm, only works if there are only two fragments in the asymmetric unit.
    // Assumes that the first half of the atoms is one molecule in the asymmetric unit,
    // the second half of the atoms the second molecule.
    // The algorithm to detect floating axes works for up to and including orthorhombic, I do not know about the other space groups.
    try
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        if ( is_odd( crystal_structure.natoms() ) )
            throw std::runtime_error( "Z' must be 2 and number of atoms must therefore be even." );
        CrystalStructure original_crystal_structure( crystal_structure );
        SpaceGroup original_space_group = crystal_structure.space_group();
        bool at_least_one_floating_axis( false );
        for ( size_t i( 0 ); i != 3; ++i )
        {
            if ( original_space_group.is_floating_axis( i ) )
            {
                std::cout << "Floating axis found " << Vector3D::index2string( i ) << std::endl;
                at_least_one_floating_axis = true;
            }
        }
        for ( size_t iSymmOp( 0 ); iSymmOp != original_space_group.nsymmetry_operators(); ++iSymmOp )
        {
            crystal_structure = original_crystal_structure;
            // Split the atoms into two. Assume first half is one molecule, second half is the other.
            // Apply each symmetry operator (including the identity) in turn to the second molecule.
            for ( size_t i( crystal_structure.natoms() / 2 ); i != crystal_structure.natoms(); ++i )
            {
                Atom new_atom( crystal_structure.atom( i ) );
                new_atom.set_position( original_space_group.symmetry_operator( iSymmOp ) * crystal_structure.atom( i ).position() );
                if ( new_atom.ADPs_type() == Atom::ANISOTROPIC )
                    new_atom.set_anisotropic_displacement_parameters( rotate_adps( new_atom.anisotropic_displacement_parameters(), original_space_group.symmetry_operator( iSymmOp ).rotation(), crystal_structure.crystal_lattice() ) );
                crystal_structure.set_atom( i, new_atom );
            }
            SpaceGroup space_group = crystal_structure.space_group();
            Vector3D com = crystal_structure.centre_of_mass( true );
            std::cout << "Centre of mass = " << std::endl;
            com.show();
            Vector3D shift; // Floating axes are set to -(c.o.m.).
            Vector3D translation_for_symmetry_operators; // Floating axes are set to 0.0.
            for ( size_t i( 0 ); i != 3; ++i )
            {
                if ( original_space_group.is_floating_axis( i ) )
                {
                    shift.set_value( i, -com.value(i) );
                }
                else
                {
                    // Not a floating axis
                    // Round to nearest 1/12
                    Fraction granularity( 1, 12 );
                    if ( ( original_crystal_structure.crystal_lattice().lattice_system() == CrystalLattice::TRICLINIC ) ||
                         ( original_crystal_structure.crystal_lattice().lattice_system() == CrystalLattice::MONOCLINIC_A ) ||
                         ( original_crystal_structure.crystal_lattice().lattice_system() == CrystalLattice::MONOCLINIC_B ) ||
                         ( original_crystal_structure.crystal_lattice().lattice_system() == CrystalLattice::MONOCLINIC_C ) ||
                         ( original_crystal_structure.crystal_lattice().lattice_system() == CrystalLattice::ORTHORHOMBIC ) )
                        granularity = Fraction( 1, 12 );
                    Fraction fraction = double2fraction( -com.value(i), granularity );
                    std::cout << "Shift rounded to a fraction = " + fraction.to_string() << std::endl;
                    shift.set_value( i, fraction.to_double() );
                    translation_for_symmetry_operators.set_value( i, shift.value( i ) );
                }
            }
            SymmetryOperator symmetry_operator( Matrix3D(), translation_for_symmetry_operators ); // "x-1/4,y,z-3/4".
            // In Mercury, if the space-group name and the set of symmetry operators do not match up,
            // the space-group name takes precedence, so we have to erase it to ensure that the
            // symmetry operators are used instead.
            space_group.set_name( "" );
            space_group.apply_similarity_transformation( symmetry_operator );
            space_group.add_inversion_at_origin();
            crystal_structure.set_space_group( space_group );
            for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
            {
                Atom new_atom( crystal_structure.atom( i ) );
                new_atom.set_position( ( crystal_structure.atom( i ).position() ) + shift );
                crystal_structure.set_atom( i, new_atom );
            }
            // Up to this point, the shift along any floating axis was based on the c.o.m. of the entire molecule.
            // Some atoms will now overlap almost perfectly, those close to chiral centres will not overlap at all.
            // The shift that we applied was distorted by these atoms near the chiral centres that do not overlap at all.
            // Here we try to calculate better shifts, based only on the atoms that overlap well.
            // First, this is only relevant for floating axes. Second, this will only work for one of the symmetry operators (c.f. the loop over iSymmOp).
            // "Almost perfect overlap" should mean a small distance between atoms of the same element and only a single unambiguous match.
            // After we have applied this better shift, some more atoms, that previously did not overlap within the tolerance because the inclusion of the atoms near
            // the chiral centres had distorted our shift, will now suddenly also match within the tolerance, so this has to be done iteratively.
            // (Well, strictly speaking we simply need an algorithm that mathes molecular topologies, but currently we do not even have topologies so that is not going to work).
            // First iteration. For each atom we find the atom that is closest. Atoms with multiple matches or with no match are discarded. If the elements do not
            // match up, then that is also fishy, so we discard those as well. From the atom pairs that are left, sort them by distance, then take the top 25% of all atoms
            // (not just of the atoms that were left) and use those to calculate better shifts.
            // The "25%" is configurable, and the tolerance for matching is also a parameter.

            // Skip the whole exercise if there were no floating axes to start with.
            if ( at_least_one_floating_axis )
            {
                double tolerance = 0.3;
                double fraction = 0.2;
                std::vector< bool > done( crystal_structure.natoms(), false ); // We only use the second half, but this is way easier to program.
                std::vector< bool > is_dodgy( crystal_structure.natoms() / 2, false );
                std::vector< size_t > best_matches( crystal_structure.natoms() / 2, 0 );
                std::vector< double > best_distances( crystal_structure.natoms() / 2, 0.0 );
                for ( size_t iAtom1( 0 ); iAtom1 != crystal_structure.natoms() / 2; ++iAtom1 )
                {
                    size_t best_match = crystal_structure.natoms() / 2;
                    double best_distance = crystal_structure.crystal_lattice().shortest_distance( crystal_structure.atom( iAtom1 ).position(), crystal_structure.atom( best_match ).position() );
                    for ( size_t iAtom2( ( crystal_structure.natoms() / 2 ) + 1 ); iAtom2 != crystal_structure.natoms(); ++iAtom2 )
                    {
                        double distance = crystal_structure.crystal_lattice().shortest_distance( crystal_structure.atom( iAtom1 ).position(), -crystal_structure.atom( iAtom2 ).position() );
                        if ( distance < best_distance )
                        {
                            best_distance = distance;
                            best_match = iAtom2;
                        }
                    }
                    best_matches[ iAtom1 ] = best_match;
                    best_distances[ iAtom1 ] = best_distance;
                    if ( best_distance > tolerance )
                        is_dodgy[ iAtom1 ] = true;
                    if ( done[ best_match ] )
                    {
                        // This iAtom2 has more than one match-all pairs involving it are dodgy.
                        for ( size_t i( 0 ); i != iAtom1 + 1; ++i )
                        {
                            if ( best_matches[ i ] == best_match )
                                is_dodgy[ i ] = true;
                        }
                    }
                    else
                    {
                        // Unique match. Check the elements.
                        if ( crystal_structure.atom( iAtom1 ).element() != crystal_structure.atom( best_match ).element() )
                            is_dodgy[ iAtom1 ] = true;
                        done[ best_match ] = true;
                    }
                }
                // Remove dodgy matches.
                std::vector< size_t > atoms_1;
                std::vector< size_t > atoms_2;
                std::vector< double > distances;
                for ( size_t i( 0 ); i != crystal_structure.natoms() / 2; ++i )
                {
                    if ( ! is_dodgy[ i ] )
                    {
                        atoms_1.push_back( i );
                        atoms_2.push_back( best_matches[ i ] );
                        distances.push_back( best_distances[ i ] );
                    }
                }
                // Has at least 25% survived?
                if ( ( static_cast<double>( atoms_1.size() ) / ( crystal_structure.natoms() / 2.0 ) ) > fraction )
                {
                    // Calculate the c.o.m.. based only on the atoms that match very well.
                    // We have two choices here: either we only use the top 25% of the matches
                    // (which must all match within tolerance), assuming that they
                    // yield the best value possible, or we use all matches up to a certain distance.
                    // We just use all of them and we do not do any iterations.
                    Vector3D com_2;
                    for ( size_t i( 0 ); i != atoms_1.size(); ++i )
                    {
                        com_2 += crystal_structure.atom( atoms_1[ i ] ).position();
                        com_2 += crystal_structure.atom( atoms_2[ i ] ).position();
                    }
                    com_2 /= ( atoms_1.size() * 2.0 );
                    std::cout << "Refined centre of mass = " << std::endl;
                    com_2.show();
                    Vector3D shift_2; // Floating axes are set to -(c.o.m.).
                    for ( size_t i( 0 ); i != 3; ++i )
                    {
                        if ( original_space_group.is_floating_axis( i ) )
                            shift_2.set_value( i, -com_2.value(i) );
                    }
                    for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
                    {
                        Atom new_atom( crystal_structure.atom( i ) );
                        new_atom.set_position( ( crystal_structure.atom( i ).position() ) + shift_2 );
                        crystal_structure.set_atom( i, new_atom );
                    }
                }
            }
            crystal_structure.save_cif( append_to_file_name( input_file_name, "_" + size_t2string( iSymmOp ) + "_inverse" ) );
            for ( size_t i( crystal_structure.natoms() / 2 ); i != crystal_structure.natoms(); ++i )
                crystal_structure.set_suppressed( i, true );
            crystal_structure.save_cif( append_to_file_name( input_file_name, "_" + size_t2string( iSymmOp ) + "_inverse_1" ) );
            for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
                crystal_structure.set_suppressed( i, !crystal_structure.suppressed( i ) );
            crystal_structure.save_cif( append_to_file_name( input_file_name, "_" + size_t2string( iSymmOp ) + "_inverse_2" ) );
        }
    MACRO_END_GAME

    try // Random numbers.
    {
        RandomNumberGenerator_integer RNG;
        for ( size_t i( 0 ); i != 10; ++i )
        {
            std::string result;
            for ( size_t j( 0 ); j != 10; ++j )
            {
                if ( is_even( RNG.next_number( 1, 20 ) ) )
                    result += 'A';
                else
                    result += 'B';
            }
            std::cout << result << std::endl;
        }
    MACRO_END_GAME

    try // MC alkanes.
    {
        size_t n = 20; // Number of carbon atoms;
        size_t nTrials( 10000000 );
        size_t noverlap( 0 );
        double exclusion_distance = 3.0;
        RandomNumberGenerator_integer rng;
        if ( false )
        {
            for ( size_t iTrial( 0 ); iTrial != nTrials; ++iTrial )
            {
                if ( ( iTrial % 1000 ) == 0 )
                    std::cout << "Trial " << iTrial << std::endl;
                std::vector< Angle > torsion_angles;
                for ( size_t j( 0 ); j != n-3; ++j )
                {
                    size_t iRandom = rng.next_number( 0, 2 );
                    torsion_angles.push_back( iRandom * Angle::from_degrees( 120.0 ) );
                }
                std::vector< Vector3D > coordinates = build_alkane( n, torsion_angles );
                if ( ( iTrial % 1000 ) == 0 )
                    save_as_xyz( coordinates, FileName( "C:\\Users\\jacco\\Documents\\MeltingPoints\\trial_" + size_t2string( iTrial ) + ".xyz" ) );
                if ( there_is_overlap( coordinates, exclusion_distance ) )
                    ++noverlap;
            }
        }
        else
        {
            // Much faster version, tests for overlap with each newly added atom and terminates as soon as overlap is detected.
            double exclusion_distance2 = square( exclusion_distance );
            Angle C_C_C = Angle::from_degrees( 113.5 );
            std::vector< Vector3D > first_three;
            first_three.push_back( Vector3D() );
            first_three.push_back( Vector3D( 1.54, 0.0, 0.0) );
            first_three.push_back( Vector3D( 1.54 + ( C_C_C - Angle::from_degrees( 90.0 ) ).sine() * 1.54, ( C_C_C - Angle::from_degrees( 90.0 ) ).cosine() * 1.54, 0.0 ) );
            for ( size_t iTrial( 0 ); iTrial != nTrials; ++iTrial )
            {
                if ( ( iTrial % 100000 ) == 0 )
                    std::cout << "Trial " << iTrial << std::endl;
                std::vector< Vector3D > coordinates = first_three;
                coordinates.reserve( n );
                for ( size_t i( 3 ); i != n; ++i )
                {
                    Vector3D new_coordinate = coordinates[i-1] + 1.54 * normalised_vector( coordinates[i-2] - coordinates[i-3] );
                    NormalisedVector3D axis = normalised_vector( coordinates[i-2] - coordinates[i-1] );
                    new_coordinate = rotate_point_about_axis( new_coordinate, coordinates[i-1], axis, rng.next_number( 0, 2 ) * Angle::from_degrees( 120.0 ) );
                    coordinates.push_back( new_coordinate );
                    // Find all interatomc distances, except those between nearest neightbours (because they are bonded).
                    bool overlap_found( false );
                    for ( size_t j( 0 ); j != i-2; ++j )
                    {
                        if ( ( coordinates[i] - coordinates[j] ).norm2() < exclusion_distance2 )
                        {
                            ++noverlap;
                            overlap_found = true;
                            break;
                        }
                    }
                    if ( overlap_found )
                        break;
                }
            }
        }
        double fraction_without_overlap = ( static_cast<double>(nTrials) - static_cast<double>(noverlap) ) / static_cast<double>(nTrials);
        std::cout << "Fraction without overlap = " << fraction_without_overlap << std::endl;
        std::cout << "Surviving conformations  = fraction * 3^(n-3)/2 = " << fraction_without_overlap * std::pow( 3.0, n-3 ) / 2.0 << std::endl;
        std::cout << "ln( fraction * 3^(n-3)/2 ) = " << ln( fraction_without_overlap * std::pow( 3.0, n-3 ) / 2.0 ) << std::endl;
    MACRO_END_GAME

    try // Analyse GRACE results.
    {
        TextFileReader_2 input_file( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\ucfx_pairs.txt" ) );
        if ( input_file.size() != 17 )
            throw std::runtime_error( "Expected 17 pairs." );
        TextFileReader_2 RESULTS_ucfr_file( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE\\RESULTS_ucfr.txt" ) );
        TextFileReader_2 RESULTS_ucfx_file( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE\\RESULTS_ucfx.txt" ) );
        std::vector< double > min_RMSCDs;
        std::vector< double > max_RMSCDs;
        std::vector< double > energy_differences;
        for ( size_t iLine( 0 ); iLine != input_file.size(); ++iLine )
        {
        //    std::cout << "iLine = " << iLine << std::endl;
//            if ( iLine != 7 ) // IPRPOL has problems with floating axes
//                continue;
            if ( iLine == 16 ) // ZUHRID: space groups are different
                continue;
            std::string line = input_file.line( iLine );
            std::vector< std::string > words = split( line );
            if ( words.size() != 4 )
                throw std::runtime_error( "Expected three items per line." );
            std::string identifier_1( words[1] );
            std::string identifier_2( words[2] );
            std::cout << identifier_1 << " " << identifier_2 << std::endl;
            size_t natoms = string2integer( words[3] );

            CrystalStructure crystal_structure_1;
            CrystalStructure crystal_structure_2;
            double RMSCD( -1.0 );
            size_t iPos( 0 );
            double energy_1( 10000.0 );
            double energy_2( 10000.0 );
            size_t shift_steps( 2 );
            Mapping mapping;
            SymmetryOperator symmetry_operator;
            std::vector< Vector3D > translations;
            bool identifier_1_has_the_lowest_RMSCD;
            
            bool as_in_paper( false );
            bool as_by_jvds( true );
            
            // RMSCD between the two exp. structures before energy minimisation.
            read_cif( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE", identifier_1 + "_Hnorm", "cif" ), crystal_structure_1 );
            read_cif( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE", identifier_2 + "_Hnorm", "cif" ), crystal_structure_2 );
            crystal_structure_1.remove_H_and_D();
            crystal_structure_2.remove_H_and_D();

            crystal_structure_1.perceive_molecules( true );
            crystal_structure_2.perceive_molecules( true );
            map( crystal_structure_2, crystal_structure_1, shift_steps, mapping, symmetry_operator, translations, false, true );
            
//            crystal_structure_2.apply_map( mapping, symmetry_operator, translations );
//            RMSCD = root_mean_square_Cartesian_displacement( crystal_structure_1, crystal_structure_2, false );
//            std::cout << identifier_1 << " " << identifier_2 << " RMSCD = " << RMSCD << std::endl;

            if ( as_in_paper )
            {
            read_cif( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE", identifier_1 + "_avguc_Hnorm", "cif" ), crystal_structure_1 );
            read_cif( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE", identifier_2 + "_avguc_Hnorm", "cif" ), crystal_structure_2 );
            crystal_structure_1.remove_H_and_D();
            crystal_structure_2.remove_H_and_D();
            crystal_structure_1.perceive_molecules( true );
            crystal_structure_2.perceive_molecules( true );
            crystal_structure_2.apply_map( mapping, symmetry_operator, translations );
            RMSCD = root_mean_square_Cartesian_displacement( crystal_structure_1, crystal_structure_2, false );
            std::cout << identifier_1 << "_avguc " << identifier_2 << "_avguc RMSCD = " << RMSCD << std::endl;
            }

//            // RMSCD between the two structures energy minimised with the unit cell fixed.
//            read_cif( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE", identifier_1 + "_Hnorm_mi_ucfx", "cif" ), crystal_structure_1 );
//            read_cif( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE", identifier_2 + "_Hnorm_mi_ucfx", "cif" ), crystal_structure_2 );
//            crystal_structure_1.remove_H_and_D();
//            crystal_structure_2.remove_H_and_D();
//            crystal_structure_1.perceive_molecules( true );
//            crystal_structure_2.perceive_molecules( true );
//            crystal_structure_2.apply_map( mapping, symmetry_operator, translations );
//            RMSCD = root_mean_square_Cartesian_displacement( crystal_structure_1, crystal_structure_2, false );
//            std::cout << identifier_1 << "_ucfx " << identifier_2 << "_ucfx RMSCD = " << RMSCD << std::endl;
//            iPos = RESULTS_ucfx_file.find_whole_word( identifier_1 + "_Hnorm" );
//            line = RESULTS_ucfx_file.line( iPos );
//            words = split( line );
//            energy_1 = string2double( words[2] );
//            iPos = RESULTS_ucfx_file.find_whole_word( identifier_2 + "_Hnorm" );
//            line = RESULTS_ucfx_file.line( iPos );
//            words = split( line );
//            energy_2 = string2double( words[2] );
//            std::cout << "Energy = " << ( energy_2 - energy_1 ) * natoms << std::endl;

            if ( as_in_paper )
            {
            read_cif( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE", identifier_1 + "_avguc_Hnorm_mi_ucfx", "cif" ), crystal_structure_1 );
            read_cif( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE", identifier_2 + "_avguc_Hnorm_mi_ucfx", "cif" ), crystal_structure_2 );
            crystal_structure_1.remove_H_and_D();
            crystal_structure_2.remove_H_and_D();
            crystal_structure_1.perceive_molecules( true );
            crystal_structure_2.perceive_molecules( true );
            crystal_structure_2.apply_map( mapping, symmetry_operator, translations );
            RMSCD = root_mean_square_Cartesian_displacement( crystal_structure_1, crystal_structure_2, false );
            std::cout << identifier_1 << "_avguc_ucfx " << identifier_2 << "_avguc_ucfx RMSCD = " << RMSCD << std::endl;
            iPos = RESULTS_ucfx_file.find_whole_word( identifier_1 + "_avguc_Hnorm" );
            line = RESULTS_ucfx_file.line( iPos );
            words = split( line );
            energy_1 = string2double( words[2] );
            iPos = RESULTS_ucfx_file.find_whole_word( identifier_2 + "_avguc_Hnorm" );
            line = RESULTS_ucfx_file.line( iPos );
            words = split( line );
            energy_2 = string2double( words[2] );
            std::cout << identifier_2 << "_avguc_ucfx - " << identifier_1 << "_avguc_ucfx Energy = " << ( energy_2 - energy_1 ) * natoms << std::endl;
            }

            // RMSCDs of the two structures energy minimised with the unit cell fixed.
//            read_cif( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE", identifier_1 + "_Hnorm", "cif" ), crystal_structure_1 );
//            read_cif( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE", identifier_1 + "_Hnorm_mi_ucfx", "cif" ), crystal_structure_2 );
//            crystal_structure_1.remove_H_and_D();
//            crystal_structure_2.remove_H_and_D();
//            RMSCD = root_mean_square_Cartesian_displacement( crystal_structure_1, crystal_structure_2, false );
//            std::cout << identifier_1 << " " << identifier_1 << "_ucfx RMSCD = " << RMSCD << std::endl;

//            read_cif( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE", identifier_2 + "_Hnorm", "cif" ), crystal_structure_1 );
//            read_cif( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE", identifier_2 + "_Hnorm_mi_ucfx", "cif" ), crystal_structure_2 );
//            crystal_structure_1.remove_H_and_D();
//            crystal_structure_2.remove_H_and_D();
//            RMSCD = root_mean_square_Cartesian_displacement( crystal_structure_1, crystal_structure_2, false );
//            std::cout << identifier_2 << " " << identifier_2 << "_ucfx RMSCD = " << RMSCD << std::endl;

            if ( as_in_paper )
            {
            read_cif( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE", identifier_1 + "_avguc_Hnorm", "cif" ), crystal_structure_1 );
            read_cif( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE", identifier_1 + "_avguc_Hnorm_mi_ucfx", "cif" ), crystal_structure_2 );
            crystal_structure_1.remove_H_and_D();
            crystal_structure_2.remove_H_and_D();
            RMSCD = root_mean_square_Cartesian_displacement( crystal_structure_1, crystal_structure_2, false );
            std::cout << identifier_1 << "_avguc " << identifier_1 << "_avguc_ucfx RMSCD = " << RMSCD << std::endl;

            read_cif( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE", identifier_2 + "_avguc_Hnorm", "cif" ), crystal_structure_1 );
            read_cif( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE", identifier_2 + "_avguc_Hnorm_mi_ucfx", "cif" ), crystal_structure_2 );
            crystal_structure_1.remove_H_and_D();
            crystal_structure_2.remove_H_and_D();
            RMSCD = root_mean_square_Cartesian_displacement( crystal_structure_1, crystal_structure_2, false );
            std::cout << identifier_2 << "_avguc " << identifier_2 << "_avguc_ucfx RMSCD = " << RMSCD << std::endl;
            }

            if ( as_in_paper )
            {
            // RMSCD between the two structures energy minimised with the unit cell free.
            read_cif( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE", identifier_1 + "_Hnorm_mi_ucfr", "cif" ), crystal_structure_1 );
            read_cif( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE", identifier_2 + "_Hnorm_mi_ucfr", "cif" ), crystal_structure_2 );
            crystal_structure_1.remove_H_and_D();
            crystal_structure_2.remove_H_and_D();
            crystal_structure_1.perceive_molecules( true );
            crystal_structure_2.perceive_molecules( true );
            crystal_structure_2.apply_map( mapping, symmetry_operator, translations );
            RMSCD = root_mean_square_Cartesian_displacement( crystal_structure_1, crystal_structure_2, false );
            std::cout << identifier_1 << "_ucfr " << identifier_2 << "_ucfr RMSCD = " << RMSCD << std::endl;
            }
            
            if ( false )
            {
            iPos = RESULTS_ucfr_file.find_whole_word( identifier_1 + "_Hnorm" );
            line = RESULTS_ucfr_file.line( iPos );
            words = split( line );
            energy_1 = string2double( words[2] );
            iPos = RESULTS_ucfr_file.find_whole_word( identifier_2 + "_Hnorm" );
            line = RESULTS_ucfr_file.line( iPos );
            words = split( line );
            energy_2 = string2double( words[2] );
            std::cout << identifier_2 << "_ucfr - " << identifier_1 << "_ucfr Energy = " << ( energy_2 - energy_1 ) * natoms << std::endl;
            }

            if ( as_in_paper )
            {
            read_cif( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE", identifier_1 + "_avguc_Hnorm_mi_ucfr", "cif" ), crystal_structure_1 );
            read_cif( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE", identifier_2 + "_avguc_Hnorm_mi_ucfr", "cif" ), crystal_structure_2 );
            crystal_structure_1.remove_H_and_D();
            crystal_structure_2.remove_H_and_D();
            crystal_structure_1.perceive_molecules( true );
            crystal_structure_2.perceive_molecules( true );
            crystal_structure_2.apply_map( mapping, symmetry_operator, translations );
            RMSCD = root_mean_square_Cartesian_displacement( crystal_structure_1, crystal_structure_2, false );
            std::cout << identifier_1 << "_avguc_ucfr " << identifier_2 << "_avguc_ucfr RMSCD = " << RMSCD << std::endl;
            iPos = RESULTS_ucfr_file.find_whole_word( identifier_1 + "_avguc_Hnorm" );
            line = RESULTS_ucfr_file.line( iPos );
            words = split( line );
            energy_1 = string2double( words[2] );
            iPos = RESULTS_ucfr_file.find_whole_word( identifier_2 + "_avguc_Hnorm" );
            line = RESULTS_ucfr_file.line( iPos );
            words = split( line );
            energy_2 = string2double( words[2] );
            std::cout << identifier_2 << "_avguc_ucfr - " << identifier_1 << "_avguc_ucfr Energy = " << ( energy_2 - energy_1 ) * natoms << std::endl;
            }

            double RMSCD_identifier_1;
            double RMSCD_identifier_2;
            if ( as_by_jvds )
            {
            // RMSCDs of the two structures energy minimised with the unit cell free.
            read_cif( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE", identifier_1 + "_Hnorm", "cif" ), crystal_structure_1 );
            read_cif( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE", identifier_1 + "_Hnorm_mi_ucfr", "cif" ), crystal_structure_2 );
            crystal_structure_1.remove_H_and_D();
            crystal_structure_2.remove_H_and_D();
            RMSCD_identifier_1 = root_mean_square_Cartesian_displacement( crystal_structure_1, crystal_structure_2, false );
         //   std::cout << identifier_1 << " " << identifier_1 << "_ucfr RMSCD = " << RMSCD << std::endl;

            read_cif( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE", identifier_2 + "_Hnorm", "cif" ), crystal_structure_1 );
            read_cif( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE", identifier_2 + "_Hnorm_mi_ucfr", "cif" ), crystal_structure_2 );
            crystal_structure_1.remove_H_and_D();
            crystal_structure_2.remove_H_and_D();
            RMSCD_identifier_2 = root_mean_square_Cartesian_displacement( crystal_structure_1, crystal_structure_2, false );
            identifier_1_has_the_lowest_RMSCD = ( RMSCD_identifier_1 < RMSCD_identifier_2 );
         //   std::cout << identifier_2 << " " << identifier_2 << "_ucfr RMSCD = " << RMSCD << std::endl;
            }

            if ( as_by_jvds )
            {
            iPos = RESULTS_ucfr_file.find_whole_word( identifier_1 + "_Hnorm" );
            line = RESULTS_ucfr_file.line( iPos );
            words = split( line );
            energy_1 = string2double( words[2] );
            iPos = RESULTS_ucfr_file.find_whole_word( identifier_2 + "_Hnorm" );
            line = RESULTS_ucfr_file.line( iPos );
            words = split( line );
            energy_2 = string2double( words[2] );
            if ( ! identifier_1_has_the_lowest_RMSCD )
            {
                std::swap( identifier_1, identifier_2 );
                std::swap( RMSCD_identifier_1, RMSCD_identifier_2 );
                std::swap( energy_1, energy_2 );
            }
            if ( identifier_1.substr( identifier_1.length()-6, 6 ) == "_trans" )
                identifier_1 = identifier_1.substr( 0, identifier_1.length()-6 );
            if ( identifier_2.substr( identifier_2.length()-6, 6 ) == "_trans" )
                identifier_2 = identifier_2.substr( 0, identifier_2.length()-6 );
            std::cout << pad( identifier_1, 8, ' ' ) << " RMSCD = " << RMSCD_identifier_1 << std::endl;
            std::cout << pad( identifier_2, 8, ' ' ) << " RMSCD = " << RMSCD_identifier_2 << std::endl;
            min_RMSCDs.push_back( RMSCD_identifier_1 );
            max_RMSCDs.push_back( RMSCD_identifier_2 );
            energy_differences.push_back( ( energy_2 - energy_1 ) * natoms );
            if ( energy_1 < energy_2 )
            {
                std::cout << "#####+++++++++++++++++++" << std::endl;
                std::cout << identifier_1 << " has the lowest  energy, it is " << ( energy_2 - energy_1 ) * natoms << " kcal/mol lower" << std::endl;
            }
            else
            {
                std::cout << "#####- - - - - - - - - " << std::endl;
                std::cout << identifier_1 << " has the highest energy, it is " << ( energy_1 - energy_2 ) * natoms << " kcal/mol higher" << std::endl;
            }
            }

            if ( as_in_paper )
            {
            read_cif( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE", identifier_1 + "_avguc_Hnorm", "cif" ), crystal_structure_1 );
            read_cif( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE", identifier_1 + "_avguc_Hnorm_mi_ucfr", "cif" ), crystal_structure_2 );
            crystal_structure_1.remove_H_and_D();
            crystal_structure_2.remove_H_and_D();
            RMSCD = root_mean_square_Cartesian_displacement( crystal_structure_1, crystal_structure_2, false );
            std::cout << identifier_1 << "_avguc " << identifier_1 << "_avguc_ucfr RMSCD = " << RMSCD << std::endl;

            read_cif( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE", identifier_2 + "_avguc_Hnorm", "cif" ), crystal_structure_1 );
            read_cif( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE", identifier_2 + "_avguc_Hnorm_mi_ucfr", "cif" ), crystal_structure_2 );
            crystal_structure_1.remove_H_and_D();
            crystal_structure_2.remove_H_and_D();
            RMSCD = root_mean_square_Cartesian_displacement( crystal_structure_1, crystal_structure_2, false );
            std::cout << identifier_2 << "_avguc " << identifier_2 << "_avguc_ucfr RMSCD = " << RMSCD << std::endl;
            }
            std::cout << std::endl;

        }
        // For ZUHRID, we cannot calculate the figures from the paper, only the
        // figures as done by JvdS

        if ( true )
        {
        std::string identifier_1 = "ZUHRID";
        std::string identifier_2 = "ZUHRID02_avguc";
        // RMSCDs of the two structures energy minimised with the unit cell free.
        CrystalStructure crystal_structure_1;
        CrystalStructure crystal_structure_2;
        read_cif( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE", identifier_1 + "_Hnorm", "cif" ), crystal_structure_1 );
        read_cif( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE", identifier_1 + "_Hnorm_mi_ucfr", "cif" ), crystal_structure_2 );
        crystal_structure_1.remove_H_and_D();
        crystal_structure_2.remove_H_and_D();
        double RMSCD_identifier_1 = root_mean_square_Cartesian_displacement( crystal_structure_1, crystal_structure_2, false );
        read_cif( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE", identifier_2 + "_Hnorm", "cif" ), crystal_structure_1 );
        read_cif( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE", identifier_2 + "_Hnorm_mi_ucfr", "cif" ), crystal_structure_2 );
        crystal_structure_1.remove_H_and_D();
        crystal_structure_2.remove_H_and_D();
        double RMSCD_identifier_2 = root_mean_square_Cartesian_displacement( crystal_structure_1, crystal_structure_2, false );
        bool identifier_1_has_the_lowest_RMSCD = ( RMSCD_identifier_1 < RMSCD_identifier_2 );
        size_t iPos = RESULTS_ucfr_file.find_whole_word( identifier_1 + "_Hnorm" );
        std::string line = RESULTS_ucfr_file.line( iPos );
        std::vector< std::string > words = split( line );
        double energy_1 = string2double( words[2] );
        iPos = RESULTS_ucfr_file.find_whole_word( identifier_2 + "_Hnorm" );
        line = RESULTS_ucfr_file.line( iPos );
        words = split( line );
        double energy_2 = string2double( words[2] );
        if ( ! identifier_1_has_the_lowest_RMSCD )
        {
            std::swap( identifier_1, identifier_2 );
            std::swap( RMSCD_identifier_1, RMSCD_identifier_2 );
            std::swap( energy_1, energy_2 );
        }
        std::cout << pad( identifier_1, 8, ' ' ) << " RMSCD = " << RMSCD_identifier_1 << std::endl;
        std::cout << pad( identifier_2, 8, ' ' ) << " RMSCD = " << RMSCD_identifier_2 << std::endl;
        min_RMSCDs.push_back( RMSCD_identifier_1 );
        max_RMSCDs.push_back( RMSCD_identifier_2 );
        energy_differences.push_back( ( energy_2 - energy_1 ) * 30 );
        if ( energy_1 < energy_2 )
        {
            std::cout << "#####+++++++++++++++++++" << std::endl;
            std::cout << identifier_1 << " has the lowest  energy, it is " << ( energy_2 - energy_1 ) * 30 << " kcal/mol lower" << std::endl;
        }
        else
        {
            std::cout << "#####- - - - - - - - - " << std::endl;
            std::cout << identifier_1 << " has the highest energy, it is " << ( energy_1 - energy_2 ) * 30 << " kcal/mol higher" << std::endl;
        }
        }
        std::cout << std::endl;
        for ( size_t i( 0 ); i != min_RMSCDs.size(); ++i )
        {
            std::cout << min_RMSCDs[i] << " " << max_RMSCDs[i] << std::endl;
        }
        std::cout << std::endl;
        for ( size_t i( 0 ); i != min_RMSCDs.size(); ++i )
        {
            std::cout << max_RMSCDs[i] << " " << energy_differences[i] << std::endl;
        }
    MACRO_END_GAME

    try // ss-NMR.
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

    try // Graeme's 50 ESI.
    {
        FileName file_name_cif( "Graeme_50_ESI.txt" );
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

    try // Generate FileList.txt .
    {
        if ( argc != 2 )
            throw std::runtime_error( "Please give the name of a directory." );
        TextFileWriter text_file_writer( FileName( argv[ 1 ], "FileList", "txt" ) );
        for ( size_t i( 0 ); i != 265; ++i )
            text_file_writer.write_line( "structure_" + size_t2string( i+1, 6, '0') + ".cif" );
    MACRO_END_GAME

    try // The original cif must have atom labels with "A" and "B", e.g. "C7A" and C7B".
    {
        if ( argc != 10 )
            throw std::runtime_error( "Usage: <name>.cif u v w y/n y/n y/n true/false (write cif) true/false (calculate XRPD)" );
        FileName input_file_name( argv[ 1 ] );
        size_t u( string2integer( argv[ 2 ] ) );
        size_t v( string2integer( argv[ 3 ] ) );
        size_t w( string2integer( argv[ 4 ] ) );
        std::vector< bool > is_correlated;
        is_correlated.push_back( string2bool( argv[ 5 ] ) );
        is_correlated.push_back( string2bool( argv[ 6 ] ) );
        is_correlated.push_back( string2bool( argv[ 7 ] ) );
        bool write_cif( string2bool( argv[ 8 ] ) );
        bool calculate_XRPD( string2bool( argv[ 9 ] ) );
        size_t ncorrelated_directions( 0 );
        for ( size_t i( 0 ); i != is_correlated.size(); ++i )
            if ( is_correlated[i] )
                ++ncorrelated_directions;
        if ( ncorrelated_directions == 3 )
            throw std::runtime_error( "The disorder is correlated in all three directions, there is therefore no disorder." );
        if ( ( u == 0 ) || ( v == 0 ) || ( w == 0 ) )
            throw std::runtime_error( "Each dimension must be at least one." );
        if ( is_correlated[0] )
        {
            if ( is_odd( u ) )
                throw std::runtime_error( "Error: correlation, but u is odd." );
            else // is_even()
            {
                if ( u != 2 )
                {
                    std::cout << "Warning: correlation and u is even but not 2, u was changed to 2." << std::endl;
                    u = 2;
                }
            }
        }
        if ( is_correlated[1] )
        {
            if ( is_odd( v ) )
                throw std::runtime_error( "Error: correlation, but v is odd." );
            else // is_even()
            {
                if ( v != 2 )
                {
                    std::cout << "Warning: correlation and v is even but not 2, v was changed to 2." << std::endl;
                    v = 2;
                }
            }
        }
        if ( is_correlated[2] )
        {
            if ( is_odd( w ) )
                throw std::runtime_error( "Error: correlation, but w is odd." );
            else // is_even()
            {
                if ( w != 2 )
                {
                    std::cout << "Warning: correlation and w is even but not 2, w was changed to 2." << std::endl;
                    w = 2;
                }
            }
        }
        CrystalStructure crystal_structure_org;
        read_cif( input_file_name, crystal_structure_org );
        crystal_structure_org.apply_space_group_symmetry( false );
        RandomNumberGenerator_double rng;
        if ( false )
        {
            for ( size_t i( 0 ); i != 10; ++i )
            {
                double x = rng.next_number();
                double y = rng.next_number();
                double z = rng.next_number();
                Vector3D new_position( x, y, z );
                bool collides_with_existing( false );
                do
                {
                    for ( size_t j( 0 ); j != crystal_structure_org.natoms(); ++j )
                    {
                        if ( crystal_structure_org.crystal_lattice().shortest_distance( crystal_structure_org.atom( j ).position(), new_position ) < 1.5 )
                        {
                            collides_with_existing = true;
                            break;
                        }
                    }
                } while ( collides_with_existing );
                Atom new_atom( Element( "C" ), new_position, "C" + size_t2string( i ) );
                crystal_structure_org.add_atom( new_atom );
            }
        }
        CrystalLattice crystal_lattice = crystal_structure_org.crystal_lattice();
        CrystalStructure crystal_structure;
        crystal_structure.set_space_group( SpaceGroup() );
        CrystalLattice new_crystal_lattice( crystal_lattice.a() * u, crystal_lattice.b() * v, crystal_lattice.c() * w, crystal_lattice.alpha(), crystal_lattice.beta(), crystal_lattice.gamma() );
        crystal_structure.set_crystal_lattice( new_crystal_lattice );
        crystal_structure.reserve_natoms( crystal_structure_org.natoms() * u * v * w );
        std::vector< std::vector< std::vector < bool > > > configuration( u, std::vector< std::vector < bool > >( v, std::vector < bool >( w, false ) ) );
        if ( ncorrelated_directions == 0 )
        {
            for ( size_t ix( 0 ); ix != u; ++ix )
            {
                for ( size_t iy( 0 ); iy != v; ++iy )
                {
                    for ( size_t iz( 0 ); iz != w; ++iz )
                    {
                        if ( rng.next_number() < 0.5 )
                            configuration[ ix ][ iy ][ iz ] = true;
                    }
                }
            }
        }
        else if ( ncorrelated_directions == 1 )
        {
            if ( is_correlated[0] )
            {
                for ( size_t iy( 0 ); iy != v; ++iy )
                {
                    for ( size_t iz( 0 ); iz != w; ++iz )
                    {
                        bool odd_is_true( rng.next_number() < 0.5 );
                        for ( size_t ix( 0 ); ix != u; ++ix )
                        {
                            configuration[ ix ][ iy ][ iz ] = odd_is_true;
                            odd_is_true = ( ! odd_is_true );
                        }
                    }
                }
            }
            else if ( is_correlated[1] )
            {
                for ( size_t ix( 0 ); ix != u; ++ix )
                {
                    for ( size_t iz( 0 ); iz != w; ++iz )
                    {
                        bool odd_is_true( rng.next_number() < 0.5 );
                        for ( size_t iy( 0 ); iy != v; ++iy )
                        {
                            configuration[ ix ][ iy ][ iz ] = odd_is_true;
                            odd_is_true = ( ! odd_is_true );
                        }
                    }
                }
            }
            else // is_correlated[2]
            {
                for ( size_t ix( 0 ); ix != u; ++ix )
                {
                    for ( size_t iy( 0 ); iy != v; ++iy )
                    {
                        bool odd_is_true( rng.next_number() < 0.5 );
                        for ( size_t iz( 0 ); iz != w; ++iz )
                        {
                            configuration[ ix ][ iy ][ iz ] = odd_is_true;
                            odd_is_true = ( ! odd_is_true );
                        }
                    }
                }
            }
        }
        else // ( ncorrelated_directions == 2 )
        {
            if ( ! is_correlated[0] )
            {
                for ( size_t ix( 0 ); ix != u; ++ix )
                {
                    bool odd_is_true( rng.next_number() < 0.5 );
                    for ( size_t iy( 0 ); iy != v; ++iy )
                    {
                        for ( size_t iz( 0 ); iz != w; ++iz )
                        {
                            if ( is_even( iy + iz ) )
                                configuration[ ix ][ iy ][ iz ] = odd_is_true;
                            else
                                configuration[ ix ][ iy ][ iz ] = ( ! odd_is_true );
                        }
                    }
                }
            }
            else if ( ! is_correlated[1] )
            {
                for ( size_t iy( 0 ); iy != v; ++iy )
                {
                    bool odd_is_true( rng.next_number() < 0.5 );
                    for ( size_t ix( 0 ); ix != u; ++ix )
                    {
                        for ( size_t iz( 0 ); iz != w; ++iz )
                        {
                            if ( is_even( ix + iz ) )
                                configuration[ ix ][ iy ][ iz ] = odd_is_true;
                            else
                                configuration[ ix ][ iy ][ iz ] = ( ! odd_is_true );
                        }
                    }
                }
            }
            else // ( ! is_correlated[2] )
            {
                for ( size_t iz( 0 ); iz != w; ++iz )
                {
                    bool odd_is_true( rng.next_number() < 0.5 );
                    for ( size_t ix( 0 ); ix != u; ++ix )
                    {
                        for ( size_t iy( 0 ); iy != v; ++iy )
                        {
                            if ( is_even( ix + iy ) )
                                configuration[ ix ][ iy ][ iz ] = odd_is_true;
                            else
                                configuration[ ix ][ iy ][ iz ] = ( ! odd_is_true );
                        }
                    }
                }
            }
        }
        for ( size_t ix( 0 ); ix != u; ++ix )
        {
            for ( size_t iy( 0 ); iy != v; ++iy )
            {
                for ( size_t iz( 0 ); iz != w; ++iz )
                {
                    bool A_or_B = configuration[ ix ][ iy ][ iz ];
                    for ( size_t iAtom( 0 ); iAtom != crystal_structure_org.natoms(); ++iAtom )
                    {
                        Atom atom( crystal_structure_org.atom( iAtom ) );
                        atom.reset_ADPs_type();
                        std::string last_character = atom.label();
                        if ( last_character.length() == 0 )
                            throw std::runtime_error( "Atom has no label." );
                        last_character = last_character.substr( last_character.length()-1, 1 );
                        if ( ( last_character == "A" ) && ( ! A_or_B ) )
                            continue;
                        else if ( ( last_character == "B" ) && A_or_B )
                            continue;
                        if ( ( last_character == "A" ) || ( last_character == "B" ) )
                        {
                            if ( ! nearly_equal( atom.occupancy(), 0.5 ) )
                                std::cout << "Warning: occupancy of disordered atom was not 0.5." << std::endl;
                        }
                        else
                        {
                            if ( ! nearly_equal( atom.occupancy(), 1.0 ) )
                                std::cout << "Warning: occupancy of ordered atom was not 1.0." << std::endl;
                        }
                        atom.set_occupancy( 1.0 );
                        Vector3D position = atom.position();
                        position = Vector3D( ( position.x() + ix ) / u, ( position.y() + iy ) / v, ( position.z() + iz ) / w );
                        atom.set_position( position );
                        atom.set_label( atom.label() + "_" + size_t2string( ix ) + "_" + size_t2string( iy ) + "_" + size_t2string( iz ) );
                        crystal_structure.add_atom( atom );
                    }
                }
            }
        }
        std::string file_appendix( "_" + bool2string( is_correlated[0] ) + bool2string( is_correlated[1] ) + bool2string( is_correlated[2] ) + "_" + size_t2string( u ) + "x" + size_t2string( v ) + "x" + size_t2string( w ) );
        if( write_cif )
            crystal_structure.save_cif( append_to_file_name( input_file_name, file_appendix ) );
        if ( calculate_XRPD )
        {
            crystal_structure.apply_space_group_symmetry();
            std::cout << "Now calculating powder pattern... " << std::endl;
            PowderPatternCalculator powder_pattern_calculator( crystal_structure );
            Angle two_theta_start( 1.0, Angle::DEGREES );
            Angle two_theta_end(  35.0, Angle::DEGREES );
            Angle two_theta_step( 0.01, Angle::DEGREES );
            powder_pattern_calculator.set_two_theta_start( two_theta_start );
            powder_pattern_calculator.set_two_theta_end( two_theta_end );
            powder_pattern_calculator.set_two_theta_step( two_theta_step );
            powder_pattern_calculator.set_FWHM( 0.10 );
            PowderPattern powder_pattern;
            powder_pattern_calculator.calculate( powder_pattern );
            powder_pattern.save_xye( append_to_file_name( replace_extension( input_file_name, "xye" ), file_appendix ), true );
        }
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
        // Cartesian coordinates in cell 0,0,0.
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
        crystal_structure.save_cif( FileName( "disorder_04.cif" ) );

        if ( false )
        {
            crystal_structure.apply_space_group_symmetry();
            std::cout << "Now calculating powder pattern... " << std::endl;
            PowderPatternCalculator powder_pattern_calculator( crystal_structure );
            Angle two_theta_start( 1.0, Angle::DEGREES );
            Angle two_theta_end(  40.0, Angle::DEGREES );
            Angle two_theta_step( 0.015, Angle::DEGREES );
            powder_pattern_calculator.set_two_theta_start( two_theta_start );
            powder_pattern_calculator.set_two_theta_end( two_theta_end );
            powder_pattern_calculator.set_two_theta_step( two_theta_step );
            powder_pattern_calculator.set_FWHM( 0.20 );
            PowderPattern powder_pattern;
            powder_pattern_calculator.calculate( powder_pattern );
            powder_pattern.save_xye( FileName( "disorder_02.xye" ), true );
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
        powder_pattern_calculator.set_two_theta_start( two_theta_start );
        powder_pattern_calculator.set_two_theta_end( two_theta_end );
        powder_pattern_calculator.set_two_theta_step( two_theta_step );
        powder_pattern_calculator.set_FWHM( 0.20 );
        PowderPattern powder_pattern;
        powder_pattern_calculator.calculate( powder_pattern );
        powder_pattern.save_xye( FileName( "C:\\Users\\jacco\\Documents\\disorder_01.xye" ), true );
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

    try // Change bond length.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        {
        std::string origin_atom_label( "C21" );
        std::string neighbour_atom_label( "H30" );
        Vector3D origin_atom_frac = crystal_structure.atom( crystal_structure.find_label( origin_atom_label ) ).position();
        Vector3D neighbour_atom_frac = crystal_structure.atom( crystal_structure.find_label( neighbour_atom_label ) ).position();
        // Find shortest distance between the two taking periodicity and space-group symmetry into account.
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
        double target_bond_length( 1.089 );
        difference_cart *= target_bond_length / distance;
        difference_frac = crystal_structure.crystal_lattice().orthogonal_to_fractional( difference_cart );
        Atom new_atom( crystal_structure.atom( crystal_structure.find_label( neighbour_atom_label ) ) );
        new_atom.set_position( origin_atom_frac + difference_frac );
    //    new_atom.set_element( Element( "F" ) );
        crystal_structure.set_atom( crystal_structure.find_label( neighbour_atom_label ), new_atom );
        }
        crystal_structure.save_cif( append_to_file_name( input_file_name, "_blc" ) );
    MACRO_END_GAME

    try // Collapse Z'=2 into disordered Z'=1.
    {
        if ( argc != 5 )
            throw std::runtime_error( "Usage: u v w <name>.cif" );
        size_t u( string2integer( argv[ 1 ] ) );
        size_t v( string2integer( argv[ 2 ] ) );
        size_t w( string2integer( argv[ 3 ] ) );
        FileName input_file_name( argv[ 4 ] );
        CrystalStructure crystal_structure;
        read_cif( input_file_name, crystal_structure );
        if ( crystal_structure.space_group().nsymmetry_operators() != 1 )
            std::cout << "Warning: the space group is not P1." << std::endl;
        for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
        {
            Atom atom( crystal_structure.atom( i ) );
            Vector3D position = atom.position();
            position.set_x( u * position.x() );
            position.set_y( v * position.y() );
            position.set_z( w * position.z() );
            atom.set_position( adjust_for_translations( position ) );
            crystal_structure.set_atom( i, atom );
        }
        CrystalLattice crystal_lattice( crystal_structure.crystal_lattice().a() / u,
                                        crystal_structure.crystal_lattice().b() / v,
                                        crystal_structure.crystal_lattice().c() / w,
                                        crystal_structure.crystal_lattice().alpha(),
                                        crystal_structure.crystal_lattice().beta(),
                                        crystal_structure.crystal_lattice().gamma() );
        crystal_structure.set_crystal_lattice( crystal_lattice );
        crystal_structure.save_cif( append_to_file_name( input_file_name, "_collapsed" ) );
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

    try // Repair XRPD pattern extracted from a .png.
    {
        MACRO_ONE_XYEFILENAME_AS_ARGUMENT
        powder_pattern.sort_two_theta();
        powder_pattern.average_if_two_theta_equal();
        std::cout << powder_pattern.average_two_theta_step() << std::endl;
        powder_pattern.recalculate_estimated_standard_deviations();
        powder_pattern.save_xye( append_to_file_name( input_file_name, "_recal" ), false );
    MACRO_END_GAME

    try // Farey.
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

    try // Materials Studio cell.xcd file.
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

    try // Calculate molecular volume and packing coefficient.
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
        powder_pattern.save_xye( FileName( "powder_pattern.xye" ), false );
    MACRO_END_GAME

    try // Calculate < Pl( cos( theta ) ) >.
    {
        // S6
        {
        std::vector< Vector3D > points;
        points.push_back( Vector3D( 1.0, 0.0, 0.0 ) );
        for ( size_t i( 1 ); i != 6; ++i )
            points.push_back( rotate_point_about_axis( points[0], Vector3D(), NormalisedVector3D( 0.0, 0.0, 1.0 ), Angle::from_degrees( 60.0 ) ) );
        std::cout << "S6 = " << orientational_order_parameter_S6( points ) << std::endl;
        }
        size_t max_l = 8;
        ++max_l;
        // Square planar.
        {
        std::cout << "Square planar" << std::endl;
        std::vector< Angle > difference_angles;
        difference_angles.push_back( Angle::from_degrees(   0.0 ) );
        difference_angles.push_back( Angle::from_degrees(   0.0 ) );
        difference_angles.push_back( Angle::from_degrees(   0.0 ) );
        difference_angles.push_back( Angle::from_degrees(   0.0 ) );
        difference_angles.push_back( Angle::from_degrees(  90.0 ) );
        difference_angles.push_back( Angle::from_degrees(  90.0 ) );
        difference_angles.push_back( Angle::from_degrees(  90.0 ) );
        difference_angles.push_back( Angle::from_degrees(  90.0 ) );
        difference_angles.push_back( Angle::from_degrees(  90.0 ) );
        difference_angles.push_back( Angle::from_degrees(  90.0 ) );
        difference_angles.push_back( Angle::from_degrees(  90.0 ) );
        difference_angles.push_back( Angle::from_degrees(  90.0 ) );
        difference_angles.push_back( Angle::from_degrees( 180.0 ) );
        difference_angles.push_back( Angle::from_degrees( 180.0 ) );
        difference_angles.push_back( Angle::from_degrees( 180.0 ) );
        difference_angles.push_back( Angle::from_degrees( 180.0 ) );
        for ( size_t l( 0 ); l != max_l; ++l )
        {
            double average = 0.0;
            for ( size_t i( 0 ); i != difference_angles.size(); ++i )
                average += Legendre_polynomial( l, difference_angles[ i ].cosine() );
            average /= difference_angles.size();
            if ( nearly_zero( average ) )
                average = 0.0;
            std::cout << "l = " << l << ", average = " << average << std::endl;
        }
        std::cout << std::endl;
        std::vector< NormalisedVector3D > orientation_vectors;
        orientation_vectors.push_back( NormalisedVector3D( -1.0, -1.0, 0.0 ) );
        orientation_vectors.push_back( NormalisedVector3D( -1.0,  1.0, 0.0 ) );
        orientation_vectors.push_back( NormalisedVector3D(  1.0, -1.0, 0.0 ) );
        orientation_vectors.push_back( NormalisedVector3D(  1.0,  1.0, 0.0 ) );
        for ( size_t l( 0 ); l != max_l; ++l )
            std::cout << "l = " << l << ", S = " << orientational_order_parameter( l, orientation_vectors ) << std::endl;
        }
        // Tetragonal.
        {
        std::cout << "Tetragonal" << std::endl;
        std::vector< Angle > difference_angles;
        difference_angles.push_back( Angle::from_degrees(   0.0 ) );
        difference_angles.push_back( Angle::from_degrees(   0.0 ) );
        difference_angles.push_back( Angle::from_degrees(   0.0 ) );
        difference_angles.push_back( Angle::from_degrees(   0.0 ) );
        difference_angles.push_back( Angle::from_degrees( 109.47 ) );
        difference_angles.push_back( Angle::from_degrees( 109.47 ) );
        difference_angles.push_back( Angle::from_degrees( 109.47 ) );
        difference_angles.push_back( Angle::from_degrees( 109.47 ) );
        difference_angles.push_back( Angle::from_degrees( 109.47 ) );
        difference_angles.push_back( Angle::from_degrees( 109.47 ) );
        difference_angles.push_back( Angle::from_degrees( 109.47 ) );
        difference_angles.push_back( Angle::from_degrees( 109.47 ) );
        difference_angles.push_back( Angle::from_degrees( 109.47 ) );
        difference_angles.push_back( Angle::from_degrees( 109.47 ) );
        difference_angles.push_back( Angle::from_degrees( 109.47 ) );
        difference_angles.push_back( Angle::from_degrees( 109.47 ) );
        for ( size_t l( 0 ); l != max_l; ++l )
        {
            double average = 0.0;
            for ( size_t i( 0 ); i != difference_angles.size(); ++i )
                average += Legendre_polynomial( l, difference_angles[ i ].cosine() );
            average /= difference_angles.size();
            if ( nearly_zero( average, 0.0001 ) )
                average = 0.0;
            std::cout << "l = " << l << ", average = " << average << std::endl;
        }
        std::cout << std::endl;
        std::vector< NormalisedVector3D > orientation_vectors;
        orientation_vectors.push_back( NormalisedVector3D( -1.0, -1.0, -1.0 ) );
        orientation_vectors.push_back( NormalisedVector3D(  1.0,  1.0, -1.0 ) );
        orientation_vectors.push_back( NormalisedVector3D( -1.0,  1.0,  1.0 ) );
        orientation_vectors.push_back( NormalisedVector3D(  1.0, -1.0,  1.0 ) );
        for ( size_t l( 0 ); l != max_l; ++l )
            std::cout << "l = " << l << ", S = " << orientational_order_parameter( l, orientation_vectors ) << std::endl;
        }
        
//Square planar
//order = 0, average = 1
//order = 1, average = 0
//order = 2, average = 0.25
//order = 3, average = 0
//order = 4, average = 0.6875
//order = 5, average = 0
//order = 6, average = 0.34375
//
//Tetragonal
//order = 0, average = 1
//order = 1, average = 0
//order = 2, average = 0
//order = 3, average = 0.555546
//order = 4, average = 0.259287
//order = 5, average = 0
//order = 6, average = 0.395034
    MACRO_END_GAME

    try // Expand all CF3 groups into disordered groups.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        // Create a list of all F atoms.
        // For each F atom, find the parent C atom.
        Element F( "F" );
        Element C( "C" );
        std::vector< size_t > F_atoms;
        std::vector< size_t > C_atoms;
        for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
        {
            if ( crystal_structure.atom( i ).element() != F )
                continue;
            F_atoms.push_back( i );
            double shortest_distance( 99999.0 );
            double C_index;
            for ( size_t j( 0 ); j != crystal_structure.natoms(); ++j )
            {
                if ( crystal_structure.atom( j ).element() != C )
                    continue;
                double distance;
                Vector3D dummy_difference_vector;
                crystal_structure.shortest_distance( crystal_structure.atom( i ).position(), crystal_structure.atom( j ).position(), distance, dummy_difference_vector );
                if ( distance < shortest_distance )
                {
                    shortest_distance = distance;
                    C_index = j;
                }
            }
            if ( shortest_distance < 1.2 )
                std::cout << "Warning: shortest C-F distance is shorter than 1.2 A." << std::endl;
            if ( shortest_distance > 1.5 )
                std::cout << "Warning: shortest C-F distance is longer than 1.5 A." << std::endl;
            C_atoms.push_back( C_index );
        }
        // If there are not at least three, nothing to do.
        if ( F_atoms.size() < 3 )
            return 0;
        // For each C atom, count how many F atoms are bonded to it.
        std::vector< size_t > nneighbours( C_atoms.size(), 0 );
        std::vector< bool > done( C_atoms.size(), false );
        for ( size_t i( 0 ); i != C_atoms.size(); ++i )
        {
            if ( done[ i ] )
                continue;
            size_t noccurrences( 1 );
            done[ i ] = true;
            for ( size_t j( i+1 ); j != C_atoms.size(); ++j )
            {
                if ( done[ j ] )
                    continue;
                if ( C_atoms[ i ] == C_atoms[ j ] )
                    ++noccurrences;
            }
            nneighbours[ i ] = noccurrences;
        }
        size_t iCF3( 0 );
        for ( size_t i( 0 ); i != C_atoms.size(); ++i )
        {
            if ( nneighbours[ i ] != 3 )
                continue;
            ++iCF3;
            size_t main_C( C_atoms[ i ] );
            double shortest_distance( 99999.0 );
            size_t second_C;
            for ( size_t j( 0 ); j != crystal_structure.natoms(); ++j )
            {
                if ( main_C == j )
                    continue;
                if ( crystal_structure.atom( j ).element() != C )
                    continue;
                double distance;
                Vector3D dummy_difference_vector;
                crystal_structure.shortest_distance( crystal_structure.atom( main_C ).position(), crystal_structure.atom( j ).position(), distance, dummy_difference_vector );
                if ( distance < shortest_distance )
                {
                    shortest_distance = distance;
                    second_C = j;
                }
            }
            if ( shortest_distance < 1.2 )
                std::cout << "Warning: shortest C-C distance is shorter than 1.2 A." << std::endl;
            if ( shortest_distance > 1.6 )
                std::cout << "Warning: shortest C-C distance is longer than 1.6 A." << std::endl;
            Vector3D C1 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( main_C ).position() );
            Vector3D C2 = crystal_structure.crystal_lattice().fractional_to_orthogonal( crystal_structure.atom( second_C ).position() );
            NormalisedVector3D n = normalised_vector( C1 - C2 );
            std::vector< size_t > F_atoms_2;
            std::vector< std::string > labels;
            for ( size_t j( 0 ); j != F_atoms.size(); ++j )
            {
                if ( C_atoms[ j ] == main_C )
                {
                    F_atoms_2.push_back( F_atoms[ j ] );
                    labels.push_back( crystal_structure.atom( F_atoms[ j ] ).label() );
                }
            }
            std::string main_C_label = crystal_structure.atom( main_C ).label();
            std::cout << "    prm Focc" + size_t2string( iCF3 )+ " 0.5" << std::endl;
            for ( size_t j( 0 ); j != F_atoms_2.size(); ++j )
            {
                size_t iAtom = F_atoms_2[ j ];
                Vector3D position = crystal_structure.atom( iAtom ).position();
                std::cout << "    site " + labels[ j ] + "a x ref_flag" + " " + double2string_pad_plus( position.x(), 5, ' ' ) +
                                                          " y ref_flag" + " " + double2string_pad_plus( position.y(), 5, ' ' ) +
                                                          " z ref_flag" + " " + double2string_pad_plus( position.z(), 5, ' ' ) +
                             " occ F  = Focc" + size_t2string( iCF3 )+ "; beq = bnonh;" << std::endl;
            }
            for ( size_t j( 0 ); j != F_atoms_2.size(); ++j )
            {
                size_t iAtom = F_atoms_2[ j ];
                Vector3D position = crystal_structure.atom( iAtom ).position();
                position = crystal_structure.crystal_lattice().fractional_to_orthogonal( position );
                position = rotate_point_about_axis( position, C1, n, Angle::from_degrees( 60.0 ) );
                position = crystal_structure.crystal_lattice().orthogonal_to_fractional( position );
                Atom new_atom( crystal_structure.atom( iAtom ) );
                new_atom.set_occupancy( 0.5 );
                std::string label = new_atom.label();
                new_atom.set_label( label + "a" );
                crystal_structure.set_atom( iAtom, new_atom );
                new_atom.set_label( label + "b" );
                new_atom.set_position( position );
                crystal_structure.add_atom( new_atom );
                std::cout << "    site " + labels[ j ] + "b x ref_flag" + " " + double2string_pad_plus( position.x(), 5, ' ' ) +
                                                          " y ref_flag" + " " + double2string_pad_plus( position.y(), 5, ' ' ) +
                                                          " z ref_flag" + " " + double2string_pad_plus( position.z(), 5, ' ' ) +
                             " occ F  = 1.0 - Focc" + size_t2string( iCF3 ) + "; beq = bnonh;" << std::endl;
            }
            for ( size_t j( 0 ); j != F_atoms_2.size(); ++j )
                std::cout << "    Distance_Restrain( " + main_C_label + " " + labels[ j ] + "a, 1.38, 0.0, bond_width, bond_weight )" << std::endl;
            for ( size_t j( 0 ); j != F_atoms_2.size(); ++j )
                std::cout << "    Distance_Restrain( " + main_C_label + " " + labels[ j ] + "b, 1.38, 0.0, bond_width, bond_weight )" << std::endl;
            for ( size_t j( 0 ); j != F_atoms_2.size(); ++j )
                std::cout << "    Angle_Restrain( " + crystal_structure.atom( second_C ).label() + " " + main_C_label + " " + labels[ j ] + "a, 112.0, 0.0, angle_width, angle_weight )" << std::endl;
            std::cout << "    Angle_Restrain( " + labels[ 0 ] + "a " + main_C_label + " " + labels[ 1 ] + "a, 106.0, 0.0, angle_width, angle_weight )" << std::endl;
            std::cout << "    Angle_Restrain( " + labels[ 0 ] + "a " + main_C_label + " " + labels[ 2 ] + "a, 106.0, 0.0, angle_width, angle_weight )" << std::endl;
            std::cout << "    Angle_Restrain( " + labels[ 1 ] + "a " + main_C_label + " " + labels[ 2 ] + "a, 106.0, 0.0, angle_width, angle_weight )" << std::endl;
            for ( size_t j( 0 ); j != F_atoms_2.size(); ++j )
                std::cout << "    Angle_Restrain( " + crystal_structure.atom( second_C ).label() + " " + main_C_label + " " + labels[ j ] + "b, 112.0, 0.0, angle_width, angle_weight )" << std::endl;
            std::cout << "    Angle_Restrain( " + labels[ 0 ] + "b " + main_C_label + " " + labels[ 1 ] + "b, 106.0, 0.0, angle_width, angle_weight )" << std::endl;
            std::cout << "    Angle_Restrain( " + labels[ 0 ] + "b " + main_C_label + " " + labels[ 2 ] + "b, 106.0, 0.0, angle_width, angle_weight )" << std::endl;
            std::cout << "    Angle_Restrain( " + labels[ 1 ] + "b " + main_C_label + " " + labels[ 2 ] + "b, 106.0, 0.0, angle_width, angle_weight )" << std::endl;

        }
        crystal_structure.save_cif( append_to_file_name( input_file_name, "_CF3_split" ) );
    MACRO_END_GAME

    try // Adjust ADPs to site symmetry.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        SymmetryOperator symmetry_operator( "x,-y,z" );
        std::vector< Matrix3D > rotation_matrices;
        rotation_matrices.push_back( Matrix3D() );
        rotation_matrices.push_back( Matrix3D(  1,  0,  0,
                                                0, -1,  0,
                                                0,  0,  1 ) );
        PointGroup point_group( rotation_matrices );
        for ( size_t i( 0 ); i != crystal_structure.natoms(); ++i )
        {
            Atom atom( crystal_structure.atom( i ) );
            if ( atom.ADPs_type() != Atom::ANISOTROPIC )
                continue;
            if ( ! nearly_equal( atom.position().y(), 0.5 ) )
                continue;
            AnisotropicDisplacementParameters adps = atom.anisotropic_displacement_parameters();
            adps = adjust_to_site_symmetry( adps, point_group, crystal_structure.crystal_lattice() );
            atom.set_anisotropic_displacement_parameters( adps );
            crystal_structure.set_atom( i, atom );
        }
        crystal_structure.save_cif( append_to_file_name( input_file_name, "_ssa" ) ); // Site-symmetry adjusted
    MACRO_END_GAME

    try // Average atoms. The original cif must have atoms labels with "A" and "B", e.g. "C7A" and C7B".
        // The atoms to be averaged are related by a symmetry operator.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        CrystalStructure crystal_structure_new;
        crystal_structure_new.set_space_group( crystal_structure.space_group() );
        crystal_structure_new.set_crystal_lattice( crystal_structure.crystal_lattice() );
        crystal_structure_new.set_name( crystal_structure.name() );
//        crystal_structure.reserve_natoms( crystal_structure_org.natoms() );
        SymmetryOperator symmetry_operator( "x,-y,z" );
        for ( size_t iAtom( 0 ); iAtom != crystal_structure.natoms(); ++iAtom )
        {
            Atom atom( crystal_structure.atom( iAtom ) );
            std::string last_character = atom.label();
            if ( last_character.length() == 0 )
                throw std::runtime_error( "Atom has no label." );
            last_character = last_character.substr( last_character.length()-1, 1 );
            if ( last_character == "B" )
                continue;
            if ( last_character == "A" )
            {
                for ( size_t bAtom( 0 ); bAtom != crystal_structure.natoms(); ++bAtom )
                {
                    std::string last_character_b = crystal_structure.atom( bAtom ).label();
                    if ( last_character_b.length() == 0 )
                        throw std::runtime_error( "Atom has no label." );
                    last_character_b = last_character_b.substr( last_character_b.length()-1, 1 );
                    if ( last_character_b == "B" )
                    {
                        Atom atom_b( crystal_structure.atom( bAtom ) );
                        if ( atom.label().substr( 0, atom.label().length()-1 ) != atom_b.label().substr( 0, atom_b.label().length()-1 ) )
                            continue;
                        atom_b.set_position( symmetry_operator * atom_b.position() );
                        if ( atom_b.ADPs_type() == Atom::ANISOTROPIC )
                            atom_b.set_anisotropic_displacement_parameters( transform_adps( atom_b.anisotropic_displacement_parameters(), symmetry_operator.rotation(), crystal_structure.crystal_lattice() ) );
                        atom = average( atom, atom_b );
                        atom.set_label( atom.label().substr( 0, atom.label().length()-1 ) );
                        break;
                    }
                }
            }
            crystal_structure_new.add_atom( atom );
        }
        crystal_structure_new.save_cif( append_to_file_name( input_file_name, "_avg" ) );
    MACRO_END_GAME

    try // For unit cell ONLY, find structure in FileList.txt using Rene de Gelder's normalised weighted cross correlations.
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
 //       target_crystal_structure.apply_space_group_symmetry();
        Angle two_theta_start( 3.0, Angle::DEGREES );
        Angle two_theta_end(  35.0, Angle::DEGREES );
        Angle two_theta_step( 0.01, Angle::DEGREES );
        double FWHM( 0.1 );
        PowderPattern target_powder_pattern;
        {
            PowderPatternCalculator powder_pattern_calculator( target_crystal_structure );
            powder_pattern_calculator.set_two_theta_start( two_theta_start );
            powder_pattern_calculator.set_two_theta_end( two_theta_end );
            powder_pattern_calculator.set_two_theta_step( two_theta_step );
            powder_pattern_calculator.set_FWHM( FWHM );
            powder_pattern_calculator.calculate_reflection_list(); // F^2 is set to 1.0 by default.
            ReflectionList reflection_list = powder_pattern_calculator.reflection_list();
            powder_pattern_calculator.calculate( reflection_list, target_powder_pattern );
        }
        // Report best match and all matches over 0.95 (sorted).
        std::vector< double > all_matches_FoMs;
        std::vector< size_t > all_matches_indices;
        FileList file_list( file_list_file_name );
        if ( file_list.empty() )
            throw std::runtime_error( std::string( "No files in file list " ) + file_list_file_name.full_name() );
        double highest_correlation( 0.0 );
        size_t highest_correlation_index( 0 );
        for ( size_t i( 0 ); i != file_list.size(); ++i )
        {
            CrystalStructure crystal_structure;
            std::cout << "Now reading cif... " + file_list.value( i ).full_name() << std::endl;
            read_cif( file_list.value( i ), crystal_structure );
        //    crystal_structure.apply_space_group_symmetry();
//            std::cout << "Now calculating powder pattern... " + size_t2string( i, 4, '0' ) << std::endl;
            PowderPatternCalculator powder_pattern_calculator( crystal_structure );
            powder_pattern_calculator.set_two_theta_start( two_theta_start );
            powder_pattern_calculator.set_two_theta_end( two_theta_end );
            powder_pattern_calculator.set_two_theta_step( two_theta_step );
            powder_pattern_calculator.set_FWHM( FWHM );
            PowderPattern powder_pattern;
            powder_pattern_calculator.calculate_reflection_list(); // F^2 is set to 1.0 by default.
            ReflectionList reflection_list = powder_pattern_calculator.reflection_list();
            powder_pattern_calculator.calculate( reflection_list, powder_pattern );
            double correlation = normalised_weighted_cross_correlation( target_powder_pattern, powder_pattern, Angle( 1.5, Angle::DEGREES ) );
//            if ( correlation > 0.95 )
//                text_file_writer.write_line( double2string( correlation ) + " " + size_t2string( i+1 ) );
//            if ( correlation > highest_correlation )
//            {
//                highest_correlation = correlation;
//                highest_correlation_index = i;
//                std::cout << "highest_correlation = " << highest_correlation << std::endl;
//                std::cout << "highest_correlation_index = " << highest_correlation_index+1 << std::endl;
//            }
            if ( correlation > 0.90 )
            {
                all_matches_FoMs.push_back( correlation );
                all_matches_indices.push_back( i );
            }
        }
        Mapping mapping = sort( all_matches_FoMs );
        std::cout << "There were " << all_matches_FoMs.size() << " matches" << std::endl;
        for ( size_t i( 0 ); i != all_matches_FoMs.size(); ++i )
        {
            std::cout << "correlation = " << all_matches_FoMs[ mapping[ i ] ] << std::endl;
            std::cout << "correlation_index = " << all_matches_indices[ mapping[ i ] ]+1 << std::endl;
            std::cout << "File name = " << file_list.value( all_matches_indices[ mapping[ i ] ] ).full_name() << std::endl;
        }
//        std::cout << "highest_correlation = " << highest_correlation << std::endl;
//        std::cout << "highest_correlation_index = " << highest_correlation_index+1 << std::endl;
    MACRO_END_GAME

    try // Crystal structures database.
    {
        CrystalStructure crystal_structure_NaCl = NaCl();
        crystal_structure_NaCl.save_cif( FileName( "C:\\Users\\jacco\\Documents\\NaCl_C.cif" ) );
    MACRO_END_GAME

    try // Normalise C-F distances.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        normalise_C_F_bonds( crystal_structure );
        crystal_structure.save_cif( append_to_file_name( input_file_name, "_CFnorm" ) );
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
            powder_pattern_calculator.set_two_theta_start( two_theta_start );
            powder_pattern_calculator.set_two_theta_end( two_theta_end );
            powder_pattern_calculator.set_two_theta_step( two_theta_step );
            powder_pattern_calculator.set_FWHM( FWHM );
            powder_pattern_calculator.calculate( target_powder_pattern );
        }
        // Report best match and all matches over 0.95 (sorted).
        std::vector< double > all_matches_FoMs;
        std::vector< size_t > all_matches_indices;
        FileList file_list( file_list_file_name );
        if ( file_list.empty() )
            throw std::runtime_error( std::string( "No files in file list " ) + file_list_file_name.full_name() );
        double highest_correlation( 0.0 );
        size_t highest_correlation_index( 0 );
    //    TextFileWriter text_file_writer( FileName( "C:\\Data_Win\\matches.txt" ) );
        std::vector< std::string > water_labels;
//        water_labels.push_back( "O0_1" );
//        water_labels.push_back( "H0_1" );
//        water_labels.push_back( "H1_1" );
//        water_labels.push_back( "O0_2" );
//        water_labels.push_back( "H0_2" );
//        water_labels.push_back( "H1_2" );
        for ( size_t i( 0 ); i != file_list.size(); ++i )
        {
            CrystalStructure crystal_structure;
            std::cout << "Now reading cif... " + file_list.value( i ).full_name() << std::endl;
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
            powder_pattern_calculator.set_two_theta_start( two_theta_start );
            powder_pattern_calculator.set_two_theta_end( two_theta_end );
            powder_pattern_calculator.set_two_theta_step( two_theta_step );
            powder_pattern_calculator.set_FWHM( FWHM );
            PowderPattern powder_pattern;
            powder_pattern_calculator.calculate( powder_pattern );
            double correlation = normalised_weighted_cross_correlation( target_powder_pattern, powder_pattern, Angle( 3.0, Angle::DEGREES ) );
//            if ( correlation > 0.95 )
//                text_file_writer.write_line( double2string( correlation ) + " " + size_t2string( i+1 ) );
            if ( correlation > highest_correlation )
            {
                highest_correlation = correlation;
                highest_correlation_index = i;
                std::cout << "highest_correlation = " << highest_correlation << std::endl;
                std::cout << "highest_correlation_index = " << highest_correlation_index+1 << std::endl;
            }
            if ( correlation > 0.95 )
            {
                all_matches_FoMs.push_back( correlation );
                all_matches_indices.push_back( i );
            }
        }
        Mapping mapping = sort( all_matches_FoMs );
        std::cout << "There were " << all_matches_FoMs.size() << " matches" << std::endl;
        for ( size_t i( 0 ); i != all_matches_FoMs.size(); ++i )
        {
            std::cout << "correlation = " << all_matches_FoMs[ mapping[ i ] ] << std::endl;
            std::cout << "correlation_index = " << all_matches_indices[ mapping[ i ] ]+1 << std::endl;
        }
        std::cout << "highest_correlation = " << highest_correlation << std::endl;
        std::cout << "highest_correlation_index = " << highest_correlation_index+1 << std::endl;
    MACRO_END_GAME

    try // Average two crystal structures weighted by their energy (for creating disorder models for fixed cell optimisations).
    {
        if ( argc != 6 )
            throw std::runtime_error( "Please give the names of two .cif files, two energies and the number of atoms." );
        FileName input_file_name_1( argv[ 1 ] );
        CrystalStructure crystal_structure_1;
        read_cif_or_cell( input_file_name_1, crystal_structure_1 );
        FileName input_file_name_2( argv[ 2 ] );
        CrystalStructure crystal_structure_2;
        read_cif_or_cell( input_file_name_2, crystal_structure_2 );
        double energy_1 = string2double( argv[ 3 ] );
        double energy_2 = string2double( argv[ 4 ] );
        size_t natoms = string2integer( argv[ 5 ] );
        double energy;
        if ( energy_1 < energy_2 )
            energy = ( energy_2 - energy_1 ) * natoms;
        else
            energy = ( energy_1 - energy_2 ) * natoms;
        double T = 300.0; // K
        double R = 8.314462; // J/mol/K
        double RT = ( ( R * T ) / 1000.0 ) / 4.184; // kcal/mol/K
        double Boltzmann_weight = exp( -energy/RT );
        std::cout << "Boltzmann weight = " << Boltzmann_weight << std::endl;
        CrystalLattice crystal_lattice;
        if ( energy_1 < energy_2 )
            crystal_lattice = average( crystal_structure_1.crystal_lattice(), crystal_structure_2.crystal_lattice(), Boltzmann_weight );
        else
            crystal_lattice = average( crystal_structure_2.crystal_lattice(), crystal_structure_1.crystal_lattice(), Boltzmann_weight );
        std::cout << "Occupancies are " + double2string( 1.0/(1.0+Boltzmann_weight), 3 ) + " and " + double2string( Boltzmann_weight/(1.0+Boltzmann_weight), 3 ) << std::endl;
        std::cout << "_cell_length_a    " + double2string( crystal_lattice.a(), 5 ) << std::endl;
        std::cout << "_cell_length_b    " + double2string( crystal_lattice.b(), 5 ) << std::endl;
        std::cout << "_cell_length_c    " + double2string( crystal_lattice.c(), 5 ) << std::endl;
        std::cout << "_cell_angle_alpha " + double2string( crystal_lattice.alpha().value_in_degrees(), 5 ) << std::endl;
        std::cout << "_cell_angle_beta  " + double2string( crystal_lattice.beta().value_in_degrees(), 5  ) << std::endl;
        std::cout << "_cell_angle_gamma " + double2string( crystal_lattice.gamma().value_in_degrees(), 5 ) << std::endl;
        std::cout << "_cell_volume      " + double2string( crystal_lattice.volume(), 5 ) << std::endl;
        crystal_structure_1.set_crystal_lattice( crystal_lattice );
        crystal_structure_2.set_crystal_lattice( crystal_lattice );
        crystal_structure_1.save_cif( append_to_file_name( input_file_name_1, "_avguc" ) );
        crystal_structure_2.save_cif( append_to_file_name( input_file_name_2, "_avguc" ) );
    MACRO_END_GAME

    try // Rescale unit-cell volume.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        CrystalLattice crystal_lattice = crystal_structure.crystal_lattice();
        crystal_lattice.rescale_volume( 1105.76, 0 );
        crystal_structure.set_crystal_lattice( crystal_lattice );
        crystal_structure.save_cif( append_to_file_name( input_file_name, "_rv" ) );
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

    try // Loop over all space groups to test constraints.
    {
        TextFileReader_2 input_file( FileName( "IT.cif" ) );
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
            do // Read the symmetry operators.
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

    try // Write FileList.txt file.
    {
        TextFileWriter text_file_writer( FileName( "C:\\Data_Win\\\\FileList.txt" ) );
        for ( size_t i( 0 ); i != 7252; ++i )
        {
            text_file_writer.write_line( "structure_" + size_t2string( i+1, 6, '0' ) + ".cif" );
        }
    MACRO_END_GAME

    try // Add methyl hydrogen atoms.
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

    try // Loop over all space groups.
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
        TextFileReader_2 input_file( FileName( "IT.cif" ) );
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
            do // Read the symmetry operators.
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
            PointGroup Laue_class = space_group.Laue_class();
            for ( size_t iPOD( 0 ); iPOD != PO_directions.size(); ++iPOD )
            {
                Vector3D PO_vector = reciprocal_lattice_point( PO_directions[iPOD], crystal_lattice );
                for ( size_t iMD_par( 0 ); iMD_par != r_values.size(); ++iMD_par )
                {
                    double r = r_values[iMD_par];
                    for ( size_t iReflection( 0 ); iReflection != reflection_list.size(); ++iReflection )
                    {
                        // Calculate equivalent reflections and check that the March-Dollase PO corrections are the same for all of them.
                        Vector3D H = reciprocal_lattice_point( reflection_list.miller_indices( iReflection ), crystal_lattice );
                        Angle alpha = angle( PO_vector, H );
                        double reference_PO = std::pow( square(r) * square(alpha.cosine()) + square(alpha.sine())/r, -3.0/2.0 );
                        for ( size_t j( 0 ); j != Laue_class.nsymmetry_operators(); ++j )
                        {
                            MillerIndices equivalent_reflection = reflection_list.miller_indices( iReflection ) * Laue_class.symmetry_operator( j );
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

    try // Analyse KONTIQ results.
    {
        TextFileReader_2 RESULTS_ucfr_file( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE\\RESULTS_ucfr.txt" ) );
        std::string line;
        std::vector< std::string > words;
        std::vector< double > energies;
        for ( size_t i( 1 ); i != 7; ++i )
        {
            size_t iPos = RESULTS_ucfr_file.find_whole_word( "KONTIQ02_H" + size_t2string( i, 2, '0' ) + "_Hnorm" );
            line = RESULTS_ucfr_file.line( iPos );
            words = split( line );
            energies.push_back( string2double( words[2] ) );
        }
        double minimum = calculate_minimum( energies );
        CrystalStructure crystal_structure_1;
        CrystalStructure crystal_structure_2;
        double RMSCD( -1.0 );
        for ( size_t i( 1 ); i != 7; ++i )
        {
            std::cout << "KONTIQ02_H" + size_t2string( i, 2, '0' ) << std::endl;
            read_cif( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE", "KONTIQ02_H" + size_t2string( i, 2, '0' ) + "_Hnorm", "cif" ), crystal_structure_1 );
            read_cif( FileName( "C:\\Users\\jacco\\Documents\\For_Data_Archive\\Results_from_GRACE", "KONTIQ02_H" + size_t2string( i, 2, '0' ) + "_Hnorm_mi_ucfr", "cif" ), crystal_structure_2 );
            RMSCD = root_mean_square_Cartesian_displacement( crystal_structure_1, crystal_structure_2, false );
            std::cout << "RMSCD = " << RMSCD << std::endl;
            std::cout << "Energy = " << ( energies[i-1] - minimum ) * 21 << std::endl;
            std::cout << std::endl;
        }
        
//KONTIQ02_H01
//RMSCD = 0.111839
//Energy = 0.38425
//
//KONTIQ02_H02
//RMSCD = 0.110141
//Energy = 0.383187
//
//KONTIQ02_H03
//RMSCD = 0.110636
//Energy = 0.383882
//
//KONTIQ02_H04
//RMSCD = 1.18543
//Energy = 3.0338
//
//KONTIQ02_H05
//RMSCD = 0.0955155
//Energy = 0.0025704
//
//KONTIQ02_H06
//RMSCD = 0.0979055
//Energy = 0        
        
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
            text_file_writer.write_line( file_list.value( i ).file_name() + " " + double2string( density ) );
        }
    MACRO_END_GAME

    try // Recalculate variable slit intensities to fixed slit intensities.
    {
        MACRO_ONE_XYEFILENAME_AS_ARGUMENT
        std::cout << powder_pattern.average_two_theta_step() << std::endl;
        powder_pattern.convert_to_fixed_slit();
        powder_pattern.save_xye( append_to_file_name( input_file_name, "_fixed_slits" ), true );
    MACRO_END_GAME

    try // Fit exponential.
    {
        TextFileReader text_file_reader( FileName( "energies.txt" ) );
        std::vector< double > energies;
        std::string line;
        double lowest_energy;
        if ( text_file_reader.get_next_line( line ) )
        {
            lowest_energy = string2double( line );
            energies.push_back( 0.0 );
        }
        else
            throw std::runtime_error( "File is empty." );
        while ( text_file_reader.get_next_line( line ) )
        {
            energies.push_back( string2double( line ) - lowest_energy );
        }
        Histogram histogram( 0.0, energies[ energies.size()-1 ], 10 );
        histogram.add_data( energies );
        histogram.plot();
        histogram.show();
        // From the first and the last bin, take the middle of the bin and the natural logarithm of the frequency
        double x1 = histogram.middle_of_bin( 0 );
        double y1 = ln( histogram.bin( 0 ) );
        double x2 = histogram.middle_of_bin( 9 );
        double y2 = ln( histogram.bin( 9 ) );
        double b = ( y2 - y1 ) / ( x2 - x1 );
        double a = exp( 0.5 * ( y1 - b*x1 + y2 - b*x2 ) );
        std::cout << "End points only " << std::endl;
        std::cout << "a = " << a << std::endl;
        std::cout << "b = " << b << std::endl;
        std::cout << "y1 = " << histogram.bin( 0 ) << " = " << a*exp(b*x1) << std::endl;
        std::cout << "y2 = " << histogram.bin( 9 ) << " = " << a*exp(b*x2) << std::endl;

        // The energy value to extrapolate to.
        double x3 = 0.33915;;
        // Now integrate between x1 and x3
        double area2 = (a/b)*exp(b*x2);
        double area3 = (a/b)*exp(b*x3);
        std::cout << "area 3 = " << energies.size() * ( area3 / area2 ) << std::endl;
        std::cout << "difference = " << energies.size() * ( area3 / area2 ) - energies.size() << std::endl;

        std::vector< double > x;
        std::vector< double > y;
        // End points only
        x.push_back( histogram.middle_of_bin( 0 ) );
        y.push_back( histogram.bin( 0 ) );
        x.push_back( histogram.middle_of_bin( 9 ) );
        y.push_back( histogram.bin( 9 ) );
        fit_exponential( x, y, a, b );
        std::cout << "End points only " << std::endl;
        std::cout << "a = " << a << std::endl;
        std::cout << "b = " << b << std::endl;
        std::cout << "y1 = " << histogram.bin( 0 ) << " = " << a*exp(b*x1) << std::endl;
        std::cout << "y2 = " << histogram.bin( 9 ) << " = " << a*exp(b*x2) << std::endl;
        x.clear();
        y.clear();
        // We deliberately leave out the two end points, because they are less reliable
        for ( size_t i( 1 ); i != histogram.size()-1; ++i )
        {
            x.push_back( histogram.middle_of_bin( i ) );
            y.push_back( histogram.bin( i ) );
        }
        fit_exponential( x, y, a, b );
        area2 = (a/b)*exp(b*x2);
        area3 = (a/b)*exp(b*x3);
        std::cout << std::endl;
        std::cout << "Regression excluding end points" << std::endl;
        std::cout << "a = " << a << std::endl;
        std::cout << "b = " << b << std::endl;
        std::cout << "y1 = " << histogram.bin( 0 ) << " = " << a*exp(b*x1) << std::endl;
        std::cout << "y2 = " << histogram.bin( 9 ) << " = " << a*exp(b*x2) << std::endl;
        std::cout << "area 3 = " << energies.size() * ( area3 / area2 ) << std::endl;
        std::cout << "difference = " << energies.size() * ( area3 / area2 ) - energies.size() << std::endl;
        // Now with the end points
        x.push_back( histogram.middle_of_bin( 0 ) );
        y.push_back( histogram.bin( 0 ) );
        x.push_back( histogram.middle_of_bin( 9 ) );
        y.push_back( histogram.bin( 9 ) );
        fit_exponential( x, y, a, b );
        area2 = (a/b)*exp(b*x2);
        area3 = (a/b)*exp(b*x3);
        std::cout << std::endl;
        std::cout << "Regression including end points" << std::endl;
        std::cout << "a = " << a << std::endl;
        std::cout << "b = " << b << std::endl;
        std::cout << "y1 = " << histogram.bin( 0 ) << " = " << a*exp(b*x1) << std::endl;
        std::cout << "y2 = " << histogram.bin( 9 ) << " = " << a*exp(b*x2) << std::endl;
        std::cout << "area 3 = " << energies.size() * ( area3 / area2 ) << std::endl;
        std::cout << "difference = " << energies.size() * ( area3 / area2 ) - energies.size() << std::endl;
    MACRO_END_GAME

    try // Loop over all space groups to test centring vectors.
    {
        TextFileReader_2 input_file( FileName( "IT.cif" ) );
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
            do // Read the symmetry operators.
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
            {
                SpaceGroup space_group( symmetry_operators, space_group_name );
                std::string centring_from_space_group_name = space_group_name.substr( 0, 1 );
                std::string centring_from_symmetry_operators = space_group.centring().centring_name();
                if ( centring_from_space_group_name != centring_from_symmetry_operators )
                    std::cout << "Error" << std::endl;
            }
            iLine += 10;
        }
    MACRO_END_GAME

    try // Make atom labels unique.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        crystal_structure.make_atom_labels_unique();
        crystal_structure.save_cif( append_to_file_name( input_file_name, "_ual" ) );
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

    try // Calculate similarity matrix.
    {
        MACRO_ONE_FILELISTNAME_OR_LIST_OF_FILES_AS_ARGUMENT
        CorrelationMatrix similarity_matrix = calculate_correlation_matrix( file_list );
        similarity_matrix.save( FileName( "SimilarityMatrix.txt" ) );
    MACRO_END_GAME

    try // Calculate similarity matrix based on unit cells.
    {
        MACRO_ONE_FILELISTNAME_OR_LIST_OF_FILES_AS_ARGUMENT
        CorrelationMatrix similarity_matrix = calculate_correlation_matrix_1( file_list );
        similarity_matrix.save( FileName( "SimilarityMatrix_1.txt" ) );
    MACRO_END_GAME

    try // Average two unit cells and normalise X-H bonds.
    {
        if ( argc != 3 )
            throw std::runtime_error( "Please give the name of two .cif or .cell files." );
        FileName file_name_1( argv[ 1 ] );
        CrystalStructure crystal_structure_1;
        read_cif_or_cell( file_name_1, crystal_structure_1 );
        FileName file_name_2( argv[ 2 ] );
        CrystalStructure crystal_structure_2;
        read_cif_or_cell( file_name_2, crystal_structure_2 );
        CrystalLattice average_crystal_lattice = average( crystal_structure_1.crystal_lattice(), crystal_structure_2.crystal_lattice() );
        normalise_X_H_bonds( crystal_structure_1 );
        normalise_X_H_bonds( crystal_structure_2 );
        crystal_structure_1.save_cif( replace_extension( append_to_file_name( file_name_1, "_Hnorm" ) , "cif" ) );
        crystal_structure_2.save_cif( replace_extension( append_to_file_name( file_name_2, "_Hnorm" ) , "cif" ) );
        crystal_structure_1.set_crystal_lattice( average_crystal_lattice );
        crystal_structure_2.set_crystal_lattice( average_crystal_lattice );
        crystal_structure_1.save_cif( replace_extension( append_to_file_name( file_name_1, "_avguc_Hnorm" ) , "cif" ) );
        crystal_structure_2.save_cif( replace_extension( append_to_file_name( file_name_2, "_avguc_Hnorm" ) , "cif" ) );
    MACRO_END_GAME

    try // Average three unit cells and normalise X-H bonds.
    {
        if ( argc != 4 )
            throw std::runtime_error( "Please give the name of three .cif or .cell files." );
        FileName file_name_1( argv[ 1 ] );
        CrystalStructure crystal_structure_1;
        read_cif_or_cell( file_name_1, crystal_structure_1 );
        FileName file_name_2( argv[ 2 ] );
        CrystalStructure crystal_structure_2;
        read_cif_or_cell( file_name_2, crystal_structure_2 );
        FileName file_name_3( argv[ 3 ] );
        CrystalStructure crystal_structure_3;
        read_cif_or_cell( file_name_3, crystal_structure_3 );
        CrystalLattice average_crystal_lattice = average( crystal_structure_1.crystal_lattice(), crystal_structure_2.crystal_lattice() );
        average_crystal_lattice = average( crystal_structure_3.crystal_lattice(), average_crystal_lattice, 2.0 );
        normalise_X_H_bonds( crystal_structure_1 );
        normalise_X_H_bonds( crystal_structure_2 );
        normalise_X_H_bonds( crystal_structure_3 );
        crystal_structure_1.save_cif( replace_extension( append_to_file_name( file_name_1, "_Hnorm" ) , "cif" ) );
        crystal_structure_2.save_cif( replace_extension( append_to_file_name( file_name_2, "_Hnorm" ) , "cif" ) );
        crystal_structure_3.save_cif( replace_extension( append_to_file_name( file_name_3, "_Hnorm" ) , "cif" ) );
        crystal_structure_1.set_crystal_lattice( average_crystal_lattice );
        crystal_structure_2.set_crystal_lattice( average_crystal_lattice );
        crystal_structure_3.set_crystal_lattice( average_crystal_lattice );
        crystal_structure_1.save_cif( replace_extension( append_to_file_name( file_name_1, "_avguc_Hnorm" ) , "cif" ) );
        crystal_structure_2.save_cif( replace_extension( append_to_file_name( file_name_2, "_avguc_Hnorm" ) , "cif" ) );
        crystal_structure_3.save_cif( replace_extension( append_to_file_name( file_name_3, "_avguc_Hnorm" ) , "cif" ) );
    MACRO_END_GAME

    try // Check space group.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        check_if_closed( crystal_structure.space_group().symmetry_operators() );
    MACRO_END_GAME

    try // doubles_as_table().
    {
        std::vector< std::string > result;
        std::string padding( "                    " );
        std::vector< double > doubles;
        Element element( "D" );
        doubles.push_back( element.solid_state_volume() );
        for ( size_t i( 1 ); i != 113; ++i )
            doubles.push_back( Element( i ).solid_state_volume() );
        
        result = doubles_as_table( doubles,
                                   2,
                                   "> ",
                                   " | ",
                                   " <",
                                   100,
                                   false,
                                   0 );
        for ( size_t i( 0 ); i != result.size(); ++i )
            std::cout << result[i] << std::endl;
            std::cout << std::endl;

        result = doubles_as_table( doubles,
                                   2,
                                   ", ",
                                   100 - padding.length() - 4,
                                   false,
                                   0 );
        for ( size_t i( 0 ); i != result.size(); ++i )
            std::cout << padding + "    " + result[i] << std::endl;
            std::cout << std::endl;

        result = generate_Gauss_Legendre_quadrature_code( -1.0, 1.0, 30, 6, 100, true );
        for ( size_t i( 0 ); i != result.size(); ++i )
            std::cout << result[i] << std::endl;
    MACRO_END_GAME

    try // Remove H atoms from a set of cif files.
    {
        MACRO_ONE_FILELISTNAME_OR_LIST_OF_FILES_AS_ARGUMENT
        for ( size_t i( 0 ); i != file_list.size(); ++i )
            remove_hydrogen_atoms( file_list.value( i ) );
    MACRO_END_GAME

    try // Loop over all space groups to determine which symmetry operators commute.
    {
        TextFileReader_2 input_file( FileName( "IT.cif" ) );
        if ( input_file.size() != 8795 )
            throw std::runtime_error( "read_cif(): symmetry line must have same number of items as specified in loop." );
        size_t iLine( 0 );
        std::vector< std::string > words;
//        Matrix3D matrix( 19.0, 13.0, 37.0, 
//                          3.0, 11.0,  7.0,
//                        113.0,  4.0, 23.0  );
        Matrix3D matrix( 19.0, 13.0, 37.0, 
                         13.0, 11.0,  7.0,
                         37.0,  7.0, 23.0  );
        Vector3D vector( 11.0, 7.0, 13.0 );
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
            do // Read the symmetry operators.
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
            {
                bool matrix_failed( false );
                bool vector_failed( false );
                for ( size_t j( 0 ); j != symmetry_operators.size(); ++j )
                {
                    {
                    Matrix3D left  = matrix * symmetry_operators[j].rotation();
                    Matrix3D right = symmetry_operators[j].rotation() * matrix;
                    if ( ! nearly_equal( left, right ) )
                        matrix_failed = true;
                    }
                    {
                    Vector3D left  = vector * symmetry_operators[j].rotation();
                    Vector3D right = symmetry_operators[j].rotation() * vector;
                    if ( ! nearly_equal( left, right ) )
                        vector_failed = true;
                    }
                }
                if ( matrix_failed )
                    std::cout << space_group_name << " Matrix" << std::endl;
                if ( vector_failed )
                    std::cout << space_group_name << " Vector" << std::endl;
                
//                SpaceGroup space_group( symmetry_operators, space_group_name );
//                std::string centring_from_space_group_name = space_group_name.substr( 0, 1 );
//                std::string centring_from_symmetry_operators = space_group.centring();
//                if ( centring_from_space_group_name != centring_from_symmetry_operators )
//                    std::cout << "Error" << std::endl;
            }
            iLine += 10;
        }
    MACRO_END_GAME

    try // Show average 2theta step.
    {
        MACRO_ONE_XYEFILENAME_AS_ARGUMENT
        std::cout << "Average 2theta step: " << powder_pattern.average_two_theta_step() << std::endl;
    MACRO_END_GAME

    try // RMSCD X.cell and X-out.cif with matching without shifts.
    {
        if ( argc != 2 )
            throw std::runtime_error( "Please give the identifier of a pair of X.cell and X-out.cif files." );
        FileName file_name_1( argv[ 1 ] );
        CrystalStructure crystal_structure_1;
        read_cif_or_cell( file_name_1, crystal_structure_1 );
        if ( to_upper( file_name_1.extension() ) == "CELL" )
            crystal_structure_1.save_cif( replace_extension( file_name_1, "cif" ) );
        else
            crystal_structure_1.reduce_to_asymmetric_unit();
        FileName file_name_2( replace_extension( append_to_file_name( file_name_1, "-out" ) , "cif" ) );
        CrystalStructure crystal_structure_2;
        read_cif_or_cell( file_name_2, crystal_structure_2 );
        if ( to_upper( file_name_2.extension() ) == "CELL" )
            crystal_structure_2.save_cif( replace_extension( file_name_2, "cif" ) );
        else
            crystal_structure_1.reduce_to_asymmetric_unit();
        double result = RMSCD_with_matching( crystal_structure_1, crystal_structure_2, 0, false, false );
        std::cout << result << std::endl;
    MACRO_END_GAME

    try // Print final energy from .castep file.
    {
        if ( argc != 2 )
            throw std::runtime_error( "Please give the name of a .castep file." );
        FileName file_name( argv[ 1 ] );
        // .castep files are HUGE, must use a TextFileReader, not a TextFileReader_2.
        TextFileReader text_file_reader( file_name );
        if ( file_name.extension() != "castep" )
            std::cout << "Warning: file extension is not .castep" << std::endl;

// BFGS: finished iteration    50 with enthalpy= -1.51074945E+004 eV
//  
// +-----------+-----------------+-----------------+------------+-----+ <-- BFGS
// | Parameter |      value      |    tolerance    |    units   | OK? | <-- BFGS
// +-----------+-----------------+-----------------+------------+-----+ <-- BFGS
// |  dE/ion   |   7.394873E-006 |   1.077880E-005 |         eV | Yes | <-- BFGS
// |  |F|max   |   2.783228E-002 |   3.036700E-002 |       eV/A | Yes | <-- BFGS
// |  |dR|max  |   2.976385E-003 |   3.000000E-003 |          A | Yes | <-- BFGS
// +-----------+-----------------+-----------------+------------+-----+ <-- BFGS
//  
// BFGS: Geometry optimization completed successfully.
        std::string target_string( " BFGS: finished iteration" );
        std::string line_with_energy;
        std::string line;
        while ( text_file_reader.get_next_line( line ) )
        {
            if ( line.length() < target_string.length() )
                continue;
            if ( line.substr( 0, target_string.length() ) != target_string )
                continue;
            line_with_energy = line;
            if ( ! text_file_reader.get_next_line( line ) )
                throw std::runtime_error( "Expected line not found 1." );
            if ( line != "  " )
                std::cout << "Unexpected line contents 1." << std::endl;
            if ( ! text_file_reader.get_next_line( line ) )
                throw std::runtime_error( "Expected line not found 2." );
            if ( line != " +-----------+-----------------+-----------------+------------+-----+ <-- BFGS" )
                std::cout << "Unexpected line contents 2." << std::endl;
            if ( ! text_file_reader.get_next_line( line ) )
                throw std::runtime_error( "Expected line not found 3." );
            if ( line != " | Parameter |      value      |    tolerance    |    units   | OK? | <-- BFGS" )
                std::cout << "Unexpected line contents 3." << std::endl;
            if ( ! text_file_reader.get_next_line( line ) )
                throw std::runtime_error( "Expected line not found 4." );
            if ( line != " +-----------+-----------------+-----------------+------------+-----+ <-- BFGS" )
                std::cout << "Unexpected line contents 4." << std::endl;
            if ( ! text_file_reader.get_next_line( line ) )
                throw std::runtime_error( "Expected line not found 5." );
            if ( line.substr( 0, 14 ) != " |  dE/ion   |" )
                std::cout << "Unexpected line contents 5." << std::endl;
            if ( ! text_file_reader.get_next_line( line ) )
                throw std::runtime_error( "Expected line not found 6." );
            if ( line.substr( 0, 14 ) != " |  |F|max   |" )
                std::cout << "Unexpected line contents 6." << std::endl;
            if ( ! text_file_reader.get_next_line( line ) )
                throw std::runtime_error( "Expected line not found 7." );
            if ( line.substr( 0, 14 ) != " |  |dR|max  |" )
                std::cout << "Unexpected line contents 7." << std::endl;
            if ( ! text_file_reader.get_next_line( line ) )
                throw std::runtime_error( "Expected line not found 8." );
            if ( line != " +-----------+-----------------+-----------------+------------+-----+ <-- BFGS" )
                std::cout << "Unexpected line contents 8." << std::endl;
            if ( ! text_file_reader.get_next_line( line ) )
                throw std::runtime_error( "Expected line not found 9." );
            if ( line != "  " )
                std::cout << "Unexpected line contents 9." << std::endl;
            if ( ! text_file_reader.get_next_line( line ) )
                throw std::runtime_error( "Expected line not found 10." );
            if ( line != " BFGS: Geometry optimization completed successfully." )
                continue;
        }
        std::vector< std::string > words = split( line_with_energy );
        if ( words.size() != 8 )
            throw std::runtime_error( "Unexpected line format." );
        double result = string2double( words[6] ); // in eV
        std::cout << "Energy = " << double2string_2( result, 5 ) << " eV / unit cell" << std::endl;
        std::cout << "Energy = " << double2string_2( result*23.061, 5 ) << " kcal/mol unit cell" << std::endl;
    MACRO_END_GAME

    try // .cell to .cif .
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

    try // add_hydrogen_atom_to_sp3_atom().
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

    try // Add hydrogen atom to sp2 ring carbon.
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

    try // Correct carbon atom in aromatic ring.
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

    try // Add hydrogen atoms to sp3 nitrogen with hydrogen bond.
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
            text_file_writer.write_line( file_list_1.value( fi ).name() + " " + file_list_2.value( fi ).name() + " " + double2string( result ) );
        }
    MACRO_END_GAME

    try // Reorder atoms by molecule.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        CrystalStructure crystal_structure_2;
        crystal_structure_2.reserve_natoms( crystal_structure.natoms() );
        crystal_structure.perceive_molecules( true );
        crystal_structure_2.set_space_group( crystal_structure.space_group() );
        crystal_structure_2.set_crystal_lattice( crystal_structure.crystal_lattice() );
        crystal_structure_2.set_name( crystal_structure.name() );
        for ( size_t i( 0 ); i != crystal_structure.nmolecules(); ++i )
        {
            MoleculeInCrystal molecule = crystal_structure.molecule_in_crystal( i );
            for ( size_t j( 0 ); j != molecule.natoms(); ++j )
            {
                Atom atom( molecule.atom( j ) );
                atom.set_label( atom.element().symbol() + size_t2string( j+1 ) + "_"  + size_t2string( i+1 ) );
                crystal_structure_2.add_atom( atom );
            }
        }
        crystal_structure_2.save_cif( append_to_file_name( input_file_name, "_reordered" ) );
    MACRO_END_GAME

    try // RMSCD of a structure and its energy-minimised counterpart. Without matching.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        std::string file_name_str = input_file_name.name();
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
        std::cout << result << " A" << std::endl;
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

    try // find_match( CrystalStructure, CrystalStructure ) crystal structure 2 (rhs) is the one that gets changed, so 1 is the target.
    {
        if ( argc < 3 )
            throw std::runtime_error( "Please give the name of two .cif files." );
        FileName file_name_1( argv[ 1 ] );
        CrystalStructure crystal_structure_1;
        read_cif_or_cell( file_name_1, crystal_structure_1 );
        FileName file_name_2( argv[ 2 ] );
        CrystalStructure crystal_structure_2;
        read_cif_or_cell( file_name_2, crystal_structure_2 );
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
        crystal_structure_2.save_cif( replace_extension( append_to_file_name( file_name_2, "_mtrans" ), "cif" ) );
    MACRO_END_GAME

    try // Add OH hydrogen atom.
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

    try // Add centre of symmetry, not at origin.
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

    try // Find unit cell with angles close to 90.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        crystal_structure.transform( crystal_structure.crystal_lattice().choose_angles_close_to_90() );
        crystal_structure.crystal_lattice().print();
    MACRO_END_GAME

    try // Recalculate ESDs XRPD pattern.
    {
        MACRO_ONE_XYEFILENAME_AS_ARGUMENT
        powder_pattern.recalculate_estimated_standard_deviations();
        std::cout << powder_pattern.average_two_theta_step() << std::endl;
        powder_pattern.save_xye( append_to_file_name( input_file_name, "_new_ESDs" ), true );
    MACRO_END_GAME

    try // RMSCDs of a list of structures and their energy-minimised counterparts. Without matching.
    {
        MACRO_ONE_FILELISTNAME_AS_ARGUMENT
        std::string minimisation_ID( "_mi_TMFF" );
        TextFileWriter text_file_writer( FileName( file_list_file_name.directory(), "RMSCDs" + minimisation_ID, "txt" ) );
        for ( size_t i( 0 ); i != file_list.size(); ++i )
        {
            FileName input_file_name( file_list.value( i ) );
            CrystalStructure crystal_structure_1;
            read_cif( input_file_name, crystal_structure_1 );
            CrystalStructure crystal_structure_2;
            FileName file_name_2( append_to_file_name( input_file_name, minimisation_ID ) );
            read_cif( file_name_2, crystal_structure_2 );
            double result = root_mean_square_Cartesian_displacement( crystal_structure_1, crystal_structure_2, false );
            std::cout << "Structure 1 = " << input_file_name.full_name() << std::endl;
            std::cout << "Structure 2 = " << file_name_2.full_name() << std::endl;
            std::cout << result << " A" << std::endl;
            text_file_writer.write_line( input_file_name.file_name() + " " + double2string( result, 5 ) );
            std::cout << std::endl;
        }
    MACRO_END_GAME

    try // RMSCDs of a list of structures and their energy-minimised counterparts. Without matching.
    {
        MACRO_ONE_FILELISTNAME_AS_ARGUMENT
        {
        std::string minimisation_ID( "_mi_ucfr" );
        TextFileWriter text_file_writer( FileName( file_list_file_name.directory(), "RMSCDs" + minimisation_ID, "txt" ) );
        for ( size_t i( 0 ); i != file_list.size(); ++i )
        {
            FileName input_file_name( file_list.value( i ) );
            CrystalStructure crystal_structure_1;
            read_cif( input_file_name, crystal_structure_1 );
//            crystal_structure_1.space_group().show();
            CrystalStructure crystal_structure_2;
            FileName file_name_2( append_to_file_name( input_file_name, minimisation_ID ) );
            read_cif( file_name_2, crystal_structure_2 );
//            crystal_structure_2.space_group().show();
            double result = root_mean_square_Cartesian_displacement( crystal_structure_1, crystal_structure_2, false );
            std::cout << "Structure 1 = " << input_file_name.full_name() << std::endl;
            std::cout << "Structure 2 = " << file_name_2.full_name() << std::endl;
            std::cout << result << " A" << std::endl;
            text_file_writer.write_line( input_file_name.file_name() + " " + double2string( result, 5 ) );
            std::cout << std::endl;
        }
        }

        {
        std::string minimisation_ID( "_PBE0_mi_ucfr" );
        TextFileWriter text_file_writer( FileName( file_list_file_name.directory(), "RMSCDs" + minimisation_ID, "txt" ) );
        for ( size_t i( 0 ); i != file_list.size(); ++i )
        {
            FileName input_file_name( file_list.value( i ) );
            CrystalStructure crystal_structure_1;
            read_cif( input_file_name, crystal_structure_1 );
//            crystal_structure_1.space_group().show();
            CrystalStructure crystal_structure_2;
            FileName file_name_2( append_to_file_name( input_file_name, minimisation_ID ) );
            read_cif( file_name_2, crystal_structure_2 );
//            crystal_structure_2.space_group().show();
            double result = root_mean_square_Cartesian_displacement( crystal_structure_1, crystal_structure_2, false );
            std::cout << "Structure 1 = " << input_file_name.full_name() << std::endl;
            std::cout << "Structure 2 = " << file_name_2.full_name() << std::endl;
            std::cout << result << " A" << std::endl;
            text_file_writer.write_line( input_file_name.file_name() + " " + double2string( result, 5 ) );
            std::cout << std::endl;
        }
        }

        {
        std::string minimisation_ID( "_PBE0_MBD_mi_ucfr" );
        TextFileWriter text_file_writer( FileName( file_list_file_name.directory(), "RMSCDs" + minimisation_ID, "txt" ) );
        for ( size_t i( 0 ); i != file_list.size(); ++i )
        {
            FileName input_file_name( file_list.value( i ) );
            CrystalStructure crystal_structure_1;
            read_cif( input_file_name, crystal_structure_1 );
//            crystal_structure_1.space_group().show();
            CrystalStructure crystal_structure_2;
            FileName file_name_2( append_to_file_name( input_file_name, minimisation_ID ) );
            read_cif( file_name_2, crystal_structure_2 );
//            crystal_structure_2.space_group().show();
            double result = root_mean_square_Cartesian_displacement( crystal_structure_1, crystal_structure_2, false );
            std::cout << "Structure 1 = " << input_file_name.full_name() << std::endl;
            std::cout << "Structure 2 = " << file_name_2.full_name() << std::endl;
            std::cout << result << " A" << std::endl;
            text_file_writer.write_line( input_file_name.file_name() + " " + double2string( result, 5 ) );
            std::cout << std::endl;
        }
        }
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

    try // Print c.o.m.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        Vector3D com = crystal_structure.centre_of_mass( true );
        std::cout << "Centre of mass = " << std::endl;
        com.show();
    MACRO_END_GAME

    try // Split cif with multiple crystal structures.
    {
        FileName file_name( "GF.cif" );
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

    try // Add two hydrogen atoms to an sp3 atom.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT

        {
        std::string neighbour_1_atom_label( "C20" );
        std::string origin_atom_label( "C21" ); // This is the atom to which the hydrogen atoms are added
        std::string neighbour_2_atom_label( "N4" );
        Vector3D neighbour_1_atom_frac = crystal_structure.atom( crystal_structure.find_label( neighbour_1_atom_label ) ).position();
        Vector3D origin_atom_frac = crystal_structure.atom( crystal_structure.find_label( origin_atom_label ) ).position();
        Vector3D neighbour_2_atom_frac = crystal_structure.atom( crystal_structure.find_label( neighbour_2_atom_label ) ).position();
        // Find shortest distance between the two taking periodicity and space-group symmetry into account.
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
        // Find shortest distance between the two taking periodicity and space-group symmetry into account.
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

    try // Convert powder pattern in .avxrd format to .xye.
    {
        // The header lines must be removed so that only the 2theta values and the intensities are left
        if ( argc != 2 )
            throw std::runtime_error( "Please give the name of a .avxrd file." );
        FileName input_file_name( argv[ 1 ] );
        PowderPattern powder_pattern;
        TextFileReader text_file_reader( input_file_name );
        std::vector< std::string > words;
        size_t iCounter = 0;
        while ( text_file_reader.get_next_line( words ) )
        {
            std::cout << "iCounter = " << iCounter << std::endl;
            ++iCounter;
            if ( words.size() < 2 )
                throw std::runtime_error( "Cannot interpret line \"" + text_file_reader.get_line() + "\"" );

            std::cout << "0{{" << words[0] << "}}" << ", length = " << words[0].length() << std::endl;
            std::cout << "1{{" << words[1] << "}}" << ", length = " << words[1].length() << std::endl;
            std::cout << "2{{" << words[2] << "}}" << ", length = " << words[2].length() << std::endl;
            std::cout << std::endl;
            powder_pattern.push_back( Angle( string2double( words[0] ), Angle::DEGREES ), string2double( words[1] ) );
        }
        std::cout << powder_pattern.average_two_theta_step() << std::endl;
        powder_pattern.save_xye( append_to_file_name( input_file_name, "_converted" ), true );
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
                new_line += " " + words[1] + " " + extract_variable_value( text_file_reader.line( text_file_reader.find( words[2], iPos_Variables ) ), splitter );
            if ( words.size() >= 5 )
                new_line += " " + words[3] + " " + extract_variable_value( text_file_reader.line( text_file_reader.find( words[4], iPos_Variables ) ), splitter );
            if ( words.size() >= 7 )
                new_line += " " + words[5] + " " + extract_variable_value( text_file_reader.line( text_file_reader.find( words[6], iPos_Variables ) ), splitter );
            if ( words.size() >= 9 )
                throw std::runtime_error( "Cannot interpret line with >7 words." );
            text_file_writer.write_line( new_line );
        }
    MACRO_END_GAME

    try // P-1 to I-1.
    {
        CrystalStructure crystal_structure;
        SpaceGroup space_group;
        space_group.add_inversion_at_origin();
        crystal_structure.set_space_group( space_group );
        Matrix3D tranformation_matrix( -1.0, -1.0, -1.0,
                                        1.0,  0.0,  0.0,
                                        0.0, -1.0,  1.0 );
        crystal_structure.transform( tranformation_matrix );
        add_centring_to_space_group_after_transformation( tranformation_matrix, space_group );
        crystal_structure.set_space_group( space_group );
        crystal_structure.space_group().show();
    MACRO_END_GAME

    try // Check space group.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        check_if_closed( crystal_structure.space_group().symmetry_operators() );
    MACRO_END_GAME

    try // Normalise X-H distances.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        normalise_X_H_bonds( crystal_structure );
        crystal_structure.save_cif( append_to_file_name( input_file_name, "_Hnorm" ) );
    MACRO_END_GAME

    try // Subtract two powder patterns.
    {
        FileName input_file_name_1( "GP_profile.xye" );
        PowderPattern powder_pattern_1;
        powder_pattern_1.read_xye( input_file_name_1 );
        FileName input_file_name_2( "GP_BKGR.xye" );
        PowderPattern powder_pattern_2;
        powder_pattern_2.read_xye( input_file_name_2 );
        powder_pattern_1 -= powder_pattern_2;

        powder_pattern_1.correct_zero_point_error( Angle::from_degrees( 0.14068 ) );
        powder_pattern_1.save_xye( append_to_file_name( input_file_name_1, "_zp" ), true );
    MACRO_END_GAME

    try // Loop over all space groups to identify all translations in standard symmetry operators.
    {
        TextFileReader_2 input_file( FileName( "IT.cif" ) );
        if ( input_file.size() != 8795 )
            throw std::runtime_error( "Incorrect number of lines." );
        size_t iLine( 0 );
        std::vector< std::string > words;
        Tally< Fraction > counts;
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
            do // Read the symmetry operators.
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
            for ( size_t j( 0 ); j != space_group.nsymmetry_operators(); ++j )
            {
                Vector3D translation_vector = space_group.symmetry_operator( j ).translation();
                for ( size_t k( 0 ); k != 3; ++k )
                {
                    double translation = translation_vector.value( k );
                    Fraction fraction = double2fraction( translation, Fraction( 1, 12 ) );
                    counts.add( fraction );
                }
            }
            iLine += 10;
        }
        counts.show();
    MACRO_END_GAME

    try // Calculate all delta_AB for ADPs for all pairs of atoms to check if rigid-body approximation for TLS is valid.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        if ( crystal_structure.natoms() == 0 )
            return 0;
        std::vector< std::string > labels;
        std::vector< double > deltas;
        double maximum_delta( 0.0 );
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
                double delta = absolute( mu_2_v_i - mu_2_v_j );
                deltas.push_back( delta );
                if ( delta > maximum_delta )
                    maximum_delta = delta;
//                std::cout << crystal_structure.atom( i ).label() << " " << crystal_structure.atom( j ).label() << " " << delta << std::endl;
            }
        }
        for ( size_t i( 0 ); i != deltas.size(); ++i )
            std::cout << labels[i] << " " << ASCII_histogram( 0, maximum_delta, deltas[i], 25, ' ' ) << deltas[i] << std::endl;
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

    try // Average two atoms into one.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        {
        std::string atom_1_label( "C39A" );
        std::string atom_2_label( "C39B" );
        Atom atom_1 = crystal_structure.atom( crystal_structure.atom( atom_1_label ) );
        Atom atom_2 = crystal_structure.atom( crystal_structure.atom( atom_2_label ) );
        Vector3D atom_1_Cartesian = crystal_structure.crystal_lattice().fractional_to_orthogonal( atom_1.position() );
        Vector3D atom_2_Cartesian = crystal_structure.crystal_lattice().fractional_to_orthogonal( atom_2.position() );
        // Check that the two atoms are the same element.
        if ( atom_1.element() != atom_2.element() )
            throw std::runtime_error( "Elements are not the same." );
        // Check that their positions are not too far apart--if they are, we should probably use the nearest copy instead.
        double shortest_distance;
        Vector3D difference_vector;
        crystal_structure.shortest_distance( atom_1.position(), atom_2.position(), shortest_distance, difference_vector );
        // Must transform to Cartesian coordinates
        if ( ! nearly_equal( shortest_distance, ( atom_2_Cartesian - atom_1_Cartesian ).length() ) )
            std::cout << "Warning: shortest distance != distance" << std::endl;
        Vector3D new_position = ( atom_1.position() + atom_2.position() ) / 2.0;
        // We have to wipe all their attributes such as ADPs or atom label.
        std::string new_label( atom_1_label );
        if ( ( ! atom_1_label.empty() ) || ( ! atom_2_label.empty() ) )
        {
            if ( ( to_upper( atom_1_label[ atom_1_label.length() - 1 ] ) == 'A' ) && ( to_upper( atom_2_label[ atom_2_label.length() - 1 ] ) == 'B' ) )
            {
                if ( atom_1_label.substr( 0, atom_1_label.length() - 1 ) == atom_2_label.substr( 0, atom_2_label.length() - 1 ) )
                    new_label = atom_1_label.substr( 0, atom_1_label.length() - 1 );
            }
        }
        Atom new_atom( atom_1.element(), new_position, new_label );
        crystal_structure.set_atom( crystal_structure.atom( atom_1_label ), new_atom );
        // Delete the second atom.
        crystal_structure.set_suppressed( crystal_structure.atom( atom_2_label ), true );
        }
        crystal_structure.save_cif( append_to_file_name( input_file_name, "_avg" ) );
    MACRO_END_GAME

    try // Average two unit cells.
    {
        if ( argc != 3 )
            throw std::runtime_error( "Please give the names of two .cif files." );
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

    try // Write lean.
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

    try // Write CASTEP input files for optimising H atoms (unit cell fixed, non-H fixed).
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        WriteCASTEPFile write_CASTEP_file( crystal_structure, input_file_name.directory(), input_file_name.name() );
        write_CASTEP_file.set_job_type( WriteCASTEPFile::H_ATOMS_ONLY );
        write_CASTEP_file.write();
    MACRO_END_GAME

    try // Test DoubleWithESD.
    {
        DoubleWithESD dwe_1( "10.81(8)" );
        double conversion_factor = 1.0/( 8.0 *CONSTANT_PI * CONSTANT_PI );
        DoubleWithESD dwe_2( conversion_factor * dwe_1.value(), conversion_factor * dwe_1.estimated_standard_deviation() );
        std::cout << dwe_2.crystallographic_style() << std::endl;
    MACRO_END_GAME

    try // Test PowderPatternCalculator::calculate_equivalent_reflections().
    {
        FileName input_file_name( "P32_No_145.cif" );
        CrystalStructure crystal_structure;
        read_cif( input_file_name, crystal_structure );
        PowderPatternCalculator powder_pattern_calculator( crystal_structure );
        powder_pattern_calculator.calculate_reflection_list();
        ReflectionList reflection_list = powder_pattern_calculator.reflection_list();
        reflection_list.save( replace_extension( input_file_name, "txt" ) );
    MACRO_END_GAME

    try // Analyse volumes in .cif file before and after minimisation.
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
            text_file_writer.write_line( exp_cif_file.name() + " " + double2string( exp_volume ) + " " + double2string( min_volume ) );
        }
    MACRO_END_GAME

    try // Test crystallographic_style().
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

    try // Test running average and ESD.
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

    try // Analyse unit-cell parameter files from Materials Studio.
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

    try // Apply space-group symmetry.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        crystal_structure.apply_space_group_symmetry();
        crystal_structure.save_cif( append_to_file_name( input_file_name, "_asgs" ) );
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

    try // De Gelder's normalised cross correlation.
    {
        if ( argc != 3 )
            throw std::runtime_error( "Please give the names of two .cif or .cell files." );
        FileName file_name_1( argv[ 1 ] );
        CrystalStructure crystal_structure_1;
        read_cif_or_cell( file_name_1, crystal_structure_1 );
        crystal_structure_1.apply_space_group_symmetry();
        FileName file_name_2( argv[ 2 ] );
        CrystalStructure crystal_structure_2;
        read_cif_or_cell( file_name_2, crystal_structure_2 );
        crystal_structure_2.apply_space_group_symmetry();
        Angle two_theta_start( 3.0, Angle::DEGREES );
        Angle two_theta_end(  35.0, Angle::DEGREES );
        Angle two_theta_step( 0.01, Angle::DEGREES );
        double FWHM( 0.1 );
        PowderPatternCalculator powder_pattern_calculator_1( crystal_structure_1 );
        powder_pattern_calculator_1.set_two_theta_start( two_theta_start );
        powder_pattern_calculator_1.set_two_theta_end( two_theta_end );
        powder_pattern_calculator_1.set_two_theta_step( two_theta_step );
        powder_pattern_calculator_1.set_FWHM( FWHM );
        PowderPattern powder_pattern_1;
        powder_pattern_calculator_1.calculate( powder_pattern_1 );
        PowderPatternCalculator powder_pattern_calculator_2( crystal_structure_2 );
        powder_pattern_calculator_2.set_two_theta_start( two_theta_start );
        powder_pattern_calculator_2.set_two_theta_end( two_theta_end );
        powder_pattern_calculator_2.set_two_theta_step( two_theta_step );
        powder_pattern_calculator_2.set_FWHM( FWHM );
        PowderPattern powder_pattern_2;
        powder_pattern_calculator_2.calculate( powder_pattern_2 );
        double nwcc = normalised_weighted_cross_correlation( powder_pattern_1, powder_pattern_2, Angle( 3.0, Angle::DEGREES ) );
        std::cout << "Normalised cross correlation = " << nwcc << std::endl;
    MACRO_END_GAME

    try // Beautify experimental structure for energy minimisation.
    {
        MACRO_ONE_CIFFILENAME_AS_ARGUMENT
        // Centred to primitive.
        if ( ! crystal_structure.space_group().centring().is_primitive() )
            crystal_structure.reduce_to_primitive();
        // Unit-cell angles as close as possible to 90.
        crystal_structure.transform( crystal_structure.crystal_lattice().choose_angles_close_to_90() );
        crystal_structure.move_com_close_to_origin();
        // Normalise X-H bonds.
        normalise_C_F_bonds( crystal_structure );
        normalise_X_H_bonds( crystal_structure );
        crystal_structure.save_cif( replace_extension( append_to_file_name( input_file_name, "_beautified" ), "cif" ) );
    MACRO_END_GAME

    try // Calculate normalised cross correlation for two cif files.
    {
        Angle two_theta_start( 5.0, Angle::DEGREES );
        Angle two_theta_end(  50.0, Angle::DEGREES );
        Angle two_theta_step( 0.01, Angle::DEGREES );
        double FWHM( 0.1 );

        CrystalStructure crystal_structure_1;
        read_cif( FileName( "C:\\Data\\Refereeing\\GF.cif" ), crystal_structure_1 );
        PowderPatternCalculator powder_pattern_calculator_1( crystal_structure_1 );
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

    try // WUBDOM.
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

    try // WIMWOE.
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
        TextFileWriter text_file_writer( FileName( "WIMWOE_output.txt" ) );
        for ( size_t i( 0 ); i != positions.size(); ++i )
        {
            Vector3D result = crystal_lattice.orthogonal_to_fractional( positions[ i ] );
            std::cout << result << std::endl;
            text_file_writer.write_line( double2string( result.x() ) + " " + double2string( result.y() ) + " " + double2string( result.z() ) );
        }
    MACRO_END_GAME

    try // TRIZIN04, crystal structures of s-triazine.
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
        crystal_structure.save_cif( FileName( "TRIZIN04.cif" ) );
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
        crystal_structure.save_cif( FileName( "TRIZIN04.cif" ) );
    MACRO_END_GAME

    try // TRIZIN05, crystal structures of s-triazine.
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
        crystal_structure.save_cif( FileName( "TRIZIN05.cif" ) );
    MACRO_END_GAME

    try // XIJKIK.
    {
        FileName input_file_name( "XIJKIK.cif" );
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

    try // Compile all headers stand-alone.
    {
        FileName file_list_file_name( "FileList.txt" );
        FileList file_list( file_list_file_name );
        file_list.set_prepend_file_name_with_basedirectory( false );
        if ( file_list.empty() )
            throw std::runtime_error( std::string( "No files in file list " ) + file_list_file_name.full_name() );
        for ( size_t i( 0 ); i != file_list.size(); ++i )
        {
            TextFileWriter text_file_writer( FileName( "", file_list.value( i ).name() + "_header", "cpp" ) );
            text_file_writer.write_line( "#include \"" + file_list.value( i ).name() + ".h\"" );
        }
    MACRO_END_GAME

    try // What is a good powder pattern.
    {
        MACRO_ONE_FILELISTNAME_AS_ARGUMENT
        FileName input_file_name( file_list_file_name );
        std::cout << "Number of structures read: " << file_list.size() << std::endl;
        // Weed out structures that look too much like another structure.
        double similarity_limit = 0.95;
        file_list = select_diverse_structures( file_list, similarity_limit );
        std::cout << "Number of structures left: " << file_list.size() << std::endl;
        bool save_xye_files( false );
        bool include_background = true;
        double background_level = 0.2;
        bool include_PO = true;
        bool include_FCJ = true;
        bool include_NaCl = false;
        bool include_noise = true;
        bool set_Uiso = true;
        double Uiso = 0.0375; // 0.0375; // 0.125;
        // We have a bit of a problem here. We simulate zero-point error by just subtracting the zero point from all
        // 2theta values, but after we have done that, we cannot compare the "experimental" powder pattern to an ideal
        // powder pattern any more because the 2theta ranges are different.
        // The easiest fix is to make the zero-point error an integer times the step size and then to shift things back and forth a couple of times.
        Angle two_theta_start( 5.0, Angle::DEGREES );
        Angle two_theta_end(  35.0, Angle::DEGREES );
        Angle two_theta_step( 0.015, Angle::DEGREES );
        size_t zero_point_error_in_step_sizes = 0;
        double zero_point_error = zero_point_error_in_step_sizes * two_theta_step.value_in_degrees(); // 0.02 // 0.06;
        two_theta_start -= zero_point_error_in_step_sizes * two_theta_step;
        two_theta_end   += zero_point_error_in_step_sizes * two_theta_step;
        double Bragg_total_signal_normalisation = 10000.0;
        double PO_extent = 0.9; // 0.9 // 0.7;
        double FWHM = 0.2; // 0.2; // ?
        double A = 10.0 / 400.0; // Finger-Cox-Jephcoat.
        double highest_peak = 10000.0; // 10000.0 // 300.0;
        bool add_impurity = false;
        double impurity_total_signal_normalisation = 0.2 * Bragg_total_signal_normalisation;
        double background_total_signal_normalisation = background_level * Bragg_total_signal_normalisation;
        double NaCl_total_signal_normalisation = 0.002 * Bragg_total_signal_normalisation;
        std::vector< PowderPattern > powder_patterns;
        // Precalcuate all perfect powder diffraction patterns.
        for ( size_t iCandidateStructure( 0 ); iCandidateStructure != file_list.size(); ++iCandidateStructure )
        {
                PowderPattern powder_pattern;
                CrystalStructure crystal_structure;
                std::cout << "Now reading cif... " + file_list.value( iCandidateStructure ).full_name() << std::endl;
                read_cif( file_list.value( iCandidateStructure ), crystal_structure );
                crystal_structure.apply_space_group_symmetry();
                PowderPatternCalculator powder_pattern_calculator( crystal_structure );
                powder_pattern_calculator.set_two_theta_start( two_theta_start );
                powder_pattern_calculator.set_two_theta_end( two_theta_end );
                powder_pattern_calculator.set_two_theta_step( two_theta_step );
                powder_pattern_calculator.set_FWHM( FWHM );
                powder_pattern_calculator.calculate( powder_pattern );
                powder_patterns.push_back( powder_pattern );
        }
        // Scan a parameter.
        for ( size_t iParVal( 0 ); iParVal != 3; ++iParVal )
        {
            std::cout << "iParVal = " << iParVal << std::endl;
            PO_extent = 0.9 - iParVal * 0.15;
            for ( size_t iCorrectStructure( 0 ); iCorrectStructure != file_list.size(); ++iCorrectStructure )
            {
                std::cout << "iCorrectStructure = " << iCorrectStructure << std::endl;
                CrystalStructure crystal_structure;
                std::cout << "Now reading cif... " + file_list.value( iCorrectStructure ).full_name() << std::endl;
                read_cif( file_list.value( iCorrectStructure ), crystal_structure );
                // ######################## CHANGE THIS ##################################
                std::string file_name_ID( "_NoNaCl" );
                MillerIndices PO_direction = MillerIndices( 0, 0, 1 );
                if ( set_Uiso )
                    crystal_structure.set_global_Uiso( Uiso );
                if ( include_PO )
                {
                    switch ( crystal_structure.crystal_lattice().lattice_system() )
                    {
                        case CrystalLattice::TRICLINIC    :
                        case CrystalLattice::ORTHORHOMBIC : {
                                                                CrystalLattice crystal_lattice = crystal_structure.crystal_lattice();
                                                                if ( crystal_lattice.a() < crystal_lattice.b() )
                                                                {
                                                                    if ( crystal_lattice.a() < crystal_lattice.c() )
                                                                        PO_direction = MillerIndices( 1, 0, 0 );
                                                                }
                                                                else // b < a
                                                                {
                                                                    if ( crystal_lattice.b() < crystal_lattice.c() )
                                                                        PO_direction = MillerIndices( 0, 1, 0 );
                                                                }
                                                                break;
                                                            }
                        case CrystalLattice::MONOCLINIC_A : {
                                                                PO_direction = MillerIndices( 1, 0, 0 );
                                                                break;
                                                            }
                        case CrystalLattice::MONOCLINIC_B : {
                                                                PO_direction = MillerIndices( 0, 1, 0 );
                                                                break;
                                                            }
                        case CrystalLattice::MONOCLINIC_C : {
                                                                PO_direction = MillerIndices( 0, 0, 1 );
                                                                break;
                                                            }
                        case CrystalLattice::TETRAGONAL   :
                        case CrystalLattice::HEXAGONAL    : break; // Nothing to do
                        case CrystalLattice::RHOMBOHEDRAL :
                        case CrystalLattice::CUBIC        : include_PO = false;
                    }
                }
                crystal_structure.apply_space_group_symmetry();
                std::cout << "Now calculating powder pattern... " << std::endl;
                PowderPatternCalculator powder_pattern_calculator( crystal_structure );
                
   //     double zero_point_error = zero_point_error_in_step_sizes * two_theta_step.value_in_degrees(); // 0.02 // 0.06;
                two_theta_start += zero_point_error_in_step_sizes * two_theta_step;
                two_theta_end   += zero_point_error_in_step_sizes * two_theta_step;
                
                powder_pattern_calculator.set_two_theta_start( two_theta_start );
                powder_pattern_calculator.set_two_theta_end( two_theta_end );
                powder_pattern_calculator.set_two_theta_step( two_theta_step );
                powder_pattern_calculator.set_FWHM( FWHM );
                if ( include_PO )
                    powder_pattern_calculator.set_preferred_orientation( PO_direction, PO_extent );
                double B = A; // Finger-Cox-Jephcoat.
                if ( include_FCJ )
                    powder_pattern_calculator.set_finger_cox_jephcoat( FingerCoxJephcoat( A, B ) );
                PowderPattern powder_pattern_Bragg_diffraction;
                powder_pattern_calculator.calculate( powder_pattern_Bragg_diffraction );
                powder_pattern_Bragg_diffraction.normalise_total_signal( Bragg_total_signal_normalisation );
                powder_pattern_Bragg_diffraction.correct_zero_point_error( Angle::from_degrees( -zero_point_error ) );
                PowderPattern result = powder_pattern_Bragg_diffraction;
                PowderPattern powder_pattern_NaCl;
                if ( include_NaCl )
                {
                    CrystalStructure crystal_structure_NaCl = NaCl();
                    crystal_structure_NaCl.apply_space_group_symmetry();
                    PowderPatternCalculator NaCl_powder_pattern_calculator( crystal_structure_NaCl );
                    NaCl_powder_pattern_calculator.set_two_theta_start( two_theta_start );
                    NaCl_powder_pattern_calculator.set_two_theta_end( two_theta_end );
                    NaCl_powder_pattern_calculator.set_two_theta_step( two_theta_step );
                    NaCl_powder_pattern_calculator.set_FWHM( 0.1 );
                    // We never include PO for NaCl (it is cubic).
                    NaCl_powder_pattern_calculator.calculate( powder_pattern_NaCl );
                    powder_pattern_NaCl.normalise_total_signal( NaCl_total_signal_normalisation );
                    powder_pattern_NaCl.correct_zero_point_error( Angle::from_degrees( -zero_point_error ) );
                    result += powder_pattern_NaCl;
                }
                PowderPattern powder_pattern_background;
                if ( include_background )
                {
                    PowderPatternCalculator background_powder_pattern_calculator( crystal_structure );
                    background_powder_pattern_calculator.set_two_theta_start( two_theta_start );
                    background_powder_pattern_calculator.set_two_theta_end( two_theta_end );
                    background_powder_pattern_calculator.set_two_theta_step( two_theta_step );
                    background_powder_pattern_calculator.set_FWHM( 5.0 );
                    // We never include PO for the amorphous background
                    background_powder_pattern_calculator.calculate( powder_pattern_background );
                    powder_pattern_background.normalise_total_signal( background_total_signal_normalisation );
                    powder_pattern_background.correct_zero_point_error( Angle::from_degrees( -zero_point_error ) );
                }
                if ( add_impurity )
                {
                    for ( size_t j( 0 ); j != file_list.size(); ++j )
                    {
                        PowderPattern powder_pattern_2;
                        CrystalStructure crystal_structure_2;
                        std::cout << "Now reading cif... " + file_list.value( j ).full_name() << std::endl;
                        read_cif( file_list.value( j ), crystal_structure_2 );
                        crystal_structure_2.apply_space_group_symmetry();
                        PowderPatternCalculator powder_pattern_calculator_2( crystal_structure_2 );
                        powder_pattern_calculator_2.set_two_theta_start( two_theta_start );
                        powder_pattern_calculator_2.set_two_theta_end( two_theta_end );
                        powder_pattern_calculator_2.set_two_theta_step( two_theta_step );
                        powder_pattern_calculator_2.set_FWHM( FWHM );
                        // Include PO?
                        powder_pattern_calculator_2.calculate( powder_pattern_2 );
                        powder_pattern_2.normalise_total_signal( impurity_total_signal_normalisation );
                        powder_pattern_2.correct_zero_point_error( Angle::from_degrees( -zero_point_error ) );
                        result += powder_pattern_2;
                    }
                }
                double scale_factor = result.normalise_highest_peak( highest_peak );
                if ( include_NaCl )
                {
                    powder_pattern_NaCl.scale( scale_factor );
                    powder_pattern_NaCl.make_counts_integer();
                    powder_pattern_NaCl.recalculate_estimated_standard_deviations();
                    if ( save_xye_files )
                        powder_pattern_NaCl.save_xye( replace_extension( append_to_file_name( input_file_name, file_name_ID + "_NaCl" ), "xye" ), true );
                }
                if ( include_background )
                {
                    powder_pattern_background.scale( scale_factor );
                    powder_pattern_background.make_counts_integer();
                    powder_pattern_background.add_constant_background( 20.0 );
                    powder_pattern_background.recalculate_estimated_standard_deviations();
                    if ( save_xye_files )
                        powder_pattern_background.save_xye( replace_extension( append_to_file_name( input_file_name, file_name_ID + "_BKGR" ), "xye" ), true );
                }
                powder_pattern_Bragg_diffraction.scale( scale_factor );
                powder_pattern_Bragg_diffraction.make_counts_integer();
                powder_pattern_Bragg_diffraction.recalculate_estimated_standard_deviations();
                if ( save_xye_files )
                    powder_pattern_Bragg_diffraction.save_xye( replace_extension( append_to_file_name( input_file_name, file_name_ID + "_Bragg" ), "xye" ), true );
                result = powder_pattern_Bragg_diffraction;
                if ( include_NaCl )
                    result += powder_pattern_NaCl;
                if ( include_background )
                    result += powder_pattern_background;
                if ( include_noise )
                {
                    PowderPattern powder_pattern_noise = calculate_Poisson_noise( result );
                    if ( save_xye_files )
                        powder_pattern_noise.save_xye( replace_extension( append_to_file_name( input_file_name, file_name_ID + "_noise" ), "xye" ), true );
                    result += powder_pattern_noise;
                }
                // They are *estimated* standard deviations, so they should be calculated *after* the noise has been introduced.
                result.recalculate_estimated_standard_deviations();
                if ( save_xye_files )
                    result.save_xye( replace_extension( append_to_file_name( input_file_name, file_name_ID + "_cal_XRPD" ), "xye" ), true );
                if ( true )
                {
                    PowderPattern estimated_background = calculate_Brueckner_background( result,
                                                                                         50, // niterations
                                                                                         round_to_int( 50.0 * ( Angle::from_degrees( 0.015 ) / result.average_two_theta_step() ) ), // window
                                                                                         true, // apply_smoothing
                                                                                         5 ); // smoothing_window
                    if ( save_xye_files )
                        estimated_background.save_xye( replace_extension( append_to_file_name( input_file_name, file_name_ID + "_Brueckner_BKGR" ), "xye" ), true );
                    result -= estimated_background;
                    if ( save_xye_files )
                        result.save_xye( replace_extension( append_to_file_name( input_file_name, file_name_ID + "_Brueckner_BKGR_subtracted" ), "xye" ), true );
                }
                
                
         //       two_theta_start -= zero_point_error_in_step_sizes * two_theta_step;
         //       two_theta_end   += zero_point_error_in_step_sizes * two_theta_step;
                
                // One of these is superfluous -- figure out which one is the one caused by sample displacement
                result.set_two_theta_start( two_theta_start + zero_point_error_in_step_sizes * two_theta_step );
                result.set_two_theta_end( two_theta_end - zero_point_error_in_step_sizes * two_theta_step );
                // When we are here, we have the distorted "experimental" powder diffraction pattern.
                size_t best_match = iCorrectStructure;
                double highest_correlation = 0.0;
                double correct_correlation = 0.0;
                for ( size_t iCandidateStructure( 0 ); iCandidateStructure != file_list.size(); ++iCandidateStructure )
                {
                    double correlation = normalised_weighted_cross_correlation( result, powder_patterns[ iCandidateStructure ]);
                    if ( iCandidateStructure == iCorrectStructure )
                        correct_correlation = correlation;
                    else
                    {
                        if ( correlation > highest_correlation )
                        {
                            highest_correlation = correlation;
                            best_match = iCandidateStructure;
                        }
                    }
                }
                std::cout << "iCorrectStructure = " << iCorrectStructure << std::endl;
                std::cout << "Correlation for iCorrectStructure = " << correct_correlation << std::endl;
                std::cout << "best match among decoys = " << best_match << std::endl;
                std::cout << "highest_correlation = " << highest_correlation << std::endl;
            }
        }
    MACRO_END_GAME

}
