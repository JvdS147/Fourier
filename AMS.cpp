#include "AMS.h"

#include "AMS_Convert_flx2xyz.h"
#include "CorrelationMatrix.h"
#include "CrystalStructure.h"
#include "EndGame.h"
#include "Histogram.h"
#include "MathFunctions.h"
#include "PowderPattern.h"
#include "PowderPatternCalculator.h"
#include "Pressure.h"
#include "ReadCif.h"
#include "SimilarityAnalysis.h"
#include "Sort.h"
#include "TextFileReader.h"
#include "TextFileWriter.h"
#include "Utilities.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

int AMS_main( int argc, char** argv )
{

    try // Report generator.
    {
        ReportGenerator report_generator( FileName( "X:\\GC\\FreeEnergyCorrection\\table.txt" ),
//        ReportGenerator report_generator( FileName( "C:\\Data_Win\\ContractResearch\\GC\\p0004\\FreeEnergyCorrection_tautomer_2\\table.txt" ),
//        ReportGenerator report_generator( FileName( "\\\\Mac\\Home\\Documents\\Data_mac\\ContractResearch\\GC\\LargeVolume_02\\FreeEnergyCorrection\\table.txt" ),
                                          Fraction( 44, 0, 1 ) );
        report_generator.generate_report();
//        report_generator.generate_some_numbers();
    MACRO_END_GAME

    try // Convert .flx to .xyz.
    {
        if ( argc != 2 )
            throw std::runtime_error( "Please give the name of a .flx file." );
        std::string input_file_name = argv[ 1 ];
        convert_flx2xyz( input_file_name );
    }
    catch ( std::exception & e )
    {
        std::cout << "An exception was thrown" << std::endl;
        std::cout << e.what() << std::endl;
        char a;
        std::cin >> a;
    }
    return 0;
        
}

namespace
{

//    item_2 {
//        forceFieldSigma : 0.3396833849
//        freeEnergySigma : 0.0064124333754743386
//        actualGenerationConvergence : 0.95003499999999996
//        actualRerankingConvergence : 0.95164863909999997

double get_value_item_keyword( const TextFileReader_2 & input_file, const size_t i, const std::string & keyword )
{
    std::string line = input_file.line( input_file.find_whole_word( keyword, input_file.find_whole_word( "item_" + size_t2string( i ) ) ) );
    std::vector< std::string > words = split( line );
    if ( words.size() != 3 )
        throw std::runtime_error( "get_value_item_keyword(): words.size() != 3 for keyword " + keyword );
    if ( words[1] != ":" )
        throw std::runtime_error( "get_value_item_keyword(): words[2] != \":\" for keyword " + keyword );
    return string2double( words[2] );
}

} // namespace

// ********************************************************************************

double bar_to_kiloPascal( double const pressure )
{
    return pressure * 100.0;
}

// ********************************************************************************

double kJ_to_kcal( double const energy )
{
    return energy / 4.185;
}

// ********************************************************************************

ReportGenerator::ReportGenerator( const FileName & table_txt_file_name, const Fraction natoms_per_molecule ):
natoms_per_molecule_(natoms_per_molecule)
{
    
// rank |       energy        |      density       |       volume       | atoms_per_asym_unit |  space_group  |         a          |         b          |         c          |       alpha        |        beta        |       gamma        |  ID
//      |   [kcal/mol/atom]   |      [g/cm3]       |        [A3]        |                     |               |        [A]         |        [A]         |        [A]         |     [degrees]      |     [degrees]      |     [degrees]      |
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//  1   | -152.55915259239009 | 1.5461175279492867 | 629.92249024307955 |         66          |      P_1      | 6.4408878052983169 | 8.687209046544055  | 11.450422511528654 | 96.344185771704446 | 83.977323219058448 | 96.476912398511757 | 34776
//  2   | -152.55348831829397 | 1.5060333128637009 | 2586.7534140061393 |         66          | P_2_1_2_1_2_1 | 11.931048062019025 | 13.461674721448558 | 16.105616120059551 |         90         |         90         |         90         | 65048
//  3   | -152.55150132620625 | 1.5675308464639806 | 1242.634817185613  |         66          |     P_2_1     | 12.534979242534275 | 8.6025362933016911 | 12.899742269317013 |         90         | 116.70516367025786 |         90         | 34777

    directory_ = table_txt_file_name.directory();
    FileName temp_file_name = FileName( directory_, "ConvergenceInformation", "txt" );
    if ( temp_file_name.exists() )
    {
        ConvergenceItem_txt_ = TextFileReader_2( temp_file_name );
        ConvergenceInformation_file_found_ = true;
    }
    else
    {
        temp_file_name = FileName( directory_, "ConvergenceInformationItem", "txt" );
        if ( temp_file_name.exists() )
        {
            ConvergenceItem_txt_ = TextFileReader_2( temp_file_name );
            ConvergenceInformation_file_found_ = true;
        }
        else
            ConvergenceInformation_file_found_ = false;
    }
    if ( natoms_per_molecule.fractional_part() == Fraction( 1, 2 ) )
    {
        table_file_ = AMS_TableFile( table_txt_file_name, (natoms_per_molecule * 2).integer_part() );
    }
    else
        table_file_ = AMS_TableFile( table_txt_file_name, natoms_per_molecule.integer_part() );
    FileName file_list_name( directory_ + "structures", "FileList", "txt" );
    TextFileWriter file_list_writer( file_list_name );
    for ( size_t i( 0 ); i != table_file_.size(); ++i )
        file_list_writer.write_line( "structure_" + size_t2string( i+1, 6, '0' ) + ".cif" );
    file_list_writer.~TextFileWriter();
    file_list_.initialise_from_file( file_list_name );
}

// ********************************************************************************

void ReportGenerator::generate_report() const
{

//   table_txt_
// rank |       energy        |      density       |       volume       | atoms_per_asym_unit |  space_group  |         a          |         b          |         c          |       alpha        |        beta        |       gamma        |  ID
//      |   [kcal/mol/atom]   |      [g/cm3]       |        [A3]        |                     |               |        [A]         |        [A]         |        [A]         |     [degrees]      |     [degrees]      |     [degrees]      |
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//  1   | -152.55915259239009 | 1.5461175279492867 | 629.92249024307955 |         66          |      P_1      | 6.4408878052983169 | 8.687209046544055  | 11.450422511528654 | 96.344185771704446 | 83.977323219058448 | 96.476912398511757 | 34776
//  2   | -152.55348831829397 | 1.5060333128637009 | 2586.7534140061393 |         66          | P_2_1_2_1_2_1 | 11.931048062019025 | 13.461674721448558 | 16.105616120059551 |         90         |         90         |         90         | 65048
//  3   | -152.55150132620625 | 1.5675308464639806 | 1242.634817185613  |         66          |     P_2_1     | 12.534979242534275 | 8.6025362933016911 | 12.899742269317013 |         90         | 116.70516367025786 |         90         | 34777

    TextFileWriter output_table( FileName( directory_, "table_2", "txt" ) );
    for ( size_t i( 0 ); i != table_file_.size(); ++i )
    {
        std::string line;
        line += double2string_2( table_file_.entry(i).density(), 3 );
        line += " ";
        line += double2string_2( ( table_file_.entry(i).energy() - table_file_.lowest_energy() ) * natoms_per_molecule_, 3 );
        line += " ";
        line += double2string_2(  table_file_.entry(i).density(), 3 );
        line += " ";
        if ( Fraction( table_file_.entry(i).natoms() ) == natoms_per_molecule_ )
        {
            line += "1";
        }
        else if ( Fraction( table_file_.entry(i).natoms() ) == ( 2 * natoms_per_molecule_ ) )
        {
            line += "2";
        }
        else
        {
            line += double2string( table_file_.entry(i).natoms() / natoms_per_molecule_ );
        }
        line += " ";
        line += remove( table_file_.entry(i).space_group(), '_' );
        line += " ";
        line += double2string_2( table_file_.entry(i).crystal_lattice().a(), 3 );
        line += " ";
        line += double2string_2( table_file_.entry(i).crystal_lattice().b(), 3 );
        line += " ";
        line += double2string_2( table_file_.entry(i).crystal_lattice().c(), 3 );
        line += " ";
        std::string alpha_str = double2string_2( table_file_.entry(i).crystal_lattice().alpha().value_in_degrees(), 3 );
        if ( alpha_str == "90.000" )
            alpha_str = "90";
        if ( alpha_str == "120.000" )
            alpha_str = "120";
        line += alpha_str;
        line += " ";
        std::string beta_str = double2string_2( table_file_.entry(i).crystal_lattice().beta().value_in_degrees(), 3 );
        if ( beta_str == "90.000" )
            beta_str = "90";
        if ( beta_str == "120.000" )
            beta_str = "120";
        line += beta_str;
        line += " ";
        std::string gamma_str = double2string_2( table_file_.entry(i).crystal_lattice().gamma().value_in_degrees(), 3 );
        if ( gamma_str == "90.000" )
            gamma_str = "90";
        if ( gamma_str == "120.000" )
            gamma_str = "120";
        line += gamma_str;
        line += " ";
        output_table.write_line( line );
    }
    
    // ############# Energy density histogram ########################
    Histogram histogram( 0.0, 1.0, 10 );
    for ( size_t i( 0 ); i != table_file_.size(); ++i )
        histogram.add_data( ( table_file_.entry( i ).energy() - table_file_.lowest_energy() ) * table_file_.natoms_per_molecule() );
    TextFileWriter text_file_writer_histogram( FileName( directory_, "EnergyHistogram", "txt" ) );
    text_file_writer_histogram.write_line( "Number of polymorphs in bottom 1 kcal/mol: " + double2string( table_file_.npolymorphs_in_bottom( 1.0 ) ) );
    for ( size_t i( 0 ); i != histogram.number_of_bins(); ++i )
        text_file_writer_histogram.write_line( size_t2string( histogram.bin( i ) ) );

    // ############# Insert the correct data name into the cif files ########################
    if ( file_list_.empty() )
        throw std::runtime_error( "No files in file list " );
    for ( size_t i( 0 ); i != file_list_.size(); ++i )
    {
        if ( to_lower( file_list_.value( i ).extension() ) != "cif" )
            std::cout << "Warning: structure does not have .cif extension." << std::endl;
        TextFileReader_2 input_file( file_list_.value( i ) );
        TextFileWriter output_file( file_list_.value( i ) );
        for ( size_t iLine( 0 ); iLine != input_file.size(); ++iLine )
        {
            if ( input_file.line( iLine ).substr( 0, 5 ) == "data_" )
                output_file.write_line( "data_" + size_t2string( i + 1, 6, '0' ) );
            else
                output_file.write_line( input_file.line( iLine ) );
        }
    }
    // ############# Similarity matrix ########################
    CorrelationMatrix correlation_matrix = calculate_correlation_matrix( file_list_ );
    correlation_matrix.save( FileName( directory_, "SimilarityMatrix", "txt" ) );
    // ############# some text / convergence item ########################
    generate_some_numbers();
}

// ********************************************************************************

void ReportGenerator::generate_some_numbers() const
{
//ConvergenceInformationItem {
//    item_0 {
//        generationStatus : CALCULATION_CONVERGED
//        rerankingStatus : CALCULATION_CONVERGED
//        freeEnergyCorrectionStatus : CALCULATION_CONVERGED
//        forceFieldSigma : 0.30161199480000001
//        freeEnergySigma : 0.0011315790116836974
//        actualGenerationConvergence : 0.99026099999999995
//        actualRerankingConvergence : 0.9948776346
//        actualGenerationWindowExtension : 3.3133849993633264
//        actualRerankingWindowExtension : 3.0000000000000031
//        nbGeneratedStructures : 1747
//        nbGeneratedStableStructures : 951
//        nbRerankedStructures : 357
//        nbFreeEnergyCorrectedStructures : 24
//        nbRerankFailures : 0
//        nbFreeEnergyCorrectionFailures : 0
//        cumulatedNormalizedCPUTime : 458597.85978200001
//    }
//    item_1 {
//        generationStatus : CALCULATION_CONVERGED
//        rerankingStatus : CALCULATION_CONVERGED
//        freeEnergyCorrectionStatus : CALCULATION_CONVERGED
//        forceFieldSigma : 0.3790919478
//        freeEnergySigma : 0.002983044594846006
//        actualGenerationConvergence : 0.95013000000000003
//        actualRerankingConvergence : 0.95001721620000001
//        actualGenerationWindowExtension : 2.7846581689544085
//        actualRerankingWindowExtension : 0.9483423526745095
//        nbGeneratedStructures : 9143
//        nbGeneratedStableStructures : 4150
//        nbRerankedStructures : 397
//        nbFreeEnergyCorrectedStructures : 15
//        nbRerankFailures : 0
//        nbFreeEnergyCorrectionFailures : 0
//        cumulatedNormalizedCPUTime : 1051550.7968919999
//    }
//    item_2 {
//        generationStatus : CALCULATION_CONVERGED
//        rerankingStatus : CALCULATION_CONVERGED
//        freeEnergyCorrectionStatus : CALCULATION_CONVERGED
//        forceFieldSigma : 0.3396833849
//        freeEnergySigma : 0.0064124333754743386
//        actualGenerationConvergence : 0.95003499999999996
//        actualRerankingConvergence : 0.95164863909999997
//        actualGenerationWindowExtension : 2.9127545594930737
//        actualRerankingWindowExtension : 0.80207921190006592
//        nbGeneratedStructures : 9890
//        nbGeneratedStableStructures : 5478
//        nbRerankedStructures : 485
//        nbFreeEnergyCorrectedStructures : 15
//        nbRerankFailures : 0
//        nbFreeEnergyCorrectionFailures : 0
//        cumulatedNormalizedCPUTime : 3297114.2460480002
//    }
//}

    if ( ! ConvergenceInformation_file_found_ )
    {
        std::cout << "ReportGenerator::generate_some_numbers(): WARNING: ConvergenceInformation.txt and ConvergenceInformationItem.txt not found, some numbers file cannot be generated. " << std::endl;
        return;
    }
    size_t iPos_item_0 = ConvergenceItem_txt_.find_whole_word( "item_0" );
    bool item_0_found = ( iPos_item_0 != std::string::npos );
    size_t iPos_item_1 = ConvergenceItem_txt_.find_whole_word( "item_1" );
    bool item_1_found = ( iPos_item_1 != std::string::npos );
    size_t iPos_item_2 = ConvergenceItem_txt_.find_whole_word( "item_2" );
    bool item_2_found = ( iPos_item_2 != std::string::npos );
    
//    item_1 {
//        generationStatus : CALCULATION_TIME_OUT
//        rerankingStatus : CALCULATION_CONVERGED
//        freeEnergyCorrectionStatus : CALCULATION_CONVERGED
//
//
//    item_2 {
//        generationStatus : CALCULATION_NOT_STARTED
//        rerankingStatus : CALCULATION_NOT_STARTED
//        freeEnergyCorrectionStatus : CALCULATION_NOT_STARTED

    size_t iPos_CALCULATION_TIME_OUT = ConvergenceItem_txt_.find_whole_word( "generationStatus : CALCULATION_TIME_OUT" );
    bool CALCULATION_TIME_OUT_found = ( iPos_CALCULATION_TIME_OUT != std::string::npos );

    size_t iPos_CALCULATION_NOT_STARTED = ConvergenceItem_txt_.find_whole_word( "generationStatus : CALCULATION_NOT_STARTED" );
    bool CALCULATION_NOT_STARTED_found = ( iPos_CALCULATION_NOT_STARTED != std::string::npos );
    
    if ( CALCULATION_NOT_STARTED_found )
    {
        if ( ! CALCULATION_TIME_OUT_found )
            std::cout << "ReportGenerator::generate_some_numbers(): Warning: generation not started but previous generation not timed out." << std::endl;
        else
        {
            if ( ( iPos_CALCULATION_NOT_STARTED > iPos_item_2 ) && ( iPos_CALCULATION_TIME_OUT > iPos_item_1 ) )
                item_2_found = false;
        }
    }

    if ( item_0_found )
    {
        if ( item_1_found )
        {
            if ( item_2_found )
            {
                TextFileWriter some_numbers_file( FileName( directory_, "SomeNumbers", "txt" ) );
                some_numbers_file.write_line( "The accuracy of the TMFF measured independently for Z'=1, Z'=2 first round and Z'=2 second round turned out to be " +
                double2string_2( get_value_item_keyword( ConvergenceItem_txt_, 0, "forceFieldSigma" ), 3 ) + " kcal/mol/atom, " +
                double2string_2( get_value_item_keyword( ConvergenceItem_txt_, 1, "forceFieldSigma" ), 3 ) + " kcal/mol/atom and " +
                double2string_2( get_value_item_keyword( ConvergenceItem_txt_, 2, "forceFieldSigma" ), 3 ) + " kcal/mol/atom, respectively. In step 1 (structure generation), " +
                size_t2string( get_value_item_keyword( ConvergenceItem_txt_, 0, "nbGeneratedStructures" ) ) + " (" +
                size_t2string( get_value_item_keyword( ConvergenceItem_txt_, 1, "nbGeneratedStructures" ) ) + ", " +
                size_t2string( get_value_item_keyword( ConvergenceItem_txt_, 2, "nbGeneratedStructures" ) ) +
                ") crystal structures were generated within a " +
                double2string_2( get_value_item_keyword( ConvergenceItem_txt_, 0, "actualGenerationWindowExtension" ), 1 ) + " (" +
                double2string_2( get_value_item_keyword( ConvergenceItem_txt_, 1, "actualGenerationWindowExtension" ), 1 ) + ", " +
                double2string_2( get_value_item_keyword( ConvergenceItem_txt_, 2, "actualGenerationWindowExtension" ), 1 ) +
                ") sigma window at a level of completeness of " +
                double2string_2( 100.0 * get_value_item_keyword( ConvergenceItem_txt_, 0, "actualGenerationConvergence" ), 1 ) + "% (" +
                double2string_2( 100.0 * get_value_item_keyword( ConvergenceItem_txt_, 1, "actualGenerationConvergence" ), 1 ) + "%, " +
                double2string_2( 100.0 * get_value_item_keyword( ConvergenceItem_txt_, 2, "actualGenerationConvergence" ), 1 ) +
                "%). The accuracy of the energy calculations in step 2 was found to be " +
                double2string_2( get_value_item_keyword( ConvergenceItem_txt_, 0, "freeEnergySigma" ), 3 ) + " (" +
                double2string_2( get_value_item_keyword( ConvergenceItem_txt_, 1, "freeEnergySigma" ), 3 ) + ", " +
                double2string_2( get_value_item_keyword( ConvergenceItem_txt_, 2, "freeEnergySigma" ), 3 ) +
                ") kcal/mol/atom. In step 2 (DFT-D reranking), " +
                size_t2string( get_value_item_keyword( ConvergenceItem_txt_, 0, "nbRerankedStructures" ) ) + " (" +
                size_t2string( get_value_item_keyword( ConvergenceItem_txt_, 1, "nbRerankedStructures" ) ) + ", " +
                size_t2string( get_value_item_keyword( ConvergenceItem_txt_, 2, "nbRerankedStructures" ) ) +
                ") structures were fully lattice-energy minimized at DFT-D level. In step 3 (final lattice-energy optimization) " +
                size_t2string( get_value_item_keyword( ConvergenceItem_txt_, 0, "nbFreeEnergyCorrectedStructures" ) ) + " (" +
                size_t2string( get_value_item_keyword( ConvergenceItem_txt_, 1, "nbFreeEnergyCorrectedStructures" ) ) + ", " +
                size_t2string( get_value_item_keyword( ConvergenceItem_txt_, 2, "nbFreeEnergyCorrectedStructures" ) ) +
                ") structures were processed. Values in brackets refer to the two sets of space groups chosen for Z'=2." );
            }
            else // item_0_found, item_1_found, ! item_2_found
            {
                TextFileWriter some_numbers_file( FileName( directory_, "SomeNumbers", "txt" ) );
                some_numbers_file.write_line( "The accuracy of the TMFF measured independently for Z'=1 and Z'=2 turned out to be " +
                double2string_2( get_value_item_keyword( ConvergenceItem_txt_, 0, "forceFieldSigma" ), 3 ) + " kcal/mol/atom and " +
                double2string_2( get_value_item_keyword( ConvergenceItem_txt_, 1, "forceFieldSigma" ), 3 ) + " kcal/mol/atom, respectively. In step 1 (structure generation), " +
                size_t2string( get_value_item_keyword( ConvergenceItem_txt_, 0, "nbGeneratedStructures" ) ) + " (" +
                size_t2string( get_value_item_keyword( ConvergenceItem_txt_, 1, "nbGeneratedStructures" ) ) +
                ") crystal structures were generated within a " +
                double2string_2( get_value_item_keyword( ConvergenceItem_txt_, 0, "actualGenerationWindowExtension" ), 1 ) + " (" +
                double2string_2( get_value_item_keyword( ConvergenceItem_txt_, 1, "actualGenerationWindowExtension" ), 1 ) +
                ") sigma window at a level of completeness of " +
                double2string_2( 100.0 * get_value_item_keyword( ConvergenceItem_txt_, 0, "actualGenerationConvergence" ), 1 ) + "% (" +
                double2string_2( 100.0 * get_value_item_keyword( ConvergenceItem_txt_, 1, "actualGenerationConvergence" ), 1 ) +
                "%). The accuracy of the energy calculations in step 2 was found to be " +
                double2string_2( get_value_item_keyword( ConvergenceItem_txt_, 0, "freeEnergySigma" ), 3 ) + " (" +
                double2string_2( get_value_item_keyword( ConvergenceItem_txt_, 1, "freeEnergySigma" ), 3 ) +
                ") kcal/mol/atom. In step 2 (DFT-D reranking), " +
                size_t2string( get_value_item_keyword( ConvergenceItem_txt_, 0, "nbRerankedStructures" ) ) + " (" +
                size_t2string( get_value_item_keyword( ConvergenceItem_txt_, 1, "nbRerankedStructures" ) ) +
                ") structures were fully lattice-energy minimized at DFT-D level. In step 3 (final lattice-energy optimization) " +
                size_t2string( get_value_item_keyword( ConvergenceItem_txt_, 0, "nbFreeEnergyCorrectedStructures" ) ) + " (" +
                size_t2string( get_value_item_keyword( ConvergenceItem_txt_, 1, "nbFreeEnergyCorrectedStructures" ) ) +
                ") structures were processed. Values in brackets refer to the space groups chosen for Z'=2." );
            }
        }
        else // item_0_found, ! item_1_found
        {
            TextFileWriter some_numbers_file( FileName( directory_, "SomeNumbers", "txt" ) );
            some_numbers_file.write_line( "The accuracy of the TMFF measured for Z'=1 turned out to be " +
            double2string_2( get_value_item_keyword( ConvergenceItem_txt_, 0, "forceFieldSigma" ), 3 ) + " kcal/mol/atom. In step 1 (structure generation), " +
            size_t2string( get_value_item_keyword( ConvergenceItem_txt_, 0, "nbGeneratedStructures" ) ) + " crystal structures were generated within a " +
            double2string_2( get_value_item_keyword( ConvergenceItem_txt_, 0, "actualGenerationWindowExtension" ), 1 ) + " sigma window at a level of completeness of " +
            double2string_2( 100.0 * get_value_item_keyword( ConvergenceItem_txt_, 0, "actualGenerationConvergence" ), 1 ) + "%. The accuracy of the energy calculations in step 2 was found to be " +
            double2string_2( get_value_item_keyword( ConvergenceItem_txt_, 0, "freeEnergySigma" ), 3 ) + " kcal/mol/atom. In step 2 (DFT-D reranking), " +
            size_t2string( get_value_item_keyword( ConvergenceItem_txt_, 0, "nbRerankedStructures" ) ) + " structures were fully lattice-energy minimized at DFT-D level. In step 3 (final lattice-energy optimization) " +
            size_t2string( get_value_item_keyword( ConvergenceItem_txt_, 0, "nbFreeEnergyCorrectedStructures" ) ) + " structures were processed." );
        }
    }
}

// ********************************************************************************

void D32GRACE()
{
    TextFileReader text_file_reader( FileName( "C:\\GD\\D3\\dftd3.3.1.0\\pars.f" ) );
    text_file_reader.set_skip_empty_lines( true );
    text_file_reader.set_allow_single_quotes( true );
    TextFileWriter text_file_writer( FileName( "C:\\GD\\D3\\C6ReferenceValues_3.txt" ) );
    Splitter splitter( "," );
    splitter.set_merge_delimiters( false );
    std::vector< std::string > words;
    std::string line;
    while ( text_file_reader.get_next_line( words ) )
    {
        line = text_file_reader.get_line();
        if ( line.length() < 55 )
            continue;
        // The first character should be a "."
        line = strip( line );
        if ( line.substr( 0, 1 ) != "." )
        {
            std::cout << "Oops 1..." << std::endl;
            std::cout << line << std::endl;
            continue;
        }
        line = line.substr( 1 );
        // The first character should now be " " or ",", either must be discarded.
        if ( ( line.substr( 0, 1 ) != " " ) && ( line.substr( 0, 1 ) != "," ) )
        {
            std::cout << "Oops 2..." << std::endl;
            std::cout << line << std::endl;
            continue;
        }
        line = line.substr( 1 );
        words = splitter.split( line );
        if ( words.size() != 5 )
        {
            std::cout << "Words.size() != 5" << std::endl;
            std::cout << line << std::endl;
            continue;
        }
        std::vector< double > values( 5, 0.0 );
        for ( size_t i( 0 ); i != words.size(); ++i )
        {
            words[i] = replace( words[i], "D", "E" );
            values[i] = string2double( words[i] );
        }
        while ( values[1] > 100.0 )
            values[1] -= 100.0;
        while ( values[2] > 100.0 )
            values[2] -= 100.0;
        double conversion_factor = 13.779294968; // Ha * a0^6 to Kcal/mol * A ^ 6
        text_file_writer.write(        size_t2string( round_to_int( values[1] ), 3, ' ' ) );
        text_file_writer.write( "  " + size_t2string( round_to_int( values[2] ), 3, ' ' ) );
        text_file_writer.write( "  " + double2string( values[3], 4 ) );
        text_file_writer.write( "  " + double2string( values[4], 4 ) );
        text_file_writer.write( "  " + double2string( conversion_factor * values[0], 8, 14, ' ' ) );
        text_file_writer.write_line();
    }
}

// ********************************************************************************

