#include "RunTests.h"
#include "TextFileReader.h"
#include "TextFileReader_2.h"
#include "TextFileWriter.h"
#include "FileName.h"
#include "FileList.h"
#include "ReadCif.h"
#include "Sort.h"
#include "CrystalStructure.h"
#include "SimilarityAnalysis.h"
#include "CorrelationMatrix.h"
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
    try // Calculate similarity matrix.
    {
        MACRO_ONE_FILELISTNAME_AS_ARGUMENT
        CorrelationMatrix similarity_matrix = calculate_correlation_matrix( file_list );
        similarity_matrix.save( FileName( file_list_file_name.directory(), append_to_file_name( file_list_file_name, "_SimilarityMatrix" ).file_name(), "txt" ) );
    MACRO_END_GAME
}
