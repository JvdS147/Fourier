#ifndef AMS_H
#define AMS_H

#include "AMS_TableFile.h"
#include "FileList.h"
#include "TextFileReader_2.h"

#include <string>

int AMS_main( int argc, char** argv );

double Celsius_to_Kelvin( double const temperature );

double bar_to_kiloPascal( double const pressure );

double kJ_to_kcal( double const energy );

// Gas constant in kcal / K / mol
double R();

void process_RESULTS_file( const FileName & input_file_name, const size_t natoms );
    
void analyse_CSP_results( const FileList & file_list );

void D32GRACE();

class ReportGenerator
{
public:

    // The FileList.txt is generated from table.txt. table.txt must include a directory.
    // The natoms_per_molecule is a fraction to allow for "hemihydrates" to really be API + 1/2 water
    // @@ We need a flag bool is_hemihydrate
    ReportGenerator( const FileName & table_txt_file_name, const Fraction natoms_per_molecule );

    void generate_report() const;
    
    void generate_some_numbers() const;

private:

    Fraction natoms_per_molecule_;
    std::string directory_; // Always ends in backslash
    AMS_TableFile table_file_;
    bool ConvergenceInformation_file_found_;
    TextFileReader_2 ConvergenceItem_txt_;
    
    FileList file_list_;
};

#endif // AMS_H

