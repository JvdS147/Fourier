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

#include "PowderPattern.h"
#include "FileName.h"
#include "MathsFunctions.h"
#include "RunningAverageAndESD.h"
#include "TextFileReader.h"
#include "TextFileReader_2.h"
#include "TextFileWriter.h"
#include "Utilities.h"
#include "Vector3D.h" // Should not have been necessary

#include <fstream>
#include <stdexcept>

#include <iostream> // for debugging

// ********************************************************************************

PowderPattern::PowderPattern() : wavelength_(1.54056)
{
}

// ********************************************************************************

PowderPattern::PowderPattern( const Angle two_theta_start, const Angle two_theta_end, const Angle two_theta_step ):
wavelength_(1.54056),
noise_is_available_(false)
{
    size_t npoints = round_to_int( ( (two_theta_end-two_theta_start) / two_theta_step ) ) + 1;
    two_theta_values_.reserve( npoints );
    for ( size_t i( 0 ); i != npoints; ++i )
        two_theta_values_.push_back( ( i * two_theta_step ) + two_theta_start );
    intensities_ = std::vector<double>( npoints, 0.0 );
    estimated_standard_deviations_ = std::vector<double>( npoints, 0.0 );
}

// ********************************************************************************

PowderPattern::PowderPattern( const FileName & file_name ):
wavelength_(1.54056),
noise_is_available_(false)
{
    read_xye( file_name );
}

// ********************************************************************************

void PowderPattern::reserve( const size_t nvalues )
{
    two_theta_values_.reserve( nvalues );
    intensities_.reserve( nvalues );
    estimated_standard_deviations_.reserve( nvalues );
}

// ********************************************************************************

void PowderPattern::push_back( const Angle two_theta, const double intensity )
{
    two_theta_values_.push_back( two_theta );
    intensities_.push_back( intensity );
    estimated_standard_deviations_.push_back( std::max( sqrt( intensity ), intensity / 100.0 ) );
}

// ********************************************************************************

void PowderPattern::push_back( const Angle two_theta, const double intensity, const double estimated_standard_deviation )
{
    two_theta_values_.push_back( two_theta );
    intensities_.push_back( intensity );
    estimated_standard_deviations_.push_back( estimated_standard_deviation );
}

// ********************************************************************************

void PowderPattern::rebin( const size_t bin_size )
{
    if ( bin_size < 2 )
        return;
    if ( empty() )
        return;
    PowderPattern result;
    RunningAverageAndESD< Angle > current_two_theta( two_theta( 0 ) );
    RunningAverageAndESD< double > current_intensity( intensity( 0 ) );
    double current_ESD( square( estimated_standard_deviation( 0 ) ) );
    for ( size_t i( 1 ); i != size() ; ++ i )
    {
        current_two_theta.add_value( two_theta( i ) );
        current_intensity.add_value( intensity( i ) );
        current_ESD += square( estimated_standard_deviation( i ) );
        if ( ( ( i+1 ) % bin_size ) == 0 )
        {
            result.push_back( current_two_theta.average(), current_intensity.average(), std::sqrt( current_ESD ) );
            current_two_theta.clear();
            current_intensity.clear();
            current_ESD = 0.0;
        }
    }
    if ( current_two_theta.nvalues() != 0 )
        result.push_back( current_two_theta.average(), current_intensity.average(), std::sqrt( current_ESD ) );
    *this = result;
}

// ********************************************************************************

size_t PowderPattern::find_two_theta( const Angle two_theta_value ) const
{
    if ( empty() )
        throw std::runtime_error( "PowderPattern::find_two_theta(): no data points." );
    if ( two_theta_value < two_theta_start() )
        return 0;
    if ( two_theta_value > two_theta_end() )
        return size()-1;
    // Initialise with guess based on uniform 2 theta step
    size_t lower_index = ( two_theta_value - two_theta_start() ) / average_two_theta_step();
    while ( ( lower_index != 0 ) && ( two_theta( lower_index ) >= two_theta_value ) )
    {
        --lower_index;
    }
    size_t upper_index = lower_index;
    while ( ( upper_index != size()-1 ) && ( two_theta( upper_index ) <= two_theta_value ) )
    {
        ++upper_index;
    }
    lower_index = upper_index;
    while ( ( lower_index != 0 ) && ( two_theta( lower_index ) >= two_theta_value ) )
    {
        --lower_index;
    }
  //  std::cout << two_theta_value << std::endl;
  //  std::cout << lower_index << std::endl;
  //  std::cout << upper_index << std::endl;
    if ( ( upper_index - lower_index ) != 1 )
        throw std::runtime_error( "PowderPattern::find_two_theta(): programming error." );
    if ( ( two_theta_value - two_theta( lower_index ) ) < ( two_theta( upper_index ) - two_theta_value ) )
        return lower_index;
    else
        return upper_index;
}

// ********************************************************************************

// Multiplies intensities and ESDs by factor
void PowderPattern::scale( const double factor )
{
    for ( size_t i( 0 ); i != size(); ++i )
    {
        intensities_[i] *= factor;
        estimated_standard_deviations_[i] *= factor;
    }
}

// ********************************************************************************

Angle PowderPattern::average_two_theta_step() const
{
    if ( empty() )
        throw std::runtime_error( "PowderPattern::average_two_theta_step(): no data points." );
    if ( size() == 1 )
        throw std::runtime_error( "PowderPattern::average_two_theta_step(): only one data point." );
    return ( ( two_theta( size()-1 ) - two_theta( 0 ) ) / ( size() - 1 ) );
}

// ********************************************************************************

Angle PowderPattern::two_theta_start() const
{
    if ( empty() )
        throw std::runtime_error( "PowderPattern::two_theta_start(): no data points." );
    return two_theta_values_[ 0 ];
}

// ********************************************************************************

Angle PowderPattern::two_theta_end() const
{
    if ( empty() )
        throw std::runtime_error( "PowderPattern::two_theta_end(): no data points." );
    return two_theta_values_[ size()-1 ];
}

// ********************************************************************************

// Uses average_two_theta_step() to add new points. Intensities and ESDs are initialised to 0.0.
void PowderPattern::set_two_theta_start( const Angle two_theta_start ) const
{
}

// ********************************************************************************

// Uses average_two_theta_step() to add new points. Intensities and ESDs are initialised to 0.0.
void PowderPattern::set_two_theta_end( const Angle two_theta_end ) const
{
}

// ********************************************************************************

double PowderPattern::cumulative_intensity() const
{
    return add_doubles( intensities_ );
}

// ********************************************************************************

double PowderPattern::cumulative_noise() const
{
    return add_doubles( noise_ );
}

// ********************************************************************************

double PowderPattern::cumulative_absolute_noise() const
{
    return add_absolute_doubles( noise_ );
}

// ********************************************************************************

double PowderPattern::cumulative_squared_noise() const
{
    return add_squared_doubles( noise_ );
}

// ********************************************************************************

void PowderPattern::read_xye( const FileName & file_name )
{
    *this = PowderPattern();
    TextFileReader text_file_reader( file_name );
    std::vector< std::string > words;
    // The first line could contain the wavelength.
    if ( text_file_reader.get_next_line( words ) )
    {
        if ( words.size() == 1 )
            wavelength_ = string2double( words[0] );
        else
            text_file_reader.push_back_last_line();
    }
    else
        return;
    while ( text_file_reader.get_next_line( words ) )
    {
        if ( ( words.size() < 2 ) || ( words.size() > 3 ) )
            throw std::runtime_error( "PowderPattern::read(): cannot interpret line \"" + text_file_reader.get_line() + "\"" );
        two_theta_values_.push_back( Angle( string2double( words[0] ), Angle::DEGREES ) );
        intensities_.push_back( string2double( words[1] ) );
        if ( words.size() == 2 )
            estimated_standard_deviations_.push_back( sqrt( intensities_[two_theta_values_.size()-1] ) );
        else
            estimated_standard_deviations_.push_back( string2double( words[2] ) );
    }
}

// ********************************************************************************

void PowderPattern::read_xrdml( const FileName & file_name )
{
    *this = PowderPattern();
    TextFileReader_2 text_file_reader( file_name );
    size_t iPos = text_file_reader.find( "<positions axis=\"2Theta\" unit=\"deg\">" );
    if ( iPos == std::string::npos )
        throw std::runtime_error( "PowderPattern::read_xrdml(): 2theta not found." );
    std::string two_theta_start_str = extract_delimited_text( text_file_reader.line( iPos + 1 ), "<startPosition>", "</startPosition>");
    std::string two_theta_end_str   = extract_delimited_text( text_file_reader.line( iPos + 2 ), "<endPosition>", "</endPosition>" );
    Angle two_theta_start = Angle::from_degrees( string2double( two_theta_start_str ) );
    Angle two_theta_end   = Angle::from_degrees( string2double( two_theta_end_str   ) );
    std::vector< std::string > divergence_corrections;
    iPos = text_file_reader.find( "<divergenceCorrections>" );
    if ( iPos != std::string::npos )
         divergence_corrections = split( extract_delimited_text( text_file_reader.line( iPos ), "<divergenceCorrections>", "</divergenceCorrections>" ) );
    std::vector< std::string > counts;
    iPos = text_file_reader.find( "<intensities unit=\"counts\">" );
    if ( iPos != std::string::npos )
         counts = split( extract_delimited_text( text_file_reader.line( iPos ), "<intensities unit=\"counts\">", "</intensities>" ) );
    else
    {
        iPos = text_file_reader.find( "<counts unit=\"counts\">" );
        if ( iPos != std::string::npos )
            counts = split( extract_delimited_text( text_file_reader.line( iPos ), "<counts unit=\"counts\">", "</counts>" ) );
        else
            throw std::runtime_error( "PowderPattern::read_xrdml(): Counts not found." );
    }
    if ( counts.empty() )
        throw std::runtime_error( "PowderPattern::read_xrdml(): no data points." );
    if ( counts.size() == 1 )
        throw std::runtime_error( "PowderPattern::read_xrdml(): only one data point." );
    reserve( counts.size() );
    Angle two_theta_step = ( two_theta_end - two_theta_start ) / ( counts.size() - 1 );
    if ( divergence_corrections.empty() )
    {
        for ( size_t i( 0 ); i != counts.size(); ++i )
            push_back( ( i * two_theta_step ) + two_theta_start, string2double( counts[i] ) );
    }
    else
    {
        if ( divergence_corrections.size() != counts.size() )
            throw std::runtime_error( "PowderPattern::read_xrdml(): number of counts and number of divergence corrections differ." );
        std::cout << "Note that the .xrdml file contains divergence corrections, which will be applied to the counts." << std::endl;
        for ( size_t i( 0 ); i != counts.size(); ++i )
            push_back( ( i * two_theta_step ) + two_theta_start, string2double( divergence_corrections[i] ) * string2double( counts[i] ) );
    }
}

// ********************************************************************************

//BANK       1    3501     350  CONST    23.20    2.30     0.0     0.0         STD
//       1       3       1       0       3       2       4       2       2       2
//       1       6       4       3       6       6       5       3       2       4
//       7      11      10      28      35      69     119     160     272     342
// The start is 0.232, the 2theta step size is 0.0230
void PowderPattern::read_raw( const FileName & file_name )
{
    *this = PowderPattern();
    // This is lab data (is that always true?), the wavelength is fine.
    TextFileReader text_file_reader( file_name );
    text_file_reader.set_skip_empty_lines( false );
    std::vector< std::string > words;
    if ( ! text_file_reader.get_next_line( words ) )
        throw std::runtime_error( "PowderPattern::read_raw(): File is empty." );
    if ( ! text_file_reader.get_next_line( words ) )
        throw std::runtime_error( "PowderPattern::read_raw(): File is empty." );
    if ( words.size() != 10 )
        throw std::runtime_error( "PowderPattern::read_raw(): unexpected format 1." );
    if ( words[0] != "BANK" )
        throw std::runtime_error( "PowderPattern::read_raw(): unexpected format 2." );
    if ( words[4] != "CONST" )
        throw std::runtime_error( "PowderPattern::read_raw(): unexpected format 3." );
    bool read_ESDs( false );
    if ( words[9] == "ESD" )
        read_ESDs = true;
    else if ( words[9] != "STD" )
        throw std::runtime_error( "PowderPattern::read_raw(): unexpected format 4." );
    size_t ndata_points = string2integer( words[2] );
    if ( ndata_points == 0 )
        throw std::runtime_error( "PowderPattern::read_raw(): No data." );
    Angle two_theta_start = Angle::from_degrees( string2double( words[5] ) / 100.0 );
    Angle two_theta_step = Angle::from_degrees( string2double( words[6] ) / 100.0 );
    size_t i( 0 );
    std::string line;
    Splitter splitter;
    splitter.split_by_length( 8 );
    while ( text_file_reader.get_next_line( line ) )
    {
        std::vector< std::string > temp_words = splitter.split( line );
        words.clear();
        bool one_word_was_empty( false );
        for ( size_t j( 0 ); j != temp_words.size(); ++j )
        {
            temp_words[j] = strip( temp_words[j] );
            if ( temp_words[j].empty() )
                one_word_was_empty = true;
            else
            {
                if ( one_word_was_empty )
                    throw std::runtime_error( "PowderPattern::read_raw(): non-empty word after empty word." );
                words.push_back( temp_words[j] );
            }
        }
        if ( read_ESDs )
        {
            if ( is_odd( words.size() ) )
                throw std::runtime_error( "PowderPattern::read_raw(): intensities plus ESDs stored, but number of values is odd." );
            for ( size_t j( 0 ); j != words.size(); j += 2 )
            {
                push_back( ( i * two_theta_step ) + two_theta_start, string2double( words[j] ), string2double( words[j+1] ) );
                ++i;
            }
        }
        else
        {
            for ( size_t j( 0 ); j != words.size(); ++j )
            {
                push_back( ( i * two_theta_step ) + two_theta_start, string2double( words[j] ) );
                ++i;
            }
        }
    }
    if ( i != ndata_points )
        std::cout << "PowderPattern::read_raw(): Warning: the number of data points in the file disagrees with the number in the header." << std::endl;
}

// ********************************************************************************

//08/06/2018 07:47:49  DIF       : t=   600s
//  0.4460 0.0460  1.0 US 1.5418     87.5240  1893
//      3      2      4      3      4      3      6      4
//      3      2      6      9      8      6      7     15
//      9     11     21     23     24     64    188    323
// The start is 0.4460, the step is 0.0460 and the end is 87.5240-0.0460. There are 1893 points
// The 2theta end value is wrong, because the last data point is a dummy data point with a value of -1, and it is *not* counted towards the number of data points,
// but it *is* counted towards the end 2theta value.
void PowderPattern::read_mdi( const FileName & file_name )
{
    *this = PowderPattern();
    // This is lab data (is that always true?), the wavelength is fine.
    TextFileReader text_file_reader( file_name );
    text_file_reader.set_skip_empty_lines( false );
    std::vector< std::string > words;
    if ( ( ! text_file_reader.get_next_line( words ) ) || ( words.size() == 0 ) )
        throw std::runtime_error( "PowderPattern::read_mdi(): First line is empty." );
    if ( ! text_file_reader.get_next_line( words ) )
        throw std::runtime_error( "PowderPattern::read_mdi(): File is empty." );
    if ( words.size() != 7 )
        throw std::runtime_error( "PowderPattern::read_mdi(): unexpected format." );
    if ( words[3] != "US" )
        throw std::runtime_error( "PowderPattern::read_mdi(): unexpected format." );
    size_t ndata_points = string2integer( words[6] );
    if ( ndata_points == 0 )
        throw std::runtime_error( "PowderPattern::read_mdi(): No data." );
    Angle two_theta_start = Angle::from_degrees( string2double( words[0] ) );
    Angle two_theta_step = Angle::from_degrees( string2double( words[1] ) );
    Angle two_theta_end = Angle::from_degrees( string2double( words[5] ) );
    size_t i( 0 );
    bool we_are_done( false );
    while ( text_file_reader.get_next_line( words ) )
    {
        for ( size_t j( 0 ); j != words.size(); ++j )
        {
            // There is one dummy data point with value -1 after the last valid data point
            if ( ( i == ndata_points ) && ( words[j] == "-1" ) && ( j == ( words.size() - 1 ) ) )
            {
                we_are_done = true;
            }
            else
            {
                if ( we_are_done )
                    throw std::runtime_error( "PowderPattern::read_mdi(): data found after last data point." );
                push_back( ( i * two_theta_step ) + two_theta_start, string2double( words[j] ) );
                ++i;
            }
        }
    }
    if ( i != ndata_points )
        std::cout << "PowderPattern::read_mdi(): Warning: the number of data points in the file (" + size_t2string( i ) + ") disagrees with the number in the header (" + size_t2string( ndata_points ) + ")." << std::endl;
    else
        std::cout << "PowderPattern::read_mdi(): the number of data points in the file (" + size_t2string( i ) + ") agrees with the number in the header (" + size_t2string( ndata_points ) + ")." << std::endl;
    std::cout << "PowderPattern::read_mdi(): two_theta_end as calculated     = " << ( i * two_theta_step ) + two_theta_start << std::endl;
    std::cout << "PowderPattern::read_mdi(): two_theta_end as read from file = " << two_theta_end << std::endl;
}

// ********************************************************************************

//      <SubScans>
//        <SubScanInfo Steps="1051" MeasuredSteps="1051" StartStepNo="0" MeasuredTimePerStep="260" PlannedTimePerStep="2" />
//      </SubScans>
//      <Datum>260,1,2,1,662</Datum>
//      <Datum>260,1,2.0409,1.0205,599</Datum>
//      <Datum>260,1,2.0819,1.0409,603</Datum>
void PowderPattern::read_brml( const FileName & file_name )
{
    *this = PowderPattern();
    // This is lab data (is that always true?), the wavelength is fine.
    TextFileReader text_file_reader( file_name );
    text_file_reader.set_skip_empty_lines( false );
    Splitter splitter( "," );
    std::string line;
    while ( text_file_reader.get_next_line( line ) )
    {
        line = strip( line );
        line = to_upper( line );
        if ( line.substr( 0, 7 ) != "<DATUM>" )
            continue;
        line = extract_delimited_text( line, "<DATUM>", "</DATUM>" );
        std::vector< std::string > words = splitter.split( line );
        push_back( Angle::from_degrees( string2double( words[2] ) ), string2double( words[4] ) );
    }
}

// ********************************************************************************

void PowderPattern::read_txt( const FileName & file_name )
{
    *this = PowderPattern();
    TextFileReader text_file_reader( file_name );
    text_file_reader.set_skip_empty_lines( true );
    std::vector< std::string > words;
    size_t stage( 1 );
    Angle two_theta_start;
    Angle two_theta_step;
    Angle two_theta_end;
    size_t i( 0 );
    while ( text_file_reader.get_next_line( words ) )
    {
        if ( words[0] == "DATE" ||
             words[0] == "ACQTIME" ||
             words[0] == "VOLTAGE" ||
             words[0] == "CURRENT" ||
             words[0] == "WAVELENGTH" ||
             words[0] == "COMMENT1" ||
             words[0] == "COMMENT2" ||
             words[0] == "COMMENT3" )
        {
            if ( stage == 1 )
                continue;
            else
                throw std::runtime_error( "PowderPattern::read_txt(): keyword out of place." );
        }
        if ( words.size() == 3 )
        {
            if ( stage != 1 )
                throw std::runtime_error( "PowderPattern::read_txt(): keyword out of place." );
            two_theta_start = Angle::from_degrees( string2double( words[0] ) );
            two_theta_step  = Angle::from_degrees( string2double( words[1] ) );
            two_theta_end   = Angle::from_degrees( string2double( words[2] ) );
            stage = 3;
            continue;
        }
        if ( stage != 3 )
            throw std::runtime_error( "PowderPattern::read_txt(): keyword out of place." );
        if ( words.size() != 1 )
            throw std::runtime_error( "PowderPattern::read_txt(): keyword out of place." );
        push_back( ( i * two_theta_step ) + two_theta_start, string2double( words[0] ) );
        ++i;
    }
}

// ********************************************************************************

void PowderPattern::save_xye( const FileName & file_name, const bool include_wave_length ) const
{
    TextFileWriter text_file_writer( file_name );
    if ( include_wave_length )
        text_file_writer.write_line( double2string( wavelength_ ) );
    for ( size_t i( 0 ); i != size(); ++i )
        text_file_writer.write_line( double2string( two_theta_values_[i].value_in_degrees(), 5 ) + "  " + double2string( intensities_[i] ) + "  " + double2string( estimated_standard_deviations_[i] ) );
}

// ********************************************************************************

void PowderPattern::generate_code( const bool include_estimated_standard_deviation ) const
{
    std::cout << "    PowderPattern powder_pattern;" << std::endl;
    if ( include_estimated_standard_deviation  )
    {
        for ( size_t i( 0 ); i != size(); ++i )
            std::cout << "    powder_pattern.push_back( Angle::from_degrees( " << two_theta( i ) << " ), " << intensity( i ) << ", " << estimated_standard_deviation( i ) << " );" << std::endl;
    }
    else
    {
        for ( size_t i( 0 ); i != size(); ++i )
            std::cout << "    powder_pattern.push_back( Angle::from_degrees( " << two_theta( i ) << " ), " << intensity( i ) << " );" << std::endl;
    }
}

// ********************************************************************************

// We assume a uniform step size.
PowderPattern & PowderPattern::operator+=( const PowderPattern & rhs )
{
    if ( ! same_range( *this, rhs ) )
        throw std::runtime_error( "PowderPattern::operator+=( const PowderPattern & ): ranges not same." );
    for ( size_t i( 0 ); i != size(); ++i )
        intensities_[i] += rhs.intensities_[i];
    return *this;
}

// ********************************************************************************

// We assume a uniform step size.
PowderPattern & PowderPattern::operator-=( const PowderPattern & rhs )
{
    if ( ! same_range( *this, rhs ) )
        throw std::runtime_error( "PowderPattern::operator-=( const PowderPattern & ): ranges not same." );
    for ( size_t i( 0 ); i != size(); ++i )
        intensities_[i] -= rhs.intensities_[i];
    return *this;
}

// ********************************************************************************

void PowderPattern::normalise_highest_peak( const double highest_peak )
{
    // Find the highest intensity
    double max_intensity( 0.0 );
    for ( size_t i( 0 ); i != size(); ++i )
    {
        if ( max_intensity < intensities_[i] )
            max_intensity = intensities_[i];
    }
    if ( nearly_zero( max_intensity ) )
        throw std::runtime_error( "PowderPattern::normalise_highest_peak(): highest peak is 0.0." );
    scale( highest_peak / max_intensity );
}

// ********************************************************************************

// Normalises the total signal = area under the pattern = cumulative_intensity() .
void PowderPattern::normalise_total_signal( const double total_signal )
{
    double current_total_signal = cumulative_intensity();
    if ( nearly_zero( current_total_signal ) )
        throw std::runtime_error( "PowderPattern::normalise_total_signal(): total signal is 0.0." );
    // Scale to total_signal
    scale( total_signal / current_total_signal );
}

// ********************************************************************************

void PowderPattern::correct_zero_point_error( const Angle two_theta_value )
{
    for ( size_t i( 0 ); i != size(); ++i )
        two_theta_values_[i] -= two_theta_value;
}

// ********************************************************************************

void PowderPattern::recalculate_estimated_standard_deviations()
{
    for ( size_t i( 0 ); i != size(); ++i )
    {
        if ( intensities_[i] < 20.0 )
        {
            estimated_standard_deviations_[i] = 4.4;
        }
        else if ( intensities_[i] > 10000.0 )
        {
            estimated_standard_deviations_[i] = intensities_[i] / 100.0;
        }
        else
        {
            estimated_standard_deviations_[i] = sqrt( intensities_[i] );
        }
    }
}

// ********************************************************************************

// I(fixed slit) = I(variable slit) / sin(theta).
void PowderPattern::convert_to_fixed_slit()
{
    for ( size_t i( 0 ); i != size(); ++i )
    {
        intensities_[i] = intensities_[i] / ( two_theta_values_[i] / 2.0 ).sine(); // @@ We should check for divide by zero
        estimated_standard_deviations_[i] = estimated_standard_deviations_[i] / ( two_theta_values_[i] / 2.0 ).sine();
    }
}

// ********************************************************************************

// I(fixed slit) = I(variable slit) / sin(theta).
void PowderPattern::convert_to_variable_slit()
{
    for ( size_t i( 0 ); i != size(); ++i )
    {
        intensities_[i] = intensities_[i] * ( two_theta_values_[i] / 2.0 ).sine();
        estimated_standard_deviations_[i] = estimated_standard_deviations_[i] * ( two_theta_values_[i] / 2.0 ).sine();
    }
}

// ********************************************************************************

void PowderPattern::add_constant_background( const double background )
{
    for ( size_t i( 0 ); i != size(); ++i )
        intensities_[i] += background;
}

// ********************************************************************************

void PowderPattern::add_Poisson_noise()
{
    for ( size_t i( 0 ); i != size(); ++i )
    {
        double old_intensity = intensities_[i];
        intensities_[i] = Poisson_distribution( round_to_int( old_intensity ) );
        noise_.push_back( intensities_[i] - old_intensity );
    }
    noise_is_available_ = true;
}

// ********************************************************************************

void PowderPattern::make_counts_integer()
{
    for ( size_t i( 0 ); i != size(); ++i )
        intensities_[i] = round_to_int( intensities_[i] );
}

// ********************************************************************************

// Check that the two patterns have the same range and 2theta step
bool same_range( const PowderPattern & lhs, const PowderPattern & rhs )
{
    if ( lhs.size() != rhs.size() )
        return false;
    if ( lhs.size() == 0 )
        return true;
    return ( nearly_equal( lhs.two_theta( 0 )           , rhs.two_theta( 0 ) ) &&
             nearly_equal( lhs.two_theta( lhs.size()-1 ), rhs.two_theta( rhs.size()-1 ) ) );
}

// ********************************************************************************

double weighted_cross_correlation( const PowderPattern & lhs, const PowderPattern & rhs, Angle l )
{
    int m = round_to_int( l / lhs.average_two_theta_step() );
    if ( m == 0 )
        m = 1;
    double result( 0.0 );
    for ( int i( 0 ); i != lhs.size(); ++i )
    {
        for ( int j( -m + 1 ); j != m; ++j )
        {
            if ( ( ( i + j ) >= 0 ) && ( ( i + j ) < lhs.size() ) )
            {
                if ( (true) )
                    result += (1.0-std::abs(j)/m) * lhs.intensity( i ) * rhs.intensity( i + j );
                else
                    result += (1.0-std::abs(j)/m) * ( lhs.intensity( i ) / lhs.estimated_standard_deviation( i ) ) * ( rhs.intensity( i + j ) / rhs.estimated_standard_deviation( i + j ) );
            }
        }
    }
    return result;
}

// ********************************************************************************

double normalised_weighted_cross_correlation( const PowderPattern & lhs, const PowderPattern & rhs, Angle l )
{
    if ( ! same_range( lhs, rhs ) )
        throw std::runtime_error( "normalised_weighted_cross_correlation( const PowderPattern &, const PowderPattern & ): ranges not same" );
    return weighted_cross_correlation( lhs, rhs, l ) / sqrt( weighted_cross_correlation( lhs, lhs, l ) * weighted_cross_correlation( rhs, rhs, l ) );
}

// ********************************************************************************

double Rwp( const PowderPattern & lhs, const PowderPattern & rhs )
{
    double numerator( 0.0 );
    double denominator( 0.0 );
    for ( size_t i( 0 ); i != lhs.size(); ++i )
    {
        numerator   += square( lhs.intensity( i ) - rhs.intensity( i ) ) / square( lhs.estimated_standard_deviation( i ) );
        denominator += square( lhs.intensity( i ) ) / square( lhs.estimated_standard_deviation( i ) );
    }
    return sqrt( numerator / denominator );
}

// ********************************************************************************

PowderPattern calculate_Brueckner_background( const PowderPattern & powder_pattern,
                                              const size_t niterations,
                                              const size_t window,
                                              const bool apply_smoothing,
                                              const size_t smoothing_window )
{
    if ( powder_pattern.size() == 0 )
        return powder_pattern;
    PowderPattern pp_new( powder_pattern );
    size_t size( powder_pattern.size() );
    if ( apply_smoothing )
    {
        for ( size_t i( 0 ); i < size; ++i )
        {
            double new_value = powder_pattern.intensity( i );
            for ( size_t j( 1 ); j <= smoothing_window; ++j )
            {
                new_value += powder_pattern.intensity( std::max( int(i)-int(j), int(0) ) );
                new_value += powder_pattern.intensity( std::min( i+j, size-1 ) );
            }
            pp_new.set_intensity( i, new_value / ( 2.0 * smoothing_window + 1.0 ) );
        }
    }
    if ( true )
    {
        RunningAverageAndESD<double> I_average;
        double I_minimum = pp_new.intensity( 0 );
        for ( size_t i( 0 ); i < size; ++i )
        {
            if ( pp_new.intensity( i ) < I_minimum )
                I_minimum = pp_new.intensity( i );
            I_average.add_value( pp_new.intensity( i ) );
        }
        for ( size_t i( 0 ); i < size; ++i )
        {
            if ( pp_new.intensity( i ) > ( I_average.average() + 2.0 * ( I_average.average() - I_minimum ) ) )
                pp_new.set_intensity( i, ( I_average.average() + 2.0 * ( I_average.average() - I_minimum ) ) );
        }
    }
    for ( size_t iter( 0 ); iter < niterations; ++iter )
    {
        PowderPattern pp_old = pp_new;
        for ( size_t i( 0 ); i < size; ++i )
        {
            double average_value( 0.0 );
            for ( size_t j( 1 ); j <= window; ++j )
            {
                average_value += pp_old.intensity( std::max( int(i)-int(j), int(0) ) );
                average_value += pp_old.intensity( std::min( i+j, size-1 ) );
            }
            average_value /= 2.0 * window;
            pp_new.set_intensity( i, average_value );
        }
        for ( size_t i( 0 ); i < size; ++i )
            pp_new.set_intensity( i, std::min( pp_old.intensity( i ), pp_new.intensity( i ) ) );
    }
    return pp_new;
}

// ********************************************************************************

PowderPattern add_powder_patterns( const std::vector< PowderPattern > & powder_patterns, const std::vector< double > & noscp2ts )
{
    if ( powder_patterns.empty() )
        throw std::runtime_error( "add_powder_patterns(): no powder patterns provided." );
    if ( powder_patterns.size() != noscp2ts.size() )
        throw std::runtime_error( "add_powder_patterns(): powder_patterns and noscp2ts not the same size." );
    // Check that they have the same wavelength and average_two_theta_step
    for ( size_t i( 1 ); i != powder_patterns.size(); ++i )
    {
        if ( ! nearly_equal( powder_patterns[0].wavelength(), powder_patterns[i].wavelength() ) )
            throw std::runtime_error( "add_powder_patterns(): wavelengths not the same." );
        if ( ! nearly_equal( powder_patterns[0].average_two_theta_step(), powder_patterns[i].average_two_theta_step() ) )
            throw std::runtime_error( "add_powder_patterns(): average_two_theta_step not the same." );
    }
    // Find the smallest and largest 2theta values
    Angle two_theta_min = powder_patterns[0].two_theta_start();
    Angle two_theta_max = powder_patterns[0].two_theta_end();
    for ( size_t i( 1 ); i != powder_patterns.size(); ++i )
    {
        if ( powder_patterns[i].two_theta_start() < two_theta_min )
            two_theta_min = powder_patterns[i].two_theta_start();
        if ( powder_patterns[i].two_theta_end() < two_theta_max )
            two_theta_max = powder_patterns[i].two_theta_end();
    }
    Angle two_theta_step = powder_patterns[0].average_two_theta_step();
    PowderPattern result( two_theta_min, two_theta_max, two_theta_step );
    for ( size_t j( 0 ); j != result.size(); ++j )
    {
        bool at_least_one_contribution( false );
        double sum_of_intensities( 0.0 );
        double sum_of_noscp2ts( 0.0 );
        for ( size_t i( 0 ); i != powder_patterns.size(); ++i )
        {
            size_t index = powder_patterns[i].find_two_theta( result.two_theta( j ) );
            if ( absolute( powder_patterns[i].two_theta( index ) - result.two_theta( j ) ) < ( two_theta_step / 2.0 ) )
            {
                at_least_one_contribution = true;
                sum_of_intensities += powder_patterns[i].intensity( index );
                sum_of_noscp2ts += noscp2ts[i];
            }
        }
        if ( at_least_one_contribution )
        {
            result.set_intensity( j, sum_of_intensities / sum_of_noscp2ts );
            result.set_estimated_standard_deviation( j, std::max( sqrt( sum_of_intensities ), sum_of_intensities / 100.0 ) / sum_of_noscp2ts );
        }
        else
            std::cout << "add_powder_patterns(): Warning, no contribution." << std::endl;
    }
    return result;
}

// ********************************************************************************

// Should not be necessary. Introduced to manipulate data from a tool that extracted a powder pattern from a bitmap picture.
void PowderPattern::sort_two_theta()
{
    if ( size() < 2 )
        return;
    bool changed( true );
    while ( changed )
    {
        changed = false;
        for ( size_t i( size()-1 ); i != 0; --i )
        {
            if ( two_theta_values_[i] < two_theta_values_[i-1] )
            {
                // @@ The following code strongly suggests that we should have used a struct to hold each triplet of values...
                std::swap( two_theta_values_[i], two_theta_values_[i-1] );
                std::swap( intensities_[i], intensities_[i-1] );
                std::swap( estimated_standard_deviations_[i], estimated_standard_deviations_[i-1] );
                changed = true;
            }
        }
    }
}

// ********************************************************************************

// Should not be necessary. Introduced to manipulate data from a tool that extracted a powder pattern from a bitmap picture.
void PowderPattern::average_if_two_theta_equal()
{
    std::vector< Angle > new_two_theta_values;
    std::vector< double > new_intensities;
    std::vector< double > new_estimated_standard_deviations;
    size_t iPos1 = 0;
    while ( iPos1 < size() )
    {
        Angle sum_two_theta_values( two_theta_values_[iPos1] );
        double sum_intensity( intensities_[iPos1] );
        double sum_estimated_standard_deviation_2( square( estimated_standard_deviations_[iPos1] ) );
        size_t iPos2 = iPos1 + 1;
        while ( ( iPos2 < size() ) && nearly_equal( two_theta_values_[iPos1], two_theta_values_[iPos2] ) )
        {
            sum_two_theta_values += two_theta_values_[iPos2];
            sum_intensity += intensities_[iPos2];
            sum_estimated_standard_deviation_2 += square( estimated_standard_deviations_[iPos2] );
            ++iPos2;
        }
        new_two_theta_values.push_back( sum_two_theta_values / ( iPos2 - iPos1 ) );
        new_intensities.push_back( sum_intensity / ( iPos2 - iPos1 ) );
        new_estimated_standard_deviations.push_back( sqrt( sum_estimated_standard_deviation_2 ) );
        iPos1 = iPos2;
    }
    two_theta_values_ = new_two_theta_values;
    intensities_ = new_intensities;
    estimated_standard_deviations_ = new_estimated_standard_deviations;
}

// ********************************************************************************

