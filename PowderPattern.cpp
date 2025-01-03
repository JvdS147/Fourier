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

#include "PowderPattern.h"
#include "FileName.h"
#include "MathsFunctions.h"
#include "RandomNumberGenerator.h"
#include "RunningAverageAndESD.h"
#include "StringFunctions.h"
#include "TextFileReader.h"
#include "TextFileReader_2.h"
#include "TextFileWriter.h"
#include "Utilities.h"
#include "Vector3D.h" // Should not have been necessary

#include <cmath>
#include <fstream>
#include <stdexcept>

#include <iostream> // for debugging

// ********************************************************************************

PowderPattern::PowderPattern()
{
}

// ********************************************************************************

PowderPattern::PowderPattern( const Angle two_theta_start, const Angle two_theta_end, const Angle two_theta_step )
{
    size_t npoints = round_to_int( ( ( two_theta_end - two_theta_start ) / two_theta_step ) ) + 1;
    if ( ( ( ( two_theta_start + (npoints-1)*two_theta_step ) - two_theta_end ) / two_theta_step ) > 0.1 )
        std::cout << "PowderPattern::PowderPattern(): Warning: start and end not commensurate with step." << std::endl;
    two_theta_values_.reserve( npoints );
    for ( size_t i( 0 ); i != npoints; ++i )
        two_theta_values_.push_back( ( i * two_theta_step ) + two_theta_start );
    intensities_ = std::vector< double >( npoints, 0.0 );
    estimated_standard_deviations_ = std::vector< double >( npoints, 0.0 );
}

// ********************************************************************************

PowderPattern::PowderPattern( const FileName & file_name )
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
    if ( intensity < 20.0 )
        estimated_standard_deviations_.push_back(  4.4 );
    else if ( intensity > 10000.0 )
        estimated_standard_deviations_.push_back( intensity / 100.0 );
    else
        estimated_standard_deviations_.push_back( sqrt( intensity ) );
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
    // Initialise with guess based on uniform 2 theta step.
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
    for ( size_t i( lower_index ); i != upper_index+1; ++i )
    {
        if ( nearly_equal( two_theta( i ), two_theta_value ) )
            return i;
    }
    if ( ( upper_index - lower_index ) != 1 )
        throw std::runtime_error( "PowderPattern::find_two_theta(): programming error." );
    if ( ( two_theta_value - two_theta( lower_index ) ) < ( two_theta( upper_index ) - two_theta_value ) )
        return lower_index;
    else
        return upper_index;
}

// ********************************************************************************

// Multiplies intensities and ESDs by factor.
void PowderPattern::scale( const double factor )
{
    for ( size_t i( 0 ); i != size(); ++i )
    {
        intensities_[i] *= factor;
        estimated_standard_deviations_[i] *= factor;
    }
}

// ********************************************************************************

Angle PowderPattern::two_theta( const size_t i ) const
{
    if ( i < size() )
        return two_theta_values_[i];
    throw std::runtime_error( "PowderPattern::two_theta(): index out of bounds." );
}

// ********************************************************************************

double PowderPattern::intensity( const size_t i ) const
{
    if ( i < size() )
        return intensities_[i];
    throw std::runtime_error( "PowderPattern::intensity(): index out of bounds." );
}

// ********************************************************************************

double PowderPattern::estimated_standard_deviation( const size_t i ) const
{
    if ( i < size() )
        return estimated_standard_deviations_[i];
    throw std::runtime_error( "PowderPattern::estimated_standard_deviation(): index out of bounds." );
}

// ********************************************************************************

void PowderPattern::set_two_theta( const size_t i, const Angle value )
{
    if ( i < size() )
        two_theta_values_[i] = value;
    else
        throw std::runtime_error( "PowderPattern::set_two_theta(): index out of bounds." );
}

// ********************************************************************************

// ESD is NOT updated.
void PowderPattern::set_intensity( const size_t i, const double value )
{
    if ( i < size() )
        intensities_[i] = value;
    else
        throw std::runtime_error( "PowderPattern::set_intensity(): index out of bounds." );
}

// ********************************************************************************

void PowderPattern::set_estimated_standard_deviation( const size_t i, const double value )
{
    if ( i < size() )
        estimated_standard_deviations_[i] = value;
    else
        throw std::runtime_error( "PowderPattern::set_estimated_standard_deviation(): index out of bounds." );
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
void PowderPattern::set_two_theta_start( const Angle two_theta_start )
{
    if ( empty() )
        throw std::runtime_error( "PowderPattern::set_two_theta_start(): no data points." );
    if ( ( two_theta_values_[ 0 ] - 0.5*average_two_theta_step() ) > two_theta_start )
        throw std::runtime_error( "PowderPattern::set_two_theta_start(): adding data points not implemented yet." );
        
}

// ********************************************************************************

// Uses average_two_theta_step() to add new points. Intensities and ESDs are initialised to 0.0.
void PowderPattern::set_two_theta_end( const Angle two_theta_end )
{
    if ( empty() )
        throw std::runtime_error( "PowderPattern::set_two_theta_end(): no data points." );
    if ( ( two_theta_values_[ size()-1 ] + 0.5*average_two_theta_step() ) < two_theta_end )
        throw std::runtime_error( "PowderPattern::set_two_theta_end(): adding data points not implemented yet." );
        
}

// ********************************************************************************

void PowderPattern::reduce_range_to( const Angle two_theta_start, const Angle two_theta_end )
{
    size_t iStart = find_two_theta( two_theta_start );
    size_t iEnd = find_two_theta( two_theta_end );
    if ( iStart != 0 )
    {
        for ( size_t i( 0 ); i != iEnd + 1 - iStart; ++i )
        {
            two_theta_values_[ i ] = two_theta_values_[ i + iStart ];
            intensities_[ i ] = intensities_[ i + iStart ];
            estimated_standard_deviations_[ i ] = estimated_standard_deviations_[ i + iStart ];
        }
    }
    two_theta_values_.resize( iEnd - iStart );
    intensities_.resize( iEnd - iStart );
    estimated_standard_deviations_.resize( iEnd - iStart );
}

// ********************************************************************************

double PowderPattern::cumulative_intensity() const
{
    return add_doubles( intensities_ );
}

// ********************************************************************************

double PowderPattern::cumulative_intensity( const Angle two_theta_start, const Angle two_theta_end ) const
{
    size_t iStart = find_two_theta( two_theta_start );
    size_t iEnd = find_two_theta( two_theta_end );
    double result( 0.0 );
    for ( size_t i( iStart ); i != iEnd + 1 ; ++i )
        result += intensities_[i];
    return result;
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
            wavelength_ = Wavelength::determine_from_wavelength( string2double( words[0] ) );
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
            // There is one dummy data point with value -1 after the last valid data point.
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

// We have a couple of problems here:
// The cif may be huge and full of other stuff that we do not need (it may even have multiple structures
// with multiple powder diffraction patterns)
//_pd_meas_2theta_range_min  0.50000
//_pd_meas_2theta_range_max  49.99654
//_pd_meas_2theta_range_inc  0.00100
//_pd_meas_number_of_points  49575
//
//loop_
//   _pd_meas_intensity_total
//   _pd_calc_intensity_total
//   _pd_proc_intensity_bkg_calc
//   _pd_proc_ls_weight
//  878.115218   0            0           0         
//  845.112609   0            0           0         
//  830.491604   0            0           0         
//  848.287054   0            0           0         
//  874.178024   0            0           0         
//  834.93539    0            0           0         
//
//
//_pd_meas_2theta_range_min     5.0019
//_pd_meas_2theta_range_inc     0.0025
//_pd_meas_2theta_range_max     40.0019
//_pd_meas_number_of_points     14001
//
//loop_
//      _pd_meas_counts_total
//32393   32393   32357   32330   32321   32406   32472   32505   32492   32484   #   5.0519   
//32484   32519   32512   32454   32442   32433   32426   32386   32380   32410   #   5.0769   
//32487   32525   32526   32508   32489   32469   32518   32530   32509   32492   #   5.1019   
//32352   32361   32416   32393   32376   32386   32388   32390   32391   32392   #   5.0269   
//32309   32331   32335   32310   32293   32289   32307   32350   32377   32370   #   5.0019   
void PowderPattern::read_cif( const FileName & file_name )
{
    *this = PowderPattern();
    TextFileReader text_file_reader( file_name );
    text_file_reader.set_skip_empty_lines( true );
    std::vector< std::string > words;
    Angle two_theta_start;
    Angle two_theta_step;
    Angle two_theta_end;
    size_t number_of_points( 0 );
    bool found_pd_meas_2theta_range_min( false );
    bool found_pd_meas_2theta_range_max( false );
    bool found_pd_meas_2theta_range_inc( false );
    bool found_pd_meas_number_of_points( false );
    std::vector< double > counts;
    while ( text_file_reader.get_next_line( words ) )
    {
        if ( words[0] == "_pd_meas_2theta_range_min" )
        {
            if ( words.size() != 2 )
                throw std::runtime_error( "PowderPattern::read_cif(): keyword _pd_meas_2theta_range_min does not have value." );
            two_theta_start = Angle::from_degrees( string2double( words[1] ) );
            found_pd_meas_2theta_range_min = true;
            continue;
        }
        if ( words[0] == "_pd_meas_2theta_range_max" )
        {
            if ( words.size() != 2 )
                throw std::runtime_error( "PowderPattern::read_cif(): keyword _pd_meas_2theta_range_max does not have value." );
            two_theta_end = Angle::from_degrees( string2double( words[1] ) );
            found_pd_meas_2theta_range_max = true;
            continue;
        }
        if ( words[0] == "_pd_meas_2theta_range_inc" )
        {
            if ( words.size() != 2 )
                throw std::runtime_error( "PowderPattern::read_cif(): keyword _pd_meas_2theta_range_inc does not have value." );
            two_theta_step = Angle::from_degrees( string2double( words[1] ) );
            found_pd_meas_2theta_range_inc = true;
            continue;
        }
        if ( words[0] == "_pd_meas_number_of_points" )
        {
            if ( words.size() != 2 )
                throw std::runtime_error( "PowderPattern::read_cif(): keyword _pd_meas_number_of_points does not have value." );
            number_of_points = string2integer( words[1] );
            found_pd_meas_number_of_points = true;
            continue;
        }
        if ( words[0] == "loop_" )
        {
            if ( words.size() != 1 )
                throw std::runtime_error( "PowderPattern::read_cif(): loop_ should be the only keyword on a line." );
            if ( ! text_file_reader.get_next_line( words ) )
                throw std::runtime_error( "PowderPattern::read_cif(): loop_ not followed by values." );
            if ( words[0] != "_pd_meas_counts_total" )
                throw std::runtime_error( "PowderPattern::read_cif(): only _pd_meas_counts_total in loop_ has been implemented." );
            std::string line;
            while ( text_file_reader.get_next_line( line ) )
            {
                line = remove_from( line, '#' );
                words = split( line );
                try
                {
                    for ( size_t i( 0 ); i != words.size(); ++i )
                        counts.push_back( string2double( words[ i ] ) );
                }
                catch ( std::exception & e )
                {
                    break;
                }
            }
            continue;
        }
    }
    if ( ! found_pd_meas_2theta_range_min )
        throw std::runtime_error( "PowderPattern::read_cif(): keyword _pd_meas_2theta_range_min not found." );
    if ( ! found_pd_meas_2theta_range_max )
        throw std::runtime_error( "PowderPattern::read_cif(): keyword _pd_meas_2theta_range_max not found." );
    if ( ! found_pd_meas_2theta_range_inc )
        throw std::runtime_error( "PowderPattern::read_cif(): keyword _pd_meas_2theta_range_inc not found." );
    // Consistency check.
    if ( found_pd_meas_number_of_points )
    {
        if ( number_of_points != counts.size() )
            throw std::runtime_error( "PowderPattern::read_cif(): _pd_meas_number_of_points inconsistent with the actual number of points." );
        Angle two_theta_step_should_be = ( two_theta_end - two_theta_start ) / ( number_of_points - 1 );
    }
    if ( counts.size() < 2 )
        throw std::runtime_error( "PowderPattern::read_cif(): there are fewer than two points in the pattern." );
    Angle two_theta_step_should_be = ( two_theta_end - two_theta_start ) / ( counts.size() - 1 );
    std::cout << "two_theta_step_should_be = " << two_theta_step_should_be << std::endl;
    std::cout << "two_theta_step           = " << two_theta_step << std::endl;
    // @@ This is a poor algorithm: two_theta_step was probably rounded and we need the exact value.
    for ( size_t i( 0 ); i != counts.size(); ++i )
    {
        push_back( ( i * two_theta_step_should_be ) + two_theta_start, counts[ i ] );
    }
}

// ********************************************************************************

void PowderPattern::save_xye( const FileName & file_name, const bool include_wave_length ) const
{
    TextFileWriter text_file_writer( file_name );
    if ( include_wave_length )
        text_file_writer.write_line( double2string( wavelength_.wavelength_1() ) );
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

double PowderPattern::normalise_highest_peak( const double highest_peak )
{
    // Find the highest intensity.
    double max_intensity = calculate_maximum( intensities_ );
    if ( nearly_zero( max_intensity ) )
        throw std::runtime_error( "PowderPattern::normalise_highest_peak(): highest peak is 0.0." );
    double scale_factor = highest_peak / max_intensity;
    scale( scale_factor );
    return scale_factor;
}

// ********************************************************************************

// Normalises the total signal = area under the pattern = cumulative_intensity() .
double PowderPattern::normalise_total_signal( const double total_signal )
{
    double current_total_signal = cumulative_intensity();
    if ( nearly_zero( current_total_signal ) )
        throw std::runtime_error( "PowderPattern::normalise_total_signal(): total signal is 0.0." );
    // Scale to total_signal.
    double scale_factor = total_signal / current_total_signal;
    scale( scale_factor );
    return scale_factor;
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
            estimated_standard_deviations_[i] = 4.4;
        else if ( intensities_[i] > 10000.0 )
            estimated_standard_deviations_[i] = intensities_[i] / 100.0;
        else
            estimated_standard_deviations_[i] = sqrt( intensities_[i] );
    }
}

// ********************************************************************************

// I(fixed slit) = I(variable slit) / sin(theta).
void PowderPattern::convert_to_fixed_slit()
{
    for ( size_t i( 0 ); i != size(); ++i )
    {
        intensities_[i] /= ( two_theta_values_[i] / 2.0 ).sine(); // @@ We should check for divide by zero
        estimated_standard_deviations_[i] /= ( two_theta_values_[i] / 2.0 ).sine();
    }
}

// ********************************************************************************

// I(fixed slit) = I(variable slit) / sin(theta).
void PowderPattern::convert_to_variable_slit()
{
    for ( size_t i( 0 ); i != size(); ++i )
    {
        intensities_[i] *= ( two_theta_values_[i] / 2.0 ).sine();
        estimated_standard_deviations_[i] *= ( two_theta_values_[i] / 2.0 ).sine();
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
        intensities_[i] = Poisson_distribution( round_to_int( intensities_[i] ) );
}

// ********************************************************************************

void PowderPattern::add_Poisson_noise_including_zero( const size_t threshold )
{
    for ( size_t i( 0 ); i != size(); ++i )
    {
        if ( intensities_[i] < threshold )
        {
            int intensity = round_to_int( intensities_[i] ) + threshold;
            intensity = Poisson_distribution( intensity ) - static_cast<int>(threshold);
            intensities_[i] = std::abs( intensity );
        }
        else
            intensities_[i] = Poisson_distribution( round_to_int( intensities_[i] ) );
    }
}

// ********************************************************************************

void PowderPattern::make_counts_integer()
{
    for ( size_t i( 0 ); i != size(); ++i )
        intensities_[i] = round_to_int( intensities_[i] );
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

// Check that the two patterns have the same range and 2theta step.
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
    if ( l < Angle() )
        throw std::runtime_error( "weighted_cross_correlation( PowderPattern, PowderPattern, Angle ): l must be non-negative." );
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
                double w = 1.0 - absolute( j ) / static_cast<double>( m );
                if ( (true) )
                    result += w * lhs.intensity( i ) * rhs.intensity( i + j );
                else
                    result += w * ( lhs.intensity( i ) / lhs.estimated_standard_deviation( i ) ) * ( rhs.intensity( i + j ) / rhs.estimated_standard_deviation( i + j ) );
            }
        }
    }
    return result;
}

// ********************************************************************************

double normalised_weighted_cross_correlation( const PowderPattern & lhs, const PowderPattern & rhs, Angle l )
{
    if ( ! same_range( lhs, rhs ) )
        throw std::runtime_error( "normalised_weighted_cross_correlation( const PowderPattern &, const PowderPattern & ): ranges not same." );
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
    if ( powder_pattern.empty() )
        return powder_pattern;
    PowderPattern result( powder_pattern );
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
            result.set_intensity( i, new_value / ( 2.0 * smoothing_window + 1.0 ) );
        }
    }
    if ( true )
    {
        RunningAverageAndESD< double > I_average;
        double I_minimum = result.intensity( 0 );
        for ( size_t i( 0 ); i < size; ++i )
        {
            if ( result.intensity( i ) < I_minimum )
                I_minimum = result.intensity( i );
            I_average.add_value( result.intensity( i ) );
        }
        for ( size_t i( 0 ); i < size; ++i )
        {
            if ( result.intensity( i ) > ( I_average.average() + 2.0 * ( I_average.average() - I_minimum ) ) )
                result.set_intensity( i, ( I_average.average() + 2.0 * ( I_average.average() - I_minimum ) ) );
        }
    }
    for ( size_t iter( 0 ); iter < niterations; ++iter )
    {
        PowderPattern pp_old = result;
        for ( size_t i( 0 ); i < size; ++i )
        {
            double average_value( 0.0 );
            for ( size_t j( 1 ); j <= window; ++j )
            {
                average_value += pp_old.intensity( std::max( int(i)-int(j), int(0) ) );
                average_value += pp_old.intensity( std::min( i+j, size-1 ) );
            }
            average_value /= 2.0 * window;
            result.set_intensity( i, average_value );
        }
        for ( size_t i( 0 ); i < size; ++i )
            result.set_intensity( i, std::min( pp_old.intensity( i ), result.intensity( i ) ) );
    }
    return result;
}

// ********************************************************************************

PowderPattern calculate_Poisson_noise( const PowderPattern & powder_pattern )
{
    PowderPattern result;
    result.set_wavelength( powder_pattern.wavelength() );
    result.reserve( powder_pattern.size() );
    for ( size_t i( 0 ); i != powder_pattern.size(); ++i )
        result.push_back( powder_pattern.two_theta( i ), Poisson_distribution( round_to_int( powder_pattern.intensity( i ) ) ) - powder_pattern.intensity( i ), 0.0 );
    return result;
}

// ********************************************************************************

// If the number of counts is less than threshold, adds threshold, then calculates the Poisson noise,
// then subtracts the threshold, then makes the remainder positive.
// If the maximum is scaled to be 10,000 counts, a good threshold value is 20.
PowderPattern calculate_Poisson_noise_including_zero( const PowderPattern & powder_pattern, const size_t threshold )
{
    PowderPattern result;
    result.set_wavelength( powder_pattern.wavelength() );
    result.reserve( powder_pattern.size() );
    for ( size_t i( 0 ); i != powder_pattern.size(); ++i )
    {
        int old_intensity = powder_pattern.intensity( i );
        int new_intensity;
        if ( old_intensity < threshold )
            new_intensity = std::abs( Poisson_distribution( old_intensity + threshold ) - static_cast<double>(threshold) );
        else
            new_intensity = Poisson_distribution( old_intensity );
        result.push_back( powder_pattern.two_theta( i ), new_intensity - old_intensity, 0.0 );
    }
    return result;
}

// ********************************************************************************

PowderPattern add_powder_patterns( const std::vector< PowderPattern > & powder_patterns, const std::vector< double > & noscp2ts )
{
    if ( powder_patterns.empty() )
        throw std::runtime_error( "add_powder_patterns(): Error: no powder patterns provided." );
    if ( powder_patterns.size() != noscp2ts.size() )
        throw std::runtime_error( "add_powder_patterns(): Error: powder_patterns and noscp2ts not the same size." );
    // Check that they have the same wavelength and average_two_theta_step.
    for ( size_t i( 1 ); i != powder_patterns.size(); ++i )
    {
        if ( ! nearly_equal( powder_patterns[0].wavelength(), powder_patterns[i].wavelength() ) )
            throw std::runtime_error( "add_powder_patterns(): Error: wavelengths not the same." );
        if ( ! nearly_equal( powder_patterns[0].average_two_theta_step(), powder_patterns[i].average_two_theta_step() ) )
            throw std::runtime_error( "add_powder_patterns(): Error: average_two_theta_step not the same." );
    }
    // Find the smallest and largest 2theta values.
    Angle two_theta_min = powder_patterns[0].two_theta_start();
    Angle two_theta_max = powder_patterns[0].two_theta_end();
    for ( size_t i( 1 ); i != powder_patterns.size(); ++i )
    {
        if ( powder_patterns[i].two_theta_start() < two_theta_min )
            two_theta_min = powder_patterns[i].two_theta_start();
        if ( powder_patterns[i].two_theta_end() > two_theta_max )
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
            // @@ The following is wrong because "result" includes the min and max 2theta, which may not exist in the current pattern.
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

std::vector< PowderPattern > split( const PowderPattern & powder_pattern, const size_t n, const bool recalculate_ESDs )
{
    std::vector< PowderPattern > result;
    for ( size_t j( 0 ); j != n; ++j )
    {
        result.push_back( PowderPattern() );
        result[j].set_wavelength( powder_pattern.wavelength() );
    }
    RandomNumberGenerator_integer rng;
    for ( size_t i( 0 ); i != powder_pattern.size(); ++i )
    {
        int old_intensity_int = round_to_int( powder_pattern.intensity( i ) );
        if ( old_intensity_int < 0 )
            throw std::runtime_error( "split(PowderPattern): Error: intensity is negative." );
        size_t old_intensity = old_intensity_int;
        double old_ESD = powder_pattern.estimated_standard_deviation( i ) / n;
        size_t target_average = round_to_int( powder_pattern.intensity( i ) / n );
        // If round_to_int( powder_pattern.intensity( i ) / n ) == 0 then Poisson_distribution( target_average ); returns 0.
        std::vector< size_t > intensities;
        size_t sum( 0 );
        if ( target_average == 0 )
        {
            if ( old_intensity > ( n / 2 ) )
            {
                sum = n;
                for ( size_t j( 0 ); j != n; ++j )
                    intensities.push_back( 1 );
            }
            else
            {
                for ( size_t j( 0 ); j != n; ++j )
                    intensities.push_back( 0 );
            }
        }
        else
        {
            for ( size_t j( 0 ); j != n; ++j )
            {
                size_t new_intensity = Poisson_distribution( target_average );
                sum += new_intensity;
                intensities.push_back( new_intensity );
            }
            if ( sum != 0 )
            {
                // We scale all intensities.
                double scale_factor = powder_pattern.intensity( i ) / sum;
                sum = 0;
                for ( size_t j( 0 ); j != n; ++j )
                {
                    intensities[j] = round_to_int( intensities[j] * scale_factor );
                    sum += intensities[j];
                }
            }
        }
        if ( sum != old_intensity )
        {
            // We randomly add / remove counts until we have reached the target.
            int step = ( sum > old_intensity ) ? -1 : 1;
            while ( sum != old_intensity )
            {
                // Randomly pick a pattern.
                size_t j = rng.next_number( 0, n-1 );
                if ( ! ( ( step == -1 ) && ( intensities[j] == 0 ) ) )
                {
                    intensities[j] += step;
                    sum += step;
                }
            }
        }
        for ( size_t j( 0 ); j != n; ++j )
        {
            if ( recalculate_ESDs )
                result[j].push_back( powder_pattern.two_theta( i ), intensities[j] );
            else
                result[j].push_back( powder_pattern.two_theta( i ), intensities[j], old_ESD );
        }
    }
    return result;
}

// ********************************************************************************

