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

#include "Histogram.h"
#include "MathsFunctions.h"
#include "StringFunctions.h"

#include <stdexcept>

#include <iostream> // For debugging

// ********************************************************************************

Histogram::Histogram( const double start, const double finish, const size_t number_of_bins ) :
start_(start),
finish_(finish),
number_of_bins_(number_of_bins),
lower_than_start_(0),
greater_than_finish_(0)
{
    if ( number_of_bins == 0 )
        throw std::runtime_error( "Histogram::Histogram(): number_of_bins cannot be 0." );
    if ( start >= finish )
        throw std::runtime_error( "Histogram::Histogram(): start must be lower than finish." );
    data_ = std::vector< size_t >( number_of_bins_, 0.0 );
}

// ********************************************************************************

void Histogram::add_data( const std::vector< double > & data )
{
    for ( size_t i( 0 ); i != data.size(); ++i )
        add_data( data[i] );
}

// ********************************************************************************

void Histogram::add_data( const double data )
{
    if ( data == start_ )
        ++data_[0];
    else if ( data < start_ )
        ++lower_than_start_;
    else if ( data > finish_ )
        ++greater_than_finish_;
    else if ( data == finish_ )
        ++data_[number_of_bins_-1];
    else
    {
        size_t index = round_to_size_t( ( number_of_bins_ * ( ( data - start_ ) / ( finish_ - start_ ) ) ) - 0.5 );
        ++data_[index];
    }
}

// ********************************************************************************

size_t Histogram::bin( const size_t i ) const
{
    if ( i < number_of_bins_ )
        return data_[i];
    throw std::runtime_error( "Histogram::bin(): index out of range." );
}

// ********************************************************************************

double Histogram::middle_of_bin( const size_t i ) const
{
    if ( i < number_of_bins_ )
    {
        double bin_size = ( finish_ - start_ ) / static_cast<double>(number_of_bins_);
        return start_ + ( bin_size * ( static_cast<double>(i) + 0.5 ) );
    }
    throw std::runtime_error( "Histogram::middle_of_bin(): index out of range." );
}

// ********************************************************************************

size_t Histogram::maximum() const
{
    return calculate_maximum( data_ );
}

// ********************************************************************************

void Histogram::show() const
{
    for ( size_t i( 0 ); i != data_.size(); ++i )
        std::cout << middle_of_bin( i ) << " " << data_[i] << std::endl;
}

// ********************************************************************************

void Histogram::plot() const
{
    double maximum = this->maximum();
    for ( size_t i( 0 ); i != data_.size(); ++i )
        std::cout << "|" << make_multiple( "#", round_to_int( ( static_cast<double>(data_[i]) / maximum ) * 80.0 ) ) << std::endl;
}

// ********************************************************************************

std::vector< double > Histogram::x_values() const
{
    std::vector< double > result;
    for ( size_t i( 0 ); i != data_.size(); ++i )
        result.push_back( middle_of_bin( i ) );
    return result;
}

// ********************************************************************************

std::vector< double > Histogram::y_values() const
{
    std::vector< double > result;
    for ( size_t i( 0 ); i != data_.size(); ++i )
        result.push_back( data_[i] );
    return result;
}

// ********************************************************************************

void Histogram::values( const size_t iStart, const size_t iEnd, std::vector< double > & x_values, std::vector< double > & y_values ) const
{
    for ( size_t i( iStart ); i != iEnd+1; ++i )
    {
        x_values.push_back( middle_of_bin( i ) );
        y_values.push_back( data_[i] );
    }
}

// ********************************************************************************

