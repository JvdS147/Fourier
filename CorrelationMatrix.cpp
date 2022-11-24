/* *********************************************
Copyright (c) 2013-2022, Cornelis Jan (Jacco) van de Streek
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

#include "CorrelationMatrix.h"
#include "TextFileWriter.h"
#include "Utilities.h"

#include <fstream>
#include <iostream>
#include <stdexcept>

// ********************************************************************************

CorrelationMatrix::CorrelationMatrix( const size_t dimension ) : dimension_(dimension)
{
    data_ptr_ = new double[ ( dimension_ * ( dimension_ - 1 ) ) / 2 ];
    value_on_diagonal_ = 1.0;
}

// ********************************************************************************

CorrelationMatrix::~CorrelationMatrix()
{
    delete[] data_ptr_;
}

// ********************************************************************************

double CorrelationMatrix::value( size_t i, size_t j ) const
{
    if ( i < j )
        std::swap( i, j );
    if ( dimension_ < (i+1) )
        throw std::runtime_error( "CorrelationMatrix::element(): out of bounds ( " + size_t2string(i) + " > " + size_t2string(dimension_) + " )" );
    if ( i == j )
        return value_on_diagonal_;
    return data_ptr_[ ((i*(i-1))/2) + j ];
}

// ********************************************************************************

void CorrelationMatrix::set_value( size_t i, size_t j, const double value )
{
    if ( i < j )
        std::swap( i, j );
    if ( i > (dimension_-1) )
        throw std::runtime_error( "CorrelationMatrix::set_element(): out of bounds ( " + size_t2string(i) + " > " + size_t2string(dimension_) + " )" );
    if ( i == j )
        return;
    data_ptr_[ ((i*(i-1))/2) + j ] = value;
}

// ********************************************************************************

// Diagonal is not included
double CorrelationMatrix::largest_value() const
{
    double result( 0.0 );
    for ( size_t i( 0 ); i != (( dimension_ * ( dimension_ - 1 ) ) / 2); ++i )
    {
        if ( data_ptr_[i] > result )
            result = data_ptr_[i];
    }
    return result;
}

// ********************************************************************************

// Diagonal is not included
double CorrelationMatrix::smallest_value() const
{
    double result( 1.0 );
    for ( size_t i( 0 ); i != (( dimension_ * ( dimension_ - 1 ) ) / 2); ++i )
    {
        if ( data_ptr_[i] < result )
            result = data_ptr_[i];
    }
    return result;
}

// ********************************************************************************

void CorrelationMatrix::save( const FileName & file_name ) const
{
    TextFileWriter text_file_writer( file_name );
    for ( size_t i( 0 ); i != size(); ++i )
    {
        for ( size_t j( 0 ); j != size(); ++j )
        {
            text_file_writer.write( double2string( value( i, j ) ) + "  " );
        }
        text_file_writer.write_line();
    }
}

// ********************************************************************************

// Divides the entries into clusters, all entries more similar than threshold are put into a cluster.
// If this leads to inconsistencies, e.g. because A = B, B = C, but A != C, then A = C.
std::vector< std::vector< size_t > > CorrelationMatrix::clusters( const double threshold ) const
{
    std::vector< std::vector< size_t > > result;
    std::vector< bool > done( size(), false );
    {
        for ( size_t i( 0 ); i != size(); ++i )
        {
            if ( done[i] )
                continue;
            std::vector< size_t > one_cluster;
            one_cluster.push_back( i );
            done[i] = true;
            for ( size_t j( i+1 ); j < size(); ++j )
            {
                if ( done[j] )
                    continue;
                if ( value( i, j ) > threshold  )
                {
                    one_cluster.push_back( j );
                    done[j] = true;
//                    // The following is to ensure that if A = B and B = C, then A = C.
//                    // This may not necessarily follow from the correlation value,
//                    // because it is possible that A and C are the beginning and the end
//                    // of a temperature series.
//                    for ( size_t k( j+1 ); k < size(); ++k )
//                    {
//                        if ( done[k] )
//                            continue;
//                        if ( ( value( j, k ) > threshold ) && ( ! ( value( i, k ) > threshold) )  )
//                        {
//                            one_cluster.push_back( k );
//                            done[k] = true;
//                        }
//                    }
                }
            }
            result.push_back( one_cluster );
        }
    }
    // @@ We may still have inconsistencies here.
    // @@ Depending on how we handle the inconsistencies, the lists of indices may not be ordered when we get here.
    return result;
}

// ********************************************************************************

std::vector< std::vector< size_t > > CorrelationMatrix::clusters( const double grey_area_threshold, const double threshold ) const
{
    std::vector< std::vector< size_t > > result;
    std::vector< bool > done( size(), false );
    {
        for ( size_t i( 0 ); i != size(); ++i )
        {
            if ( done[i] )
                continue;
            std::vector< size_t > one_cluster;
            one_cluster.push_back( i );
            done[i] = true;
            for ( size_t j( i+1 ); j < size(); ++j )
            {
                if ( done[j] )
                    continue;
                if ( value( i, j ) > threshold  )
                {
                    one_cluster.push_back( j );
                    done[j] = true;
                }
                else if ( value( i, j ) > grey_area_threshold  )
                {
                    std::cout << "Warning: similarity " << i << ", " << j << " is " << value( i, j ) << std::endl;
                }
            }
            result.push_back( one_cluster );
        }
    }
    // @@ We may still have inconsistencies here.
    // @@ Depending on how we handle the inconsistencies, the lists of indices may not be ordered when we get here.
    return result;
}

// ********************************************************************************

