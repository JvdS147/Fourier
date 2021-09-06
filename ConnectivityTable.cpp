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
    * Neither the name of the University of Copenhagen nor the
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

#include "ConnectivityTable.h"
#include "Utilities.h"

#include <iostream>
#include <stdexcept>

// ********************************************************************************

ConnectivityTable::ConnectivityTable( const size_t natoms ) : dimension_(natoms)
{
    data_ = std::vector< size_t >( ( dimension_ * ( dimension_ - 1 ) ) / 2, 0 );
}

// ********************************************************************************

size_t ConnectivityTable::value( size_t i, size_t j ) const
{
    if ( i < j )
        std::swap( i, j );
    if ( i < dimension_ )
    {
        if ( i == j )
            return 1;
        return data_[ ((i*(i-1))/2) + j ];
    }
    else
        throw std::runtime_error( "ConnectivityTable::value(): out of bounds ( " + size_t2string(i) + " > " + size_t2string(dimension_) + " )" );
}

// ********************************************************************************

void ConnectivityTable::set_value( size_t i, size_t j, const size_t value )
{
    if ( i < j )
        std::swap( i, j );
    if ( i < dimension_ )
    {
        if ( i == j )
            return;
        data_[ ((i*(i-1))/2) + j ] = value;
    }
    else
        throw std::runtime_error( "ConnectivityTable::set_value(): out of bounds ( " + size_t2string(i) + " > " + size_t2string(dimension_) + " )" );
}

// ********************************************************************************

void ConnectivityTable::show() const
{
    for ( size_t i( 0 ); i != size(); ++i )
    {
        for ( size_t j( 0 ); j != size(); ++j )
        {
            std::cout << value( i, j ) << " ";
        }
        std::cout << std::endl;
    }
}

// ********************************************************************************

std::vector< bool > initialise_row( const ConnectivityTable & connectivity_table, const size_t iRow )
{
    std::vector< bool > result( connectivity_table.size(), false ); 
    for ( size_t j( 0 ); j != connectivity_table.size(); ++j )
        result[j] = ( connectivity_table.value( iRow, j ) > 0 );
    return result;
}

// ********************************************************************************

void show( const std::vector< bool > & input )
{
    for ( size_t i( 0 ); i != input.size(); ++i )
    {
        if ( input[i] )
            std::cout << "1";
        else
            std::cout << "0";        
    }
}

// ********************************************************************************

std::vector< std::vector< size_t > > split( const ConnectivityTable & connectivity_table )
{
    std::vector< std::vector< size_t > > result;
    std::vector< bool > done( connectivity_table.size(), false );
    for ( size_t i( 0 ); i != connectivity_table.size(); ++i )
    {
        if ( done[ i ] )
            continue;
        std::vector< bool > this_molecule = initialise_row( connectivity_table, i );
        done[ i ] = true;
        bool change( false );
        do
        {
            change = false;
            for ( size_t j( 0 ); j != connectivity_table.size(); ++j )
            {
                if ( done[ j ] )
                    continue;
                if ( this_molecule[j] )
                {
                    // Add row j to the current molecule
                    for ( size_t k( 0 ); k != connectivity_table.size(); ++k )
                    {
                        if ( done[ k ] )
                            continue;
                        if ( ( connectivity_table.value( j, k ) > 0 ) && ( ! this_molecule[k] ) )
                        {
                            this_molecule[k] = true;
                            change = true;
                        }
                    }                    
                    done[ j ] = true;
                }
            }
        }
        while ( change );
        std::vector< size_t > this_molecule_2;
        for ( size_t j( 0 ); j != connectivity_table.size(); ++j )
        {
            if ( this_molecule[j] )
            {
                this_molecule_2.push_back( j );             
            }
        }
        result.push_back( this_molecule_2 );
    }
    return result;
}

// ********************************************************************************

