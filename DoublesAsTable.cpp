/* *********************************************
Copyright (c) 2013-2024, Cornelis Jan (Jacco) van de Streek
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

#include "DoublesAsTable.h"
#include "BasicMathsFunctions.h"
#include "MathsFunctions.h"
#include "StringFunctions.h"
#include "Utilities.h"

//#include <iostream>
//#include <stdexcept>

// ********************************************************************************

std::vector< std::string > doubles_as_table( const std::vector< double > & doubles,
                                             const size_t ndecimals,
                                             const std::string & separator_begin,
                                             const std::string & separator_middle,
                                             const std::string & separator_end,
                                             const size_t target_length,
                                             const bool force_minus_sign_padding,
                                             size_t ncolumns )
{
    std::vector< std::string > result;
    if ( doubles.empty() )
        return result;
    // Establish the maximum length.
    bool need_minus_sign( force_minus_sign_padding );
    size_t maximum_length( 0 );
    for ( size_t i( 0 ); i != doubles.size(); ++i )
    {
        size_t length;
        std::string word = double2string_2( doubles[i], ndecimals );
        length = word.length();
        if ( doubles[i] < 0.0 )
        {
            need_minus_sign = true;
            --length;
        }
        if ( length > maximum_length )
            maximum_length = length;
    }
    if ( need_minus_sign )
        ++maximum_length;
    // We cannot have less than one column.
    // Whatever is left we can fill with values. Each value has a length maximum_length + separator_middle.length().
    // If ncolumns has been supplied, use it, otherwise calculate it based on target_length.
    if ( ncolumns == 0 )
    {
        if ( target_length < ( separator_begin.length() + maximum_length + separator_end.length() ) )
            ncolumns = 1;
        else
            ncolumns = 1 + round_to_size_t( ( target_length - separator_begin.length() - maximum_length - separator_end.length() ) / ( maximum_length + separator_middle.length() ) );
    }
    size_t nrows = ( doubles.size() + ncolumns - 1 ) / ncolumns;
    size_t k( 0 );
    for ( size_t i( 0 ); i != nrows; ++i )
    {
        std::string line( separator_begin );
        for ( size_t j( 0 ); j != ncolumns-1; ++j )
        {
            if ( k < doubles.size() )
                line += double2string( doubles[k], ndecimals, maximum_length );
            else
                line += make_multiple( " ", maximum_length );
            line += separator_middle;
            ++k;
        }
        if ( k < doubles.size() )
            line += double2string( doubles[k], ndecimals, maximum_length );
        else
            line += make_multiple( " ", maximum_length );
        line += separator_end;
        ++k;
        result.push_back( line );
    }
    return result;
}

// ********************************************************************************

std::vector< std::string > doubles_as_table( const std::vector< double > & doubles,
                                             const size_t ndecimals,
                                             const std::string & separator,
                                             const size_t target_length,
                                             const bool force_minus_sign_padding,
                                             size_t ncolumns )
{
    std::vector< std::string > result;
    if ( doubles.empty() )
        return result;
    // Establish the maximum length.
    bool need_minus_sign( force_minus_sign_padding );
    size_t maximum_length( 0 );
    for ( size_t i( 0 ); i != doubles.size(); ++i )
    {
        size_t length;
        std::string word = double2string_2( doubles[i], ndecimals );
        length = word.length();
        if ( doubles[i] < 0.0 )
        {
            need_minus_sign = true;
            --length;
        }
        if ( length > maximum_length )
            maximum_length = length;
    }
    if ( need_minus_sign )
        ++maximum_length;
    // The terminating separator on each line must have its trailing whitespace removed.
    std::string separator_2 = remove_trailing_whitespace( separator );
    // We cannot have less than one column.
    // Whatever is left we can fill with values. Each value has a length maximum_length + separator_middle.length().
    // If ncolumns has been supplied, use it, otherwise calculate it based on target_length.
    if ( ncolumns == 0 )
    {
        if ( target_length < ( maximum_length + separator_2.length() ) )
            ncolumns = 1;
        else
            ncolumns = 1 + round_to_size_t( ( target_length - maximum_length - separator_2.length() ) / ( maximum_length + separator.length() ) );
    }
    size_t nrows = ( doubles.size() + ncolumns - 1 ) / ncolumns;
    size_t k( 0 );
    for ( size_t i( 0 ); i != nrows; ++i )
    {
        std::string line;
        for ( size_t j( 0 ); j != ncolumns; ++j )
        {
            if ( k < doubles.size() )
            {
                line += double2string( doubles[k], ndecimals, maximum_length );
                if ( ( k + 1 ) < doubles.size() )
                {
                    if ( ( k || ncolumns ) == 0 )
                        line += separator_2;
                    else
                        line += separator;
                }
            }
            ++k;
        }
        result.push_back( line );
    }
    return result;    
}

// ********************************************************************************

std::vector< std::string > generate_Gauss_Legendre_quadrature_code( const double x1,
                                                                    const double x2,
                                                                    const size_t npoints,
                                                                    const size_t ndecimals,
                                                                    const size_t target_length,
                                                                    const bool force_minus_sign_padding,
                                                                    const size_t ncolumns )
{
    std::vector< std::string > result;
    std::vector< double > x;
    std::vector< double > w;
    std::string padding( "                    " );
    Gauss_Legendre_quadrature( x1, x2, npoints, x, w );
    result.push_back( std::string( "static const double x[" + size_t2string( npoints ) + "] =" ) );
    result.push_back( std::string( padding + "{" ) );
    std::vector< std::string > values = doubles_as_table( x, ndecimals, ", ", target_length - padding.length() - 4, force_minus_sign_padding, ncolumns );
    for ( size_t i( 0 ); i != values.size(); ++i )
        result.push_back( padding + "    " + values[i] );
    result.push_back( std::string( padding + "};" ) );
    result.push_back( std::string( "" ) );
    result.push_back( std::string( "static const double w[" + size_t2string( npoints ) + "] =" ) );
    result.push_back( std::string( padding + "{" ) );
    values = doubles_as_table( w, ndecimals, ", ", target_length - padding.length() - 4, force_minus_sign_padding, ncolumns );
    for ( size_t i( 0 ); i != values.size(); ++i )
        result.push_back( padding + "    " + values[i] );
    result.push_back( std::string( padding + "};" ) );
    return result;
}

// ********************************************************************************

