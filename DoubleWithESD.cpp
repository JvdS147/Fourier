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

#include "DoubleWithESD.h"
#include "StringFunctions.h"
#include "Utilities.h"

#include <cmath>
#include <stdexcept>
#include <sstream>

// ********************************************************************************

// At the moment, can only cope with "0.12345(6)"
DoubleWithESD::DoubleWithESD( std::string input )
{
    if ( input.empty() )
        throw std::runtime_error( "DoubleWithESD::DoubleWithESD( const std::string & ): input is empty." );
    // Check if there is a "("
    size_t iPos = input.find_first_of( '(' );
    // Check if last character is a ")"
    bool ends_in_parenthesis = ( input[ input.length()-1 ] == ')' );
    if ( iPos == std::string::npos )
    {
        if ( ends_in_parenthesis )
            throw std::runtime_error( "DoubleWithESD::DoubleWithESD( const std::string & ): no opening parenthesis." );
        no_esd_ = true;
        estimated_standard_deviation_ = 0.0;
    }
    else
    {
        if ( ! ends_in_parenthesis )
            throw std::runtime_error( "DoubleWithESD::DoubleWithESD( const std::string & ): no closing parenthesis." );
        no_esd_ = false;
        std::string esd_str = input.substr( iPos+1, (input.length()-1)-(iPos+1) );
        input.erase( iPos );
        iPos = input.find_first_of( '.' );
        if ( iPos == std::string::npos )
            throw std::runtime_error( "DoubleWithESD::DoubleWithESD( const std::string & ): value must be floating point." );
        estimated_standard_deviation_ = string2integer( esd_str ) * std::pow( 10.0, -static_cast<double>(input.length() - 1 - iPos) );
    }
    value_ = string2double( input );
}

// ********************************************************************************

std::string DoubleWithESD::crystallographic_style() const
{
    if ( no_esd_ )
        return double2string( value_, 6 );
    return ::crystallographic_style( value_, estimated_standard_deviation_ );
}

// ********************************************************************************

std::string crystallographic_style( const double value, const double estimated_standard_deviation )
{
    if ( estimated_standard_deviation <= 0.0 )
        throw std::runtime_error( "crystallographic_style(): estimated_standard_deviation must be > 0.0." );
//    if ( std::abs(value) <= estimated_standard_deviation )
//        return double2string( value, 6 ) + "(99)";
    std::string result;
    std::string esd_str;
    {
    std::stringstream string_stream;
    string_stream << std::scientific;
    string_stream.precision( 1 );
    string_stream << estimated_standard_deviation;
    esd_str = string_stream.str();
    }
    // esd_str is now for example "1.9e-010" / "4.0e-010"
    // esd_str is now for example "1.9e+000" / "4.0e+000"
    // esd_str is now for example "1.9e+010" / "4.0e+010"
    if ( esd_str[0] != '1' )
    {
        std::stringstream string_stream;
        string_stream << std::scientific;
        string_stream.precision( 0 );
        string_stream << estimated_standard_deviation;
        esd_str = string_stream.str();
    }
    // esd_str is now for example "1.9e-010" / "4e-010"
    // esd_str is now for example "1.9e+000" / "4e+000"
    // esd_str is now for example "1.9e+010" / "4e+010"
    // Or "1e-010"!!!
    if ( ( esd_str[0] == '1' ) && ( esd_str[1] != '.' ) )
        esd_str = "1.0" + esd_str.substr(1);
    std::string precision_str;
    if ( esd_str[0] == '1' )
        precision_str = esd_str.substr(5,3);
    else
        precision_str = esd_str.substr(3,3);
    if ( estimated_standard_deviation < 2.0 )
    {
        // @@ does estimated_standard_deviation = "1.9e+000" work?
        size_t precision = string2integer( precision_str );
        std::stringstream string_stream;
        string_stream << std::fixed;
        if ( esd_str[0] == '1' )
            ++precision;
        string_stream.precision( precision );
        string_stream << value;
        result = string_stream.str() + "(" + esd_str.substr(0,1);
        if ( esd_str[0] == '1' )
            result += esd_str.substr(2,1);
        result += ")";
    }
    else // estimated_standard_deviation >= 2.0
    {
        // esd_str is now "4e+000", "1.9e+010" or "4e+010" or 1e+010
        std::stringstream string_stream;
        string_stream << std::fixed;
        string_stream.precision( 0 );
        string_stream << value;
        result = string_stream.str();
        size_t esd_pow = string2integer( esd_str.substr( esd_str.length()-3, std::string::npos ) );
        // @@ The following is wrong, it should round, not cut
        result = result.substr( 0, result.length() - esd_pow );
        result += make_multiple( "0", esd_pow );
        if ( esd_str[0] == '1' )
            esd_str = esd_str.substr( 0, 1 ) + esd_str.substr( 2, 1 ) + make_multiple( "0", esd_pow-1 );
        else
            esd_str = esd_str[0] + make_multiple( "0", esd_pow );
        result += "(" + esd_str + ")";
    }
    return result;
}

// ********************************************************************************

std::string crystallographic_style_pad_plus( const double value, const double estimated_standard_deviation )
{
    throw std::runtime_error( "crystallographic_style_pad_plus(): not yet implemented." );
}

// ********************************************************************************

