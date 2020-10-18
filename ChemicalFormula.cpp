/* *********************************************
Copyright (c) 2013-2020, Cornelis Jan (Jacco) van de Streek
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

#include "ChemicalFormula.h"
#include "Utilities.h"

#include <stdexcept>
#include <algorithm>

// ********************************************************************************

// @@ There are all sorts of special cases, e.g. C10C10H10. We do not deal with any of them.
ChemicalFormula::ChemicalFormula( const std::string & input ) :
sort_by_atomic_number_(false)
{
    if ( input.empty() )
        return;
    if ( ! is_upper_case_letter( input[0] ) )
        throw std::runtime_error( "ChemicalFormula::ChemicalFormula(): start of new element is not an upper case letter." );
    size_t iPos( 0 );
    while ( iPos < input.length() )
    {
        // Current character is an upper case letter
        // Find the next upper case letter
        size_t iPos2( iPos + 1 );
        while ( ( iPos2 < input.length() ) && ( ! is_upper_case_letter( input[iPos2] ) ) )
            ++iPos2;
        // We must now have C, Cc, Cdd or Ccdd, where C is upper case letter, c is lower case letter, d is digit.
        Element element;
        size_t count( 1 );
        if ( ( iPos2 > ( iPos + 1 ) ) && is_lower_case_letter( input[iPos+1] ) )
        {
            element = Element( input.substr( iPos, 2 ) );
            iPos += 2;
        }
        else
        {
            element = Element( input.substr( iPos, 1 ) );
            ++iPos;
        }
        if ( ( iPos2 - iPos ) > 0 )
            count = string2integer( input.substr( iPos, iPos2-iPos ) );
        if ( contains( element ) )
            throw std::runtime_error( "ChemicalFormula::ChemicalFormula(): element is present more than once." );
        elements_[ element ] = count;
        iPos = iPos2;
    }
}

// ********************************************************************************

void ChemicalFormula::add_element( const Element element )
{
    if ( contains( element ) )
        ++elements_[ element ];
    else
        elements_[ element ] = 1;
}

// ********************************************************************************

void ChemicalFormula::multiply( const Fraction fraction )
{
    for ( std::map< Element, size_t >::iterator pos( elements_.begin() ); pos != elements_.end(); ++pos )
    {
        Fraction result = pos->second * fraction;
        if ( ! result.is_integer() )
            throw std::runtime_error( "ChemicalFormula::multiply(): result is not an integer." );
        pos->second = result.integer_part();
    }
}

// ********************************************************************************

// See sort_by_atomic_number()
Element ChemicalFormula::element( const size_t i ) const
{
    if ( sort_by_atomic_number_ )
    {
        std::map< Element, size_t >::const_iterator pos = elements_.begin();
        size_t j( 0 );
        while ( j != i )
        {
            ++j;
            ++pos;
        }
        return pos->first;
    }
    else
    {
        std::vector< Element > elements;
        for ( std::map< Element, size_t >::const_iterator pos( elements_.begin() ); pos != elements_.end(); ++pos )
        {
            elements.push_back( pos->first );
        }
        std::sort( elements.begin(), elements.end(), elements_less );
        return elements[i];
    }
}

// ********************************************************************************

// See sort_by_atomic_number()
size_t ChemicalFormula::number( const size_t i ) const
{
    if ( sort_by_atomic_number_ )
    {
        std::map< Element, size_t >::const_iterator pos = elements_.begin();
        size_t j( 0 );
        while ( j != i )
        {
            ++j;
            ++pos;
        }
        return pos->second;
    }
    else
    {
        std::vector< Element > elements;
        for ( std::map< Element, size_t >::const_iterator pos( elements_.begin() ); pos != elements_.end(); ++pos )
        {
            elements.push_back( pos->first );
        }
        std::sort( elements.begin(), elements.end(), elements_less );
        for ( std::map< Element, size_t >::const_iterator pos( elements_.begin() ); pos != elements_.end(); ++pos )
        {
            if ( pos->first == elements[i] )
                return pos->second;
        }
    }
    throw std::runtime_error( "ChemicalFormula::number( size_t ): ." );
}

// ********************************************************************************

bool ChemicalFormula::contains( const Element element ) const
{
    return ( elements_.find( element ) != elements_.end() );
}

// ********************************************************************************

double ChemicalFormula::molecular_weight() const
{
    double result( 0.0 );
    for ( std::map< Element, size_t >::const_iterator pos( elements_.begin() ); pos != elements_.end(); ++pos )
        result += pos->first.atomic_weight() * pos->second;
    return result;
}

// ********************************************************************************

double ChemicalFormula::solid_state_volume() const
{
    double result( 0.0 );
    for ( std::map< Element, size_t >::const_iterator pos( elements_.begin() ); pos != elements_.end(); ++pos )
        result += pos->first.solid_state_volume() * pos->second;
    return result;
}

// ********************************************************************************

size_t ChemicalFormula::nelectrons() const
{
    size_t result( 0 );
    for ( std::map< Element, size_t >::const_iterator pos( elements_.begin() ); pos != elements_.end(); ++pos )
        result += pos->first.atomic_number() * pos->second;
    return result;
}

// ********************************************************************************

std::string ChemicalFormula::to_string( const bool insert_spaces, const bool explicit_1 ) const
{
    std::string result;
    for ( size_t i( 0 ); i != size(); ++i )
    {
        if ( insert_spaces && ( i != 0 ) )
            result += " ";
        result += element( i ).symbol();
        if ( ( number( i ) != 1 ) || explicit_1 )
            result += size_t2string( number( i ) );
    }
    return result;
}

// ********************************************************************************

