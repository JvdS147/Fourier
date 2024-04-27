#ifndef TALLY_H
#define TALLY_H

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

//#include <string>
//#include <iosfwd>

#include <map>
#include <iostream>

/*
 * Counts number of items.
 * 
 * Basically gives a C++ container an old-fashioned index-based interface, probably not good.
 */
template <class T>
class Tally
{
public:

    Tally() {}

    void add( const T & t, const size_t counts = 1 )
    {
        if ( contains( t ) )
            counts_[ t ] += counts;
        else
            counts_[ t ] = counts;
    }

    size_t size() const { return counts_.size(); }

    bool contains( const T & t )
    {
        return ( counts_.find( t ) != counts_.end() );
    }

    std::vector< T > values() const
    {
        std::vector< T > result;
        for ( typename std::map< T, size_t >::const_interator it( counts_.begin() ); it != counts_.end(); ++it )
        {
            result.push_back( it->first );
        }
        return result;
    }

//    T value( const size_t i ) const
//    {
//        
//    }

    size_t counts( const T & t ) const
    {
        return contains( t ) ? counts_[ t ] : 0;
    }

    void show() const
    {
        for ( typename std::map< T, size_t >::const_iterator it( counts_.begin() ); it != counts_.end(); ++it )
        {
            std::cout << it->first << " " << it->second << std::endl;
        }
    }

//    Fraction operator+( const Fraction & rhs ) const { return Fraction(*this) += rhs; }
//    Fraction operator-( const Fraction & rhs ) const { return Fraction(*this) -= rhs; }
//    Fraction operator*( const Fraction & rhs ) const { return Fraction(*this) *= rhs; }
//    Fraction operator/( const Fraction & rhs ) const { return Fraction(*this) /= rhs; }

//    Fraction operator+() const { return Fraction( *this ); }

//    Fraction & operator+=( const Fraction & rhs );
//    Fraction & operator-=( const Fraction & rhs );
//    Fraction & operator*=( const Fraction & rhs );
//    Fraction & operator/=( const Fraction & rhs );

//    Fraction & operator++();    // Prefix
//    Fraction   operator++(int); // Postfix
//    Fraction & operator--();    // Prefix
//    Fraction   operator--(int); // Postfix

//    bool operator==( const Fraction & rhs ) const;
//    bool operator!=( const Fraction & rhs ) const { return ! ( *this == rhs ); }
//    bool operator< ( const Fraction & rhs ) const;
//    bool operator> ( const Fraction & rhs ) const { return ( rhs < *this ); }
//    bool operator>=( const Fraction & rhs ) const { return ! ( *this < rhs ); }
//    bool operator<=( const Fraction & rhs ) const { return ! ( rhs < *this ); }

    // Returns the Fraction in string form, e.g. "0", "3", "1/2", "2 + 2/3", "-9 + -1/4"
//    std::string to_string() const;

private:

    std::map< T, size_t > counts_;
};

#endif // TALLY_H

