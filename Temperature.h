#ifndef TEMPERATURE_H
#define TEMPERATURE_H

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

#include <cmath>
#include <iosfwd>

/*

*/
class Temperature
{
public:

    enum TemperatureUnit { KELVIN, CELSIUS };

    // Default constructor: 0.0 GPa
    Temperature(): temperature_(0.0) {}

    Temperature( const double value, const TemperatureUnit temperature_unit );

    // Named constructors
    static Temperature from_Kelvin( const double value ) { return Temperature( value ); }
    static Temperature from_Celsius( const double value ) { return Temperature( value + 273.15 ); }

    static Temperature ambient_temperature() { return from_Celsius( 20.0 ); }

    double value_in_Kelvin() const { return temperature_; }
    double value_in_Celsius() const { return temperature_ - 273.15; }

    Temperature operator+( const Temperature rhs ) const { return Temperature(*this) += rhs; }
    Temperature operator-( const Temperature rhs ) const { return Temperature(*this) -= rhs; }
    Temperature operator*( const double rhs ) const { return Temperature(*this) *= rhs; }
    Temperature operator/( const double rhs ) const { return Temperature(*this) /= rhs; }
    double operator/( const Temperature rhs ) const { return temperature_ / rhs.temperature_; }

    Temperature operator+() const { return Temperature( *this ); }

    Temperature operator-() const { return Temperature( -temperature_ ); }

    Temperature & operator+=( const Temperature rhs )
    {
        temperature_ += rhs.temperature_;
        return *this;
    }

    Temperature & operator-=( const Temperature rhs )
    {
        temperature_ -= rhs.temperature_;
        return *this;
    }

    Temperature & operator*=( const double rhs )
    {
        temperature_ *= rhs;
        return *this;
    }

    Temperature & operator/=( const double rhs )
    {
        temperature_ /= rhs;
        return *this;
    }

    // Of course, operator== for floating point values is nonsense
    bool operator==( const Temperature rhs ) const { return ( this->temperature_ == rhs.temperature_ ); }
    bool operator!=( const Temperature rhs ) const { return ! ( *this == rhs ); }
    bool operator< ( const Temperature rhs ) const { return ( this->temperature_ < rhs.temperature_ ); }
    bool operator> ( const Temperature rhs ) const { return ( rhs < *this ); }
    bool operator>=( const Temperature rhs ) const { return ! ( *this < rhs ); }
    bool operator<=( const Temperature rhs ) const { return ! ( rhs < *this ); }

    friend Temperature absolute( const Temperature temperature );

private:
    double temperature_; // The temperature, stored in GigaPascal

    // Private constructor for named constructors
    explicit Temperature( const double value ) : temperature_(value) {}

};

std::ostream & operator<<( std::ostream & os, const Temperature temperature );

inline Temperature absolute( const Temperature temperature ) { return Temperature( std::abs( temperature.temperature_ ) ); }

inline Temperature operator*( const double lhs, const Temperature rhs ) { return Temperature::from_Kelvin( rhs.value_in_Kelvin() * lhs ); }

inline bool nearly_equal( const Temperature lhs, const Temperature rhs, const Temperature tolerance = Temperature::from_Kelvin( 0.0000001 ) )
{
    return ( absolute( rhs - lhs ) < tolerance );
}

#endif // TEMPERATURE_H

