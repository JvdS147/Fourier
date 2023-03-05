#ifndef PRESSURE_H
#define PRESSURE_H

/* *********************************************
Copyright (c) 2013-2023, Cornelis Jan (Jacco) van de Streek
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
class Pressure
{
public:

    enum PressureUnit { GPa, MPa, bar, atm };

    // Default constructor: 0.0 GPa
    Pressure(): pressure_(0.0) {}

    Pressure( const double value, const PressureUnit pressure_unit );

    // Named constructors
    static Pressure from_GPa( const double value ) { return Pressure( value ); }
    static Pressure from_MPa( const double value ) { return Pressure( value / 1000.0 ); }
    static Pressure from_kPa( const double value ) { return Pressure( value / 1000000.0 ); }
    static Pressure from_Pa(  const double value ) { return Pressure( value / 1000000000.0 ); }

    static Pressure ambient_pressure() { return Pressure( 101.325 / 1000000.0 ); }

    double value_in_GPa() const { return pressure_; }
    double value_in_MPa() const { return pressure_ * 1000.0; }
    double value_in_kPa() const { return pressure_ * 1000000.0; }
    double value_in_Pa() const { return pressure_ * 1000000000.0; }

    Pressure operator+( const Pressure rhs ) const { return Pressure(*this) += rhs; }
    Pressure operator-( const Pressure rhs ) const { return Pressure(*this) -= rhs; }
    Pressure operator*( const double rhs ) const { return Pressure(*this) *= rhs; }
    Pressure operator/( const double rhs ) const { return Pressure(*this) /= rhs; }
    double operator/( const Pressure rhs ) const { return pressure_ / rhs.pressure_; }

    Pressure operator+() const { return Pressure( *this ); }

    Pressure operator-() const { return Pressure( -pressure_ ); }

    Pressure & operator+=( const Pressure rhs )
    {
        pressure_ += rhs.pressure_;
        return *this;
    }

    Pressure & operator-=( const Pressure rhs )
    {
        pressure_ -= rhs.pressure_;
        return *this;
    }

    Pressure & operator*=( const double rhs )
    {
        pressure_ *= rhs;
        return *this;
    }

    Pressure & operator/=( const double rhs )
    {
        pressure_ /= rhs;
        return *this;
    }

    // Of course, operator== for floating point values is nonsense
    bool operator==( const Pressure rhs ) const { return ( this->pressure_ == rhs.pressure_ ); }
    bool operator!=( const Pressure rhs ) const { return ! ( *this == rhs ); }
    bool operator< ( const Pressure rhs ) const { return ( this->pressure_ < rhs.pressure_ ); }
    bool operator> ( const Pressure rhs ) const { return ( rhs < *this ); }
    bool operator>=( const Pressure rhs ) const { return ! ( *this < rhs ); }
    bool operator<=( const Pressure rhs ) const { return ! ( rhs < *this ); }

    friend Pressure absolute( const Pressure pressure );

private:
    double pressure_; // The pressure, stored in GigaPascal

    // Private constructor for named constructors
    explicit Pressure( const double value ) : pressure_(value) {}

};

std::ostream & operator<<( std::ostream & os, const Pressure pressure );

inline Pressure absolute( const Pressure pressure ) { return Pressure( std::abs( pressure.pressure_ ) ); }

inline Pressure operator*( const double lhs, const Pressure rhs ) { return Pressure::from_GPa( rhs.value_in_GPa() * lhs ); }

inline bool nearly_equal( const Pressure lhs, const Pressure rhs, const Pressure tolerance = Pressure::from_GPa( 0.0000001 ) )
{
    return ( absolute( rhs - lhs ) < tolerance );
}

#endif // PRESSURE_H

