#ifndef ANGLE_H
#define ANGLE_H

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

#include <cmath>
#include <iosfwd>

class Angle;

inline Angle absolute( const Angle angle );

const double CONSTANT_PI = 4.0*atan(1.0);

const double degrees2radians = CONSTANT_PI / 180.0;
const double radians2degrees = 180.0 / CONSTANT_PI;

/*
  A lightweight angle class, which hides the distinction between radians and degrees.
  
  An angle can only be constructed if it is specified whether the value was intended to be degrees or radians.
  
  The angle is *not* normalised to the interval [0,360>, because this class is meant to have as
  little overhead as possible with respect to using a plain double, and to allow code like:

    Angle sum;
    sum += Angle::from_radians( 50.0 );
    sum += Angle::from_radians( 60.0 );
    sum += Angle::from_radians( 70.0 );
    Angle average = sum / 3.0;
*/
class Angle
{
public:

    enum AngleType { RADIANS, DEGREES };

    // Default constructor: 0.0 radians.
    Angle(): angle_(0.0) {}

    // 
    Angle( const double value, const AngleType radians_or_degrees )
    {
        angle_ = (radians_or_degrees == RADIANS) ? value : value * degrees2radians;
    }

    // Named constructors, the names can be improved.
    static Angle from_radians( const double value ) { return Angle( value ); }
    static Angle from_degrees( const double value ) { return Angle( value * degrees2radians ); }

    static Angle angle_30_degrees()  { return Angle( CONSTANT_PI/6.0 ); }
    static Angle angle_45_degrees()  { return Angle( CONSTANT_PI/4.0 ); }
    static Angle angle_60_degrees()  { return Angle( CONSTANT_PI/3.0 ); }
    static Angle angle_90_degrees()  { return Angle( 0.5*CONSTANT_PI ); }
    static Angle angle_120_degrees() { return Angle( (2.0/3.0)*CONSTANT_PI ); }
    static Angle angle_180_degrees() { return Angle( CONSTANT_PI ); }
    static Angle angle_360_degrees() { return Angle( 2.0*CONSTANT_PI ); }

    double value_in_radians() const { return angle_; }
    double value_in_degrees() const { return angle_ * radians2degrees; }

    bool nearly_zero( const Angle tolerance = Angle::from_radians( 0.000001 ) ) const { return ::absolute( *this ) < tolerance; }
    bool nearly_90( const Angle tolerance = Angle::from_radians( 0.000001 ) ) const { return ::absolute( *this - angle_90_degrees() ) < tolerance; }

    Angle operator+( const Angle rhs ) const { return Angle(*this) += rhs; }
    Angle operator-( const Angle rhs ) const { return Angle(*this) -= rhs; }
    Angle operator*( const double rhs ) const { return Angle(*this) *= rhs; }
    Angle operator/( const double rhs ) const { return Angle(*this) /= rhs; }
    double operator/( const Angle rhs ) const { return angle_ / rhs.angle_; }

    Angle operator+() const { return Angle( *this ); }

    Angle operator-() const { return Angle( -angle_ ); }

    Angle & operator+=( const Angle rhs )
    {
        angle_ += rhs.angle_;
        return *this;
    }

    Angle & operator-=( const Angle rhs )
    {
        angle_ -= rhs.angle_;
        return *this;
    }

    Angle & operator*=( const double rhs )
    {
        angle_ *= rhs;
        return *this;
    }

    Angle & operator/=( const double rhs )
    {
        angle_ /= rhs;
        return *this;
    }

    // Of course, operator== for floating point values is nonsense.
    bool operator==( const Angle rhs ) const { return ( this->angle_ == rhs.angle_ ); }
    bool operator!=( const Angle rhs ) const { return ! ( *this == rhs ); }
    bool operator< ( const Angle rhs ) const { return ( this->angle_ < rhs.angle_ ); }
    bool operator> ( const Angle rhs ) const { return ( rhs < *this ); }
    bool operator>=( const Angle rhs ) const { return ! ( *this < rhs ); }
    bool operator<=( const Angle rhs ) const { return ! ( rhs < *this ); }

    double sine() const { return sin( angle_ ); }
    double cosine() const { return cos( angle_ ); }
    void sincos( double & sine, double & cosine ) const { sine = this->sine(); cosine = this->cosine(); }
    double tangent() const { return tan( angle_ ); }
    double cotangent() const { return 1.0/this->tangent(); }
    
    void absolute() { if ( angle_ < 0.0 ) angle_ = -angle_; }

    friend Angle absolute( const Angle angle );

private:
    double angle_; // The angle, stored in radians.
    // The angle is stored in radians because radians are used for calculations,
    // i.e. where speed is important, and degrees are only necessary when user output
    // is required, which takes a lot of time anyway, so the time for the conversion
    // is negligible.
    
    // Private constructor for named constructors.
    explicit Angle( const double value ): angle_(value) {}
};

inline Angle absolute( const Angle angle ) { return ( angle.angle_ < 0.0 ) ? Angle( -angle.angle_ ) : angle; }

std::ostream & operator<<( std::ostream & os, const Angle angle );

inline Angle arcsine( const double value ) { return Angle::from_radians( asin( value ) ); }
inline Angle arccosine( const double value ) { return Angle::from_radians( acos( value ) ); }

// It is more efficient to calculate sine and cosine of the same angle simultaneously.
// The algorithm that is used is an APPROXIMATION.
void sincos( Angle angle, double & sine, double & cosine );

inline Angle arctangent( const double x ) { return Angle::from_radians( atan( x ) ); }

// cot-1(x) = pi/2 - tan-1(x) for any x.
inline Angle arccotangent( const double x ) { return Angle::angle_90_degrees() - arctangent( x ); }

inline Angle operator*( const double lhs, const Angle rhs ) { return Angle::from_radians( rhs.value_in_radians() * lhs ); }

inline Angle ATAN2( const double y, const double x ) { return Angle::from_radians( atan2( y, x ) ); }

// To calculate the average of four values:
// Angle average = average( value_1, value_2 );
// average = average( value_3, average, 2.0 );
// average = average( value_4, average, 3.0 );
// Even this works:
//    Angle prev_estimate; // No need to initialise...
//    Angle next_estimate; // No need to initialise...
//    size_t iStep( 0 );
//    while ( iStep < 1000000 )
//    {
//        prev_estimate = next_estimate;
//        Angle current_value = some_function();
//        next_estimate = average( current_value, prev_estimate, iStep );
//        ++iStep;
//    }
inline Angle average( const Angle lhs, const Angle rhs, const double weight = 1.0 )
{
    return ( lhs + weight * rhs ) / ( 1.0 + weight );
}

inline bool nearly_equal( const Angle lhs, const Angle rhs, const Angle tolerance = Angle::from_radians( 0.000001 ) )
{
    return ( absolute( rhs - lhs ) < tolerance );
}

bool triquality( const Angle x1, const Angle x2, const Angle x3, const Angle tolerance = Angle::from_radians( 0.000001 ) );

#endif // ANGLE_H

