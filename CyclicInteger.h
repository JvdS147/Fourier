#ifndef CYCLICINTEGER_H
#define CYCLICINTEGER_H

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

#include <limits>

/*
  A cyclic integer. It has a start (usually 0) and an end, and the value of the class is always between these two (inclusive).
  
  Remarkably little protection against underflow or overflow of intermediate results...

  @@@@@@@@@@@@@ what if start=end????
  
*/
class CyclicInteger
{
public:

    // Default constructor: [INT_MIN,INT_MAX], value = 0.
    CyclicInteger();

    // Constructor to ensure behaviour mimics an int.
    // Could decide not to make this one explicit.
    explicit CyclicInteger( const int value );

    // Constructor
    CyclicInteger( const int start, const int end, const int value );

    // Does NOT increment value by 1.
    int current_value() const;

    // Increments value by 1.
    int next_value() const;

    // Returns current value + n, adjusted for cyclicity, leaves current value unchanged
    int plus_n( const int n ) const;

    // Implicit conversion
//    operator int() const { return value_; }

    CyclicInteger operator+( const CyclicInteger rhs ) const { return CyclicInteger(*this) += rhs; }
    CyclicInteger operator-( const CyclicInteger rhs ) const { return CyclicInteger(*this) -= rhs; }
    CyclicInteger operator*( const int rhs ) const { return CyclicInteger(*this) *= rhs; }
    CyclicInteger operator/( const int rhs ) const { return CyclicInteger(*this) /= rhs; }
//    double operator/( const CyclicInteger rhs ) const { return angle_ / rhs.value_; }

    CyclicInteger operator+() const { return CyclicInteger( *this ); }

    CyclicInteger operator-() const { return CyclicInteger( -value_ ); }

    CyclicInteger & operator+=( const CyclicInteger rhs )
    {
        value_ += rhs.value_;
        adjust();
        return *this;
    }

    CyclicInteger & operator-=( const CyclicInteger rhs )
    {
        value_ -= rhs.value_;
        adjust();
        return *this;
    }

    CyclicInteger & operator*=( const int rhs )
    {
        value_ *= rhs;
        adjust();
        return *this;
    }

    CyclicInteger & operator/=( const int rhs )
    {
        value_ /= rhs;
        adjust();
        return *this;
    }

    bool operator==( const CyclicInteger rhs ) const { return ( this->value_ == rhs.value_ ); }
    bool operator!=( const CyclicInteger rhs ) const { return ! ( *this == rhs ); }
    bool operator< ( const CyclicInteger rhs ) const { return ( this->value_ < rhs.value_ ); }
    bool operator> ( const CyclicInteger rhs ) const { return ( rhs < *this ); }
    bool operator>=( const CyclicInteger rhs ) const { return ! ( *this < rhs ); }
    bool operator<=( const CyclicInteger rhs ) const { return ! ( rhs < *this ); }


    CyclicInteger & operator++();    // Prefix
    CyclicInteger   operator++(int); // Postfix
    CyclicInteger & operator--();    // Prefix
    CyclicInteger   operator--(int); // Postfix

//    void absolute() { angle_ = std::abs( angle_ ); }

    friend CyclicInteger absolute( const CyclicInteger cyclic_integer );
//    friend CyclicInteger square( const CyclicInteger cyclic_integer );

private:

//    inline void negate() { integer_part_ = -integer_part_; numerator_ = -numerator_; }

    int offset_;
    unsigned int range_;
    mutable int value_;

    void initialise( const int start, const int end );
    void adjust() const;
    int adjust( const int value ) const;

};

inline bool operator==( const CyclicInteger lhs, const int rhs ) { return ( lhs.current_value() == rhs ); }
inline bool operator!=( const CyclicInteger lhs, const int rhs ) { return ! ( lhs.current_value() == rhs ); }
inline bool operator< ( const CyclicInteger lhs, const int rhs ) { return ( lhs.current_value() < rhs ); }
inline bool operator< ( const int lhs, const CyclicInteger rhs ) { return ( lhs < rhs.current_value() ); }
inline bool operator> ( const CyclicInteger lhs, const int rhs ) { return ( rhs < lhs ); }
inline bool operator> ( const int lhs, const CyclicInteger rhs ) { return ( rhs < lhs ); }
inline bool operator>=( const CyclicInteger lhs, const int rhs ) { return ! ( lhs.current_value() < rhs ); }
inline bool operator<=( const CyclicInteger lhs, const int rhs ) { return ! ( rhs < lhs.current_value() ); }
inline bool operator==( const int lhs, const CyclicInteger rhs ) { return ( lhs == rhs.current_value() ); }
inline bool operator!=( const int lhs, const CyclicInteger rhs ) { return ! ( lhs == rhs ); }
inline bool operator>=( const int lhs, const CyclicInteger rhs ) { return ! ( lhs < rhs ); }
inline bool operator<=( const int lhs, const CyclicInteger rhs ) { return ! ( rhs < lhs ); }

//inline CyclicInteger absolute( const CyclicInteger cyclic_integer ) { return CyclicInteger( cyclic_integer.start_, cyclic_integer.end_, std::abs( cyclic_integer.value_ ) ); }

//std::ostream & operator<<( std::ostream & os, const CyclicInteger cyclic_integer );

//inline CyclicInteger operator*( const double lhs, const CyclicInteger rhs ) { return CyclicInteger::from_radians( rhs.value_in_radians() * lhs ); }

//inline CyclicInteger square( const CyclicInteger cyclic_integer )
//{
//    return CyclicInteger( cyclic_integer.start_, cyclic_integer.end_, cyclic_integer.value_ * cyclic_integer.value_ );
//}

#endif // CYCLICINTEGER_H

