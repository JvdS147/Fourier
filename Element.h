#ifndef ELEMENT_H
#define ELEMENT_H

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

#include <string>
#include <cstddef> // For definition of size_t

// H and D are two distinct elements
class Element
{
public:

    Element() : id_(6) {}

    // atomic_number must be smaller than 113.
    // atomic number 1 defaults to hydrogen, deuterium cannot be constructed with this constructor,
    // use the std::string based constructor instead.
    explicit Element( const size_t atomic_number );

    explicit Element( std::string symbol );

    // D = 0, others are atomic number.
    size_t id() const { return id_; }

    size_t atomic_number() const;

    std::string symbol() const;

    double Van_der_Waals_radius() const;

    bool is_H_or_D() const { return ( id_ == 0 ) || ( id_ == 1 ); }

    double atomic_weight() const;

    // Uses Hofmann's values.
    double solid_state_volume() const;

    double scattering_factor( const double sine_theta_over_lambda ) const;

    bool operator<( const Element & rhs ) const { return ( id_ < rhs.id_ ); }
    
private:
    // We really need to be able to distinguish between D and H,
    // so just the atomic number is not enough
    size_t id_;

};

// The ususal order of elements: C, H, D, rest alphabetical.
bool elements_less( const Element & lhs, const Element & rhs );

// The distance is the distance squared, to avoid the necessity of an expensive sqrt().
bool are_bonded( const Element & lhs, const Element & rhs, const double distance2 );

// Helper function, assumes atom label of the type "C3  ".
// "C_3" and "Ow14" are also interpreted correctly.
Element element_from_atom_label( std::string label );

// H and D test as different.
inline bool operator==( const Element & lhs, const Element & rhs ) { return ( lhs.id() == rhs.id() ); }

// H and D test as different.
inline bool operator!=( const Element & lhs, const Element & rhs ) { return ! ( lhs == rhs ); }

#endif // ELEMENT_H

