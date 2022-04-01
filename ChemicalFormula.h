#ifndef CHEMICALFORMULA_H
#define CHEMICALFORMULA_H

/* *********************************************
Copyright (c) 2013-2022, Cornelis Jan (Jacco) van de Streek
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

class Fraction;

#include "Element.h"

#include <map>
#include <string>

/*
*/
class ChemicalFormula
{
public:

    ChemicalFormula() : sort_by_atomic_number_(false) {}
    
    // There will always be ambiguity, e.g. "BI" can mean Bismut or Boron and Iodine.
    // A solution would be to insist on writing "BI" as "B1I1", but that is unusal
    // and the whole point of a constructor is to do the work for you.
    // The chemical formula must be well formed, so HCl *must* be HCl, not HCL. C10C10H10 is not possible either.
    explicit ChemicalFormula( const std::string & input );

    void add_element( const Element element );
    
    // Elements can be sorted by atomic number or in the "usual" order: C, H, D, rest alphabetical.
    // The default is: C, H, D, rest alphabetical.
    bool sort_by_atomic_number() const { return sort_by_atomic_number_; }
    
    void set_sort_by_atomic_number( const bool sort_by_atomic_number ) { sort_by_atomic_number_ = sort_by_atomic_number;}

    // Multiplies the chemical formula by an integer, e.g. to convert from an asymmetric unit to a whole molecule
    // This probably should not be a member function
    void multiply( const Fraction fraction );

    size_t size() const { return elements_.size(); }
    
    // See sort_by_atomic_number()
    // Currently extremely wasteful but works and is safe
    Element element( const size_t i ) const;

    // See sort_by_atomic_number()
    // Currently extremely wasteful but works and is safe
    size_t number( const size_t i ) const;

    bool contains( const Element element ) const;
    
    double molecular_weight() const;

    // Uses Hofmann's values
    double solid_state_volume() const;

    size_t nelectrons() const;

    std::string to_string( const bool insert_spaces = false, const bool explicit_1 = false ) const;

private:
    std::map< Element, size_t > elements_;
    bool sort_by_atomic_number_;
};

#endif // CHEMICALFORMULA_H

