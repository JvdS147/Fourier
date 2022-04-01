#ifndef REFLECTIONLIST_H
#define REFLECTIONLIST_H

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

class FileName;

#include "Mapping.h"
#include "MillerIndices.h"

#include <vector>

// We deliberately store d-spacing and not 2 theta, because d-spacing is not wavelength dependent.
// We store F^2 and not "intensity", because intensity also includes contributions from the
// Debije-Waller factor and the LP-factor.
class ReflectionList
{
public:

    ReflectionList();

    void push_back( const MillerIndices & miller_indices, const double F_squared, const double d_spacing, const size_t multiplicity );

    void reserve( const size_t nvalues );
    size_t size() const { return miller_indices_.size(); }

    // Does not understand anything about equivalence, so (100) is not the same as (-100).
    // Returns size() if not found.
    size_t index( const MillerIndices & miller_indices );

    // The index is zero-based
    // We don't actually sort the lists, but create a sorted map
    MillerIndices miller_indices( const size_t i ) const { return miller_indices_[ sorted_map_[ i ] ];   }
    double        F_squared(      const size_t i ) const { return F_squared_[ sorted_map_[ i ] ];      }
    double        d_spacing(      const size_t i ) const { return d_spacings_[ sorted_map_[ i ] ]; }
    size_t        multiplicity(   const size_t i ) const { return multiplicity_[ sorted_map_[ i ] ]; }

    void set_miller_indices( const size_t i, const MillerIndices & miller_indices ) { miller_indices_[ sorted_map_[ i ] ] = miller_indices; }
    void set_F_squared(      const size_t i, const double F_squared ) { F_squared_[ sorted_map_[ i ] ] = F_squared; }
    void set_d_spacing(      const size_t i, const double d_spacing ) { d_spacings_[ sorted_map_[ i ] ] = d_spacing; sort_by_d_spacing(); }
    void set_multiplicity(   const size_t i, const size_t multiplicity ) { multiplicity_[ sorted_map_[ i ] ] = multiplicity; }

    // A ReflectionsList and a SHELX .hkl file are quite different, so this is a bit of an abuse of the class...
    void read_hkl( const FileName & file_name );

    // For debugging
    void show() const;

    void save( const FileName & file_name ) const;

private:
    std::vector< MillerIndices > miller_indices_;
    std::vector< double >        F_squared_;
    std::vector< double >        d_spacings_;
    std::vector< size_t >        multiplicity_;
    // We don't actually sort the lists, but create a sorted map
    Mapping sorted_map_;

    // We don't actually sort the lists, but create a sorted map
    void sort_by_d_spacing();
};

#endif // REFLECTIONLIST_H

