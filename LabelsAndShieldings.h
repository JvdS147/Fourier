#ifndef LABELSANDSHIELDINGS_H
#define LABELSANDSHIELDINGS_H

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

class FileName;

#include "Element.h"
#include "Mapping.h"

#include <vector>
#include <string>

class LabelsAndShieldings
{
public:

    LabelsAndShieldings() {}

    // Invalidates the sorting of the list
    void push_back( const std::string & label, const double shielding );

    size_t size() const { return labels_.size(); }

    // The index is zero-based
    // We don't actually sort the lists, but create a sorted map
    std::string label(     const size_t i ) const { return labels_[ sorted_map_[ i ] ]; }
    double      shielding( const size_t i ) const { return shieldings_[ sorted_map_[ i ] ]; }

    // For debugging
    void show() const;

    void save( const FileName & file_name ) const;

private:
    std::vector< std::string > labels_;
    std::vector< Element > elements_;
    std::vector< double > shieldings_;
    // We don't actually sort the lists, but create a sorted map
    mutable Mapping sorted_map_;

    // We don't actually sort the lists, but create a sorted map
    // We sort first by element, then by shieldings in descending order
    void sort() const;
};

#endif // LABELSANDSHIELDINGS_H

