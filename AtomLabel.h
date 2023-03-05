#ifndef ATOMLABEL_H
#define ATOMLABEL_H

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

#include <cstddef> // For definition of size_t.
#include <string>

// Assumes that an atom label is e.g. C12b_5_0_1_0
// This is element, index, "a", "A", "b" or "B", symmetry operator index, u, v, w (zero based)
class AtomLabel
{
public:

    AtomLabel( const std::string & label );

    std::string label() const { return label_; };

    std::string element() const { return element_; };

    size_t index() const { return index_; };

    std::string subindex() const { return subindex_; };

    size_t symmetry_operator_index() const { return symmetry_operator_index_; };

    size_t u() const { return u_; };

    size_t v() const { return v_; };

    size_t w() const { return w_; };

private:
    std::string label_;
    std::string element_;
    size_t index_;
    std::string subindex_;
    size_t symmetry_operator_index_;
    size_t u_;
    size_t v_;
    size_t w_;
};

#endif // ATOMLABEL_H

