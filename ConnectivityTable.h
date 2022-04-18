#ifndef CONNECTIVITYTABLE_H
#define CONNECTIVITYTABLE_H

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
    * Neither the name of the University of Copenhagen nor the
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

#include <cstddef> // For definition of size_t
#include <vector>

/*
  A connectivity table in its simplest form: for N atoms a NxN table.
  
  Could make the entry the type of bond (single, double), but currently it is just "bonded or not bonded".
  
  We will probably have things like Molecule3D, Molecule3D and ChemicalCompound (= Molecule2D objects + stoichiometry + racemic yes / no)
  in the future.
*/
class ConnectivityTable
{
public:

    explicit ConnectivityTable( const size_t natoms = 0 );

    size_t size() const { return dimension_; }

    // In keeping with the silly C++ convention: zero-based
    // Returns 1 for i == j
    size_t value( size_t i, size_t j ) const;

    void set_value( size_t i, size_t j, const size_t value );
    
    void show() const;

private:
    std::vector< size_t > data_; // Could have used bool, but I assume it will store the bond type one day.
    size_t dimension_;

};

std::vector< std::vector< size_t > > split( const ConnectivityTable & connectivity_table );

#endif // CONNECTIVITYTABLE_H

