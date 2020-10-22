#ifndef POWDERMATCHTABLE_H
#define POWDERMATCHTABLE_H

/* *********************************************
Copyright (c) 2013-2020, Cornelis Jan (Jacco) van de Streek
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

#include "CrystalLattice.h"
//#include "Fraction.h"
#include "MillerIndices.h"

#include <string>
#include <vector>

/*

*/
class PowderMatchTable
{
public:

    // Default constructor
    PowderMatchTable();

    PowderMatchTable( const FileName & file_name );
    
    std::string name( const size_t i ) const { return names_[i]; }
    double density( const size_t i ) const { return densities_[i]; }
    std::string space_group_name( const size_t i ) const { return space_group_names_[i]; }
    double figure_of_merit( const size_t i ) const { return figures_of_merit_[i]; }
    double Biso( const size_t i ) const { return Biso_values_[i]; }
    double MarchDollase( const size_t i ) const { return MarchDollase_values_[i]; }
    MillerIndices PO_direction( const size_t i ) const { return PO_directions_[i]; }
    double cell_deformation( const size_t i ) const { return cell_deformations_[i]; }
    CrystalLattice crystal_lattice( const size_t i ) const { return crystal_lattices_[i]; }

private:
    std::vector< std::string > names_;
//    std::vector<  > status
//    std::vector<  > comment
//    std::vector<  > energy
    std::vector< double > densities_;
    std::vector< std::string > space_group_names_;
//    std::vector< Fraction > number_of_dof
    std::vector< double > figures_of_merit_;
    std::vector< double > Biso_values_;
    std::vector< double > MarchDollase_values_;
    std::vector< MillerIndices > PO_directions_;
    std::vector< double > cell_deformations_;
    std::vector< CrystalLattice > crystal_lattices_;
};

#endif // POWDERMATCHTABLE_H

