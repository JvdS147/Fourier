#ifndef VOIDSFINDER_H
#define VOIDSFINDER_H

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

class CrystalStructure;
/*

*/

// Returns the void volume
// Uses the old algorithm, based on first finding the set of spheres that determine the voids followed by
// determining the volume of the set of spheres.
double find_voids( const CrystalStructure & crystal_structure, const double probe_radius = 1.2 );

// DO NOT USE
// Currently only returns the percentage voids
// Uses the new algorithm that samples the unit cell and for each sample point determines if it is
// part of a void or not.
// This gives the wrong answer, we only find the points where the centre of a sphere would fit, or we find *all* points outside all molecules,
// including the tiny voids where no solvent can ever reach.
double find_voids_2( const CrystalStructure & crystal_structure );

// Probe size is 0.0, this is useful for calculating e.g. the packing coefficient.
double void_volume( const CrystalStructure & crystal_structure );

#endif // VOIDSFINDER_H

