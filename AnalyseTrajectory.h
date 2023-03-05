#ifndef ANALYSETRAJECTORY_H
#define ANALYSETRAJECTORY_H

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

#include "CrystalLattice.h"
#include "FileList.h"
#include "RunningAverageAndESD.h"
#include "SpaceGroup.h"
#include "Vector3D.h"

#include <vector>

/*
  Analyses an MD trajectory which must be provided as a set of cif files in the correct order.
  Each cif file must have the same number of atoms, in the same order.
  u, v, w are the dimensions of the supercell with respect to the original unit cell.
  Collapse supercell, assume order *in the unit cell* (not in the molecule) can be trusted
  (if there are n atoms in a unit cell, then atom n+1 corresponds to atom 1 in unit cell 1).
  The method assumes that you have not repositioned atoms to lie within the unit cell.
  The method assumes that we are dealing with a solid where atomic coordinates are fairly constant and symmetry-related copies can be
  easily identified from simple geometric considerations.
  Because it is assumed to be a solid, there is no provision for rotational drift.
  Space-group symmetry is also used to group *all* symmetry copies of each atom--so make sure that the space group is correct,
  i.e. that no phase transition has taken place. Of course, if the space group is P1 there is no problem.
  You'll get horrible results when applying the space-group symmetry if you have manually repositioned the molecules so that all molecules were comfortably within the unit cell.
  Based on experience, Z' = 2 x 1/2 does not preserve the order of the atoms, so you cannot specify a space group in that case.
*/
class AnalyseTrajectory
{
public:

    // There is no default constructor
    //AnalyseTrajectory();

    // Make sure the space group name is set properly: it is written to the cif file.
    // transformation does not work
    explicit AnalyseTrajectory( const FileList file_list,
                                const size_t u = 1,
                                const size_t v = 1,
                                const size_t w = 1,
                                const SpaceGroup & space_group = SpaceGroup(),
                                const Matrix3D & transformation = Matrix3D() );

    enum DriftCorrection { NONE, USE_FIRST_FRAME, USE_VECTOR };

//    void set_u_v_w( const size_t u, const size_t v, const size_t w ) { u_ = u; v_ = v; w_ = w; }

//    void set_space_group( const SpaceGroup & space_group ) { space_group_ = space_group; }

//    void set_drift_correction( const DriftCorrection drift_correction ) { drift_correction_ = drift_correction; }
    
//    void set_drift_correction_vector( const Vector3D & drift_correction_vector ) { drift_correction_vector_ = drift_correction_vector; }

    // Define a unit-cell transformation, atomic coordinates are also transformed
//    void set_transformation( const Matrix3D & transformation ) { transformation_ = transformation; apply_transformation_ = true; }

    // Not yet implemented. Applied *after* drift correction
//    void set_translation( const Vector3D & translation ) { translation_ = translation; }

    CrystalLattice average_crystal_lattice() const;

    std::vector< Vector3D > centres_of_mass() const { return centres_of_mass_; }
    
    // Convenience function for lazy people. Writes centres of mass to file in same directory as FileList.
    void save_centres_of_mass() const;
    
//  ADPs / ESDs / averages

private:
    FileList file_list_;
    size_t u_;
    size_t v_;
    size_t w_;
    SpaceGroup space_group_;
    Matrix3D transformation_;
    Vector3D translation_;
    bool write_lean_;
    bool write_average_;
    bool write_average_noH_;
    bool write_average_ESDs_;
    bool write_sum_;
    DriftCorrection drift_correction_;
    Vector3D drift_correction_vector_;
    RunningAverageAndESD<double> average_a_;
    RunningAverageAndESD<double> average_b_;
    RunningAverageAndESD<double> average_c_;
    RunningAverageAndESD<Angle> average_alpha_;
    RunningAverageAndESD<Angle> average_beta_;
    RunningAverageAndESD<Angle> average_gamma_;
    RunningAverageAndESD<double> average_volume_;
    std::vector< Vector3D > centres_of_mass_; // Monitors the drift.

    void analyse();
};

#endif // ANALYSETRAJECTORY_H

