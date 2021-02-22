#ifndef ANALYSERINGS_H
#define ANALYSERINGS_H

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

//class CollectionOfPoints;
class Vector3D;

#include <cstddef> // For definition of size_t
#include <vector>

class FiveMemberedRingAnalyser
{
public:
    
    enum GeometryType { NONE, AXIAL, EQUATORIAL };

    FiveMemberedRingAnalyser();

    explicit FiveMemberedRingAnalyser( const std::vector< Vector3D > & points );

    void analyse( const std::vector< Vector3D > & points );
    
    std::vector< double > results() const { return results_; }
    bool is_planar() const { return is_planar_; }
    bool is_envelope() const { return is_envelope_; }
    bool is_double_envelope() const { return is_double_envelope_; }
    size_t unique_envelope_point_1() const { return unique_envelope_point_1_; }
    size_t unique_envelope_point_2() const { return unique_envelope_point_2_; }

    // Returns NONE if the ring is not a single envelope
    // This should be done in a more sophisticated manner. "AXIAL" is unambiguous,
    // but equatorial is harder to establish. So the proper way to do this would
    // be to supply *both* bonded atoms and to determine which is axial,
    // the other is then equatorial. For the moment, this has been solved
    // by returning "equatorial" if it is not axial.
    GeometryType axial_or_equatorial( const Vector3D & point ) const;

    std::vector< double > distances_from_plane_;
    std::vector< double > rmsds_from_mean_plane_;
    std::vector< size_t > sorted_map_;

private:
    std::vector< Vector3D > points_;
    bool is_planar_;
    bool is_envelope_;
    bool is_double_envelope_;
    double root_mean_square_devation_from_mean_plane_;
    double distance_;
    size_t unique_envelope_point_1_; // Index of the point that forms the "flap" of the envelope
    size_t unique_envelope_point_2_; // Index of the point that forms the "flap" of the envelope

    std::vector< double > results_;

};

// The difference between boat and twisted boat is a gliding scale.
class SixMemberedRingAnalyser
{
public:
    
    enum GeometryType { NONE, AXIAL, EQUATORIAL };

    SixMemberedRingAnalyser();

    explicit SixMemberedRingAnalyser( const std::vector< Vector3D > & points );

    void analyse( const std::vector< Vector3D > & points );

    bool is_planar() const { return is_planar_; }
    bool is_chair() const { return is_chair_; }
    
    // The difference between boat and twisted boat is a gliding scale.
    bool is_boat() const { return is_boat_; }
    bool is_twisted_boat() const { return is_twisted_boat_; }

    // Returns NONE if the ring is not a chair
    // This should be done in a more sophisticated manner. "AXIAL" is unambiguous,
    // but equatorial is harder to establish. So the proper way to do this would
    // be to supply *both* bonded atoms and to determine which is axial,
    // the other is then equatorial. For the moment, this has been solved
    // by returning "equatorial" if it is not axial.
    GeometryType axial_or_equatorial( const Vector3D & point ) const;
    
private:
    std::vector< Vector3D > points_;
    bool is_planar_;
    bool is_chair_;
    bool is_boat_;
    bool is_twisted_boat_;

};

// How parallel / co-planar are they, what is the distance between their centres of mass
class TwoRingsRingAnalyser
{
public:
    
    explicit TwoRingsRingAnalyser( const std::vector< Vector3D > & points );

    bool is_planar() const { return is_planar_; }
    bool is_chair() const { return is_chair_; }
    bool is_boat() const { return is_boat_; }

private:
    std::vector< Vector3D > points_;
    bool is_planar_;
    bool is_chair_;
    bool is_boat_;

};

#endif // ANALYSERINGS_H

