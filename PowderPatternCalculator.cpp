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

#include "PowderPatternCalculator.h"
#include "3DCalculations.h"
#include "Angle.h"
#include "CrystallographicCalculations.h"
#include "CrystalStructure.h"
#include "MathsFunctions.h"
#include "PointGroup.h"
#include "PowderPattern.h"
#include "ReflectionList.h"

#include <cmath>
#include <stdexcept>

#include <iostream> // for debugging

// This is of course not flexible enough. This should be a base class "PeakShapeFunction" with
// instantiations like "pseudo_Voigt".

namespace

{
// ********************************************************************************

// Precalculating the peak shape is a speed optimisation.
// Technically, it is an approximation, because the peak centre
// is not exactly at one of the i*2theta_step points, so there is a slight shift of
// on average 1/4 of a 2theta_step.
// This approximation should be fine if 2theta_step is small, such as 0.01 or 0.02.
std::vector< double > peak_shape( const Angle two_theta_step, const double FWHM )
{
    // Find number of points out to the right that need to be calculated to reach 0.1% of intensity at 0.0.
    std::vector<double> values;
    values.reserve( 137 ); // 44 for 1% of intensity, but the value also depends on eta.
    double I100 = pseudo_Voigt( 0.0, FWHM );
    values.push_back( I100 );
    size_t i( 0 );
    double value;
    do
    {
        ++i;
        value = pseudo_Voigt( i * two_theta_step.value_in_degrees(), FWHM );
        values.push_back( value );
    }
    while ( (I100/1000.0) < value );
    // The number of points is now 2*i + 1 (the 1 is 0.0).
    std::vector<double> result;
    result.reserve( 2*i+1 );
    for ( size_t j(0); j < i; ++j )
        result.push_back( values[i-j] );
    result.push_back( values[0] );
    for ( size_t j(0); j < i; ++j )
        result.push_back( values[j+1] );
    return result;
}

} // namespace

// ********************************************************************************

PowderPatternCalculator::PowderPatternCalculator( const CrystalStructure & crystal_structure ):
wavelength_(1.54056),
two_theta_start_(5.0,Angle::DEGREES),
two_theta_end_(30.0,Angle::DEGREES),
two_theta_step_(0.01,Angle::DEGREES),
FWHM_(0.1),
include_zero_point_error_(false),
include_preferred_orientation_(false),
preferred_orientation_direction_( 0, 0, 0 ),
r_(1.0),
include_finger_cox_jephcoat_(false),
finger_cox_jephcoat_( 0.0001, 0.0001 ),
crystal_structure_(crystal_structure)
{
    if ( ! crystal_structure.space_group_symmetry_has_been_applied() )
        std::cout << "PowderPatternCalculator::PowderPatternCalculator( CrystalStructure ): Warning: space-group symmetry has not been applied for input crystal structure."<< std::endl;
//    if ( ! crystal_structure.space_group_symmetry_has_been_applied() )
//        throw std::runtime_error( "PowderPatternCalculator::PowderPatternCalculator( CrystalStructure ): Error: space-group symmetry has not been applied for input crystal structure." );
    Laue_class_ = crystal_structure_.space_group().Laue_class();
}

// ********************************************************************************

void PowderPatternCalculator::set_two_theta_step( const Angle two_theta_step )
{
    two_theta_step_ = two_theta_step;
    if ( two_theta_step_ < Angle::from_degrees( TOLERANCE ) )
         throw std::runtime_error( "PowderPatternCalculator::set_two_theta_step(): Error: value must be positive." );
    // There is an approximation in the calculation of the powder pattern that expects the 2theta step to be small.
    if ( two_theta_step_ > Angle::from_degrees( 0.05 ) )
        std::cout << "PowderPatternCalculator::set_two_theta_step(): Warning: because of an internal approximation, 2theta step is expected to be small." << std::endl;
}

// ********************************************************************************

void PowderPatternCalculator::set_zero_point_error( const Angle zero_point_error )
{
    include_zero_point_error_ = true;
    zero_point_error_ = zero_point_error;
    if ( zero_point_error_.nearly_zero() )
        std::cout << "PowderPatternCalculator::set_zero_point_error(): Warning: zero-point error has been set to 0.0 ." << std::endl;
}

// ********************************************************************************

void PowderPatternCalculator::set_preferred_orientation( const MillerIndices & miller_indices, const double r )
{
    include_preferred_orientation_ = true;
    preferred_orientation_direction_ = miller_indices;
    r_ = r;
    // Check that the PO direction is commensurate with the space-group symmetry.
    MillerIndices reflection( 37, -23, 3 );
    Vector3D PO_vector = reciprocal_lattice_point( preferred_orientation_direction_, crystal_structure_.crystal_lattice() );
    Vector3D H = reciprocal_lattice_point( reflection, crystal_structure_.crystal_lattice() );
    double reference_dot_product = absolute( PO_vector * H );
    for ( size_t i( 0 ); i != Laue_class_.nsymmetry_operators(); ++i )
    {
        MillerIndices equivalent_reflection = reflection * Laue_class_.symmetry_operator( i );
        // Now check that the March-Dollase PO corrections are the same for all of them
        Vector3D H = reciprocal_lattice_point( equivalent_reflection, crystal_structure_.crystal_lattice() );
        double current_dot_product = absolute( PO_vector * H );
        if ( ! nearly_equal( current_dot_product, reference_dot_product ) )
        {
            std::cout << "PowderPatternCalculator::set_preferred_orientation(): Warning: PO direction is not commensurate with space-group symmetry." << std::endl;
            return;
        }
    }
}

// ********************************************************************************

void PowderPatternCalculator::set_finger_cox_jephcoat( const FingerCoxJephcoat & finger_cox_jephcoat )
{
    include_finger_cox_jephcoat_ = true;
    finger_cox_jephcoat_ = finger_cox_jephcoat;
}

// ********************************************************************************

void PowderPatternCalculator::calculate( PowderPattern & powder_pattern )
{
    calculate_reflection_list();
    calculate_structure_factors();
    calculate( reflection_list_, powder_pattern );
}

// ********************************************************************************

void PowderPatternCalculator::calculate_reflection_list( const bool exact )
{
    // Get a list of all reflections.
    // As in Mercury, we ignore two_theta_start_ here.
//    std::cout << "Now generating reflection list... " << std::endl;
    // We add a little extra at the end to avoid cut-off effects.
    Angle theta_end( ( two_theta_end_ / 2.0 ) + Angle::from_degrees( 1.0 ) );
    double one_over_d_min = ( 2.0 * theta_end.sine() ) / wavelength_.wavelength_1();
    double l_max_z = one_over_d_min / crystal_structure_.crystal_lattice().c_star_vector().z();
    double k_max_y = one_over_d_min / crystal_structure_.crystal_lattice().b_star_vector().y();
    double l_max_y = -k_max_y * ( crystal_structure_.crystal_lattice().b_star_vector().z() / crystal_structure_.crystal_lattice().c_star_vector().z() );
    double h_max_x = one_over_d_min / crystal_structure_.crystal_lattice().a_star_vector().x();
    double k_max_x = -h_max_x * ( crystal_structure_.crystal_lattice().a_star_vector().y() / crystal_structure_.crystal_lattice().b_star_vector().y() );
    double l_max_x = -h_max_x * ( crystal_structure_.crystal_lattice().a_star_vector().z() / crystal_structure_.crystal_lattice().c_star_vector().z() )
                     -k_max_x * ( crystal_structure_.crystal_lattice().b_star_vector().z() / crystal_structure_.crystal_lattice().c_star_vector().z() );
    double h_max =           fabs( h_max_x );
    double k_max = std::max( fabs( k_max_x ),           fabs( k_max_y ) );
    double l_max = std::max( fabs( l_max_x ), std::max( fabs( l_max_y ), fabs( l_max_z ) ) );
    // @@ This can be made more efficient because -h, -k, -l is never needed if programmed properly.
    int h_upper = round_to_int( h_max ) + 1;
    int h_lower = -h_upper;
    int k_upper = round_to_int( k_max ) + 1;
    int k_lower = -k_upper;
    int l_upper = round_to_int( l_max ) + 1;
    int l_lower = -l_upper;
    // We can precalculate how many reflections there will be to dimension the vector properly
    // This can be improved (current numbers are wrong)
    size_t cube_size = (2*h_upper+1) * (2*k_upper+1) * (l_upper+1);
    size_t sphere_size = static_cast<size_t>( cube_size/2.0 );
    reflection_list_.reserve( sphere_size );
    for ( int h(h_lower); h <= h_upper; ++h )
    {
        for ( int k(k_lower); k <= k_upper; ++k )
        {
            for ( int l(l_lower); l <= l_upper; ++l )
            {
                MillerIndices current_reflection( h, k, l );
                if ( current_reflection.is_000() )
                    continue;
                if ( is_systematic_absence( current_reflection ) )
                    continue;
                // Get a list of all equivalent reflections
                std::set< MillerIndices > equivalent_reflections = calculate_equivalent_reflections( current_reflection );
                // Only keep the first one
                if ( ! ( current_reflection == *equivalent_reflections.begin() ) )
                    continue;
                Vector3D H = reciprocal_lattice_point( current_reflection, crystal_structure_.crystal_lattice() );
                double d = 1.0 / ( H.length() );
                // Some of the reflections that are generated lead to asin( x ) with x > 1.0, which is an ERROR.
                if ( wavelength_.wavelength_1() > 2.0 * d )
                    continue;
                Angle two_theta = 2.0 * arcsine( wavelength_.wavelength_1() / ( 2.0 * d ) );
                if ( exact )
                {
                    if ( ( two_theta >= two_theta_start_ ) && ( two_theta <= two_theta_end_ ) )
                        reflection_list_.push_back( current_reflection, 1.0, d, equivalent_reflections.size() );
                }
                else
                {
                    if ( two_theta < ( two_theta_end_ + Angle::from_degrees( 0.1 ) ) )
                        reflection_list_.push_back( current_reflection, 1.0, d, equivalent_reflections.size() );
                }
            }
        }
    }
//    reflection_list_.save( "C:\\Data\\for_testing\\PY110.hkl" );
}

// ********************************************************************************

void PowderPatternCalculator::calculate_structure_factors()
{
    if ( ! crystal_structure_.space_group_symmetry_has_been_applied() )
        throw std::runtime_error( "PowderPatternCalculator::calculate_structure_factors(): Error: space-group symmetry has not been applied for input crystal structure." );
//    std::cout << "Now calculating F^2 values... " << std::endl;
    // For each reflection, calculate an intensity.
    for ( size_t i( 0 ); i != reflection_list_.size(); ++i )
    {
        MillerIndices miller_indices( reflection_list_.miller_indices( i ) );
        int h = miller_indices.h();
        int k = miller_indices.k();
        int l = miller_indices.l();
        double d = reflection_list_.d_spacing( i );
        double sine_theta_over_lambda = 1.0 / ( 2.0 * d );
        double cosine_term( 0.0 );
        double sine_term( 0.0 );
        // This can be sped up by sorting the list into lists of atoms of the same element.
        double f0_H  = Element(  1 ).scattering_factor( sine_theta_over_lambda );
        double f0_C  = Element(  6 ).scattering_factor( sine_theta_over_lambda );
        double f0_N  = Element(  7 ).scattering_factor( sine_theta_over_lambda );
        double f0_O  = Element(  8 ).scattering_factor( sine_theta_over_lambda );
        double f0_Cl = Element( 17 ).scattering_factor( sine_theta_over_lambda );
        double f0_I  = Element( 53 ).scattering_factor( sine_theta_over_lambda );
        for ( size_t j( 0 ); j != crystal_structure_.natoms(); ++j )
        {
            double x = crystal_structure_.atom(j).position().x();
            double y = crystal_structure_.atom(j).position().y();
            double z = crystal_structure_.atom(j).position().z();
            Angle argument = Angle::from_radians( 2.0 * CONSTANT_PI * ( h*x + k*y + l*z ) );
            // Get the scattering factor from 2theta, the wavelength and the element of the atom.
            double f0;
            switch ( crystal_structure_.atom(j).element().atomic_number() )
            {
                case  1 : f0 = f0_H; break;
                case  6 : f0 = f0_C; break;
                case  7 : f0 = f0_N; break;
                case  8 : f0 = f0_O; break;
                case 17 : f0 = f0_Cl; break;
                case 53 : f0 = f0_I; break;
                default : f0 = crystal_structure_.atom(j).element().scattering_factor( sine_theta_over_lambda );
            }
            f0 *= crystal_structure_.atom(j).occupancy();
//            if ( ! nearly_equal( crystal_structure_.atom(j).occupancy(), 1.0 ) )
//                std::cout << "occ = " << crystal_structure_.atom(j).occupancy() << std::endl;
            // The Debije-Waller factor
            double T( 1.0 );
            if ( crystal_structure_.atom(j).ADPs_type() == Atom::ANISOTROPIC )
            {
//                std::cout << "Anisotropic ADPs used for atom " << j << std::endl;
                AnisotropicDisplacementParameters ADPs = crystal_structure_.atom(j).anisotropic_displacement_parameters();
                T = exp( -2.0 * square( CONSTANT_PI ) * ( miller_indices * ADPs.U_star( crystal_structure_.crystal_lattice() ) * miller_indices ) );
            }
            else if ( crystal_structure_.atom(j).ADPs_type() == Atom::ISOTROPIC )
            {
                T = exp( -8.0 * square( CONSTANT_PI ) * crystal_structure_.atom(j).Uiso() * square( sine_theta_over_lambda ) );
            }
            else
            {
                // This is what Mercury does according to the manual.
                if ( crystal_structure_.atom(j).element().atomic_number() == 1 )
                    T = exp( -8.0 * square( CONSTANT_PI ) * 0.06 * square( sine_theta_over_lambda ) );
                else
                    T = exp( -8.0 * square( CONSTANT_PI ) * 0.05 * square( sine_theta_over_lambda ) );
            }
            double sine;
            double cosine;
            bool use_sincos = true;
            if ( use_sincos )
                sincos( argument, sine, cosine );
            else
            {
                sine   = argument.sine();
                cosine = argument.cosine();
            }
            sine_term += T * f0 * sine;
            cosine_term += T * f0 * cosine;
        }
        double F_squared = square( cosine_term ) + square( sine_term );
        reflection_list_.set_F_squared( i, F_squared );
    }
//    reflection_list_.save( FileName( "C:\\Data_Win\\ReflectionList_Cpp.hkl" ) );
}

// ********************************************************************************

void PowderPatternCalculator::calculate( const ReflectionList & reflection_list, PowderPattern & powder_pattern )
{
    powder_pattern = PowderPattern( two_theta_start_, two_theta_end_, two_theta_step_ );
    // Calculate one peak with area 1.0.
    std::vector< double > peak_points = peak_shape( two_theta_step_, FWHM_ );
    Vector3D PO_vector;
    if ( include_preferred_orientation_ )
        PO_vector = reciprocal_lattice_point( preferred_orientation_direction_, crystal_structure_.crystal_lattice() );
    // For each reflection, convolute it with a peak shape.
    for ( size_t i( 0 ); i != reflection_list.size(); ++i )
    {
        double multiplicity( 0.0 );
        if ( include_preferred_orientation_ )
        {
      //      std::cout << "PO is included" << std::endl;
            std::set< MillerIndices > equivalent_reflections = calculate_equivalent_reflections( reflection_list.miller_indices( i ) );
            for ( std::set< MillerIndices >::const_iterator it( equivalent_reflections.begin() ); it != equivalent_reflections.end(); ++it )
            {
                Vector3D H = reciprocal_lattice_point( *it, crystal_structure_.crystal_lattice() );
                Angle alpha = angle( PO_vector, H );
                multiplicity += std::pow( square(r_) * square( alpha.cosine() ) + square( alpha.sine() ) / r_, -3.0/2.0 );
            }
        }
        else
            multiplicity = reflection_list.multiplicity( i );
        double peak_intensity = reflection_list.F_squared( i ) * multiplicity;
        Angle theta = arcsine( wavelength_.wavelength_1() / ( 2.0 * reflection_list.d_spacing( i ) ) );
        Angle two_theta = 2.0 * theta;
        if ( include_zero_point_error_ )
        {
            two_theta += zero_point_error_;
            theta = two_theta / 2.0;
        }
        // Multiply by the LP factor.
        double LP_factor = ( 1.0 + square( two_theta.cosine() ) ) / ( 2.0 * two_theta.sine() * theta.sine() );
        peak_intensity *= LP_factor;
        if ( include_finger_cox_jephcoat_ && ( two_theta < Angle::angle_45_degrees() ) )
        {
            // Peak asymmetry just goes on and on, we essentially have to start at 2theta = 0.0.
            Angle current_2theta = two_theta_start_;
            Angle maximum_2theta_value_for_FCJ = two_theta + static_cast<double>(((static_cast<int>(peak_points.size())-1)/2)) * two_theta_step_;
            std::vector< Angle > two_phi_values;
            size_t j( 0 );
            do
            {
                current_2theta = two_theta_start_ + j * two_theta_step_;
                if ( j < powder_pattern.size() )
                    two_phi_values.push_back( current_2theta );
                ++j;
            }
            while ( current_2theta < maximum_2theta_value_for_FCJ );
            std::vector< double > peak_points_2 = finger_cox_jephcoat_.asymmetric_peak( two_theta, two_phi_values, FWHM_ );
            for ( size_t index( 0 ); index != peak_points_2.size(); ++index )
            {
                double intensity = powder_pattern.intensity( index );
                intensity += peak_intensity * peak_points_2[index];
                powder_pattern.set_intensity( index, intensity );
            }
        }
        else
        {
            int peak_centre_index = round_to_int( ( two_theta - two_theta_start_ ) / two_theta_step_ );
            int peak_start_index = peak_centre_index - ( ( static_cast<int>( peak_points.size() ) - 1 ) / 2 );
            for ( size_t j ( 0 ); j != peak_points.size(); ++j )
            {
                // Calculate new index.
                int index = peak_start_index + j;
                if ( ( index < 0 ) || ( index >= powder_pattern.size() ) )
                    continue;
                double intensity = powder_pattern.intensity( index );
                intensity += peak_intensity * pseudo_Voigt( ( powder_pattern.two_theta( index ) - two_theta ).value_in_degrees(), FWHM_ );
                powder_pattern.set_intensity( index, intensity );
            }
        }
    }
    powder_pattern.normalise_highest_peak();
    powder_pattern.recalculate_estimated_standard_deviations();
    powder_pattern.set_wavelength( wavelength_ );
}

// ********************************************************************************

bool PowderPatternCalculator::is_systematic_absence( const MillerIndices H ) const
{
    // Note that we skip the first symmetry operator, which is guaranteed to be the identity.
    for ( size_t i( 1 ); i != crystal_structure_.space_group().nsymmetry_operators(); ++i )
    {
        if ( H * crystal_structure_.space_group().symmetry_operator(i).rotation() == H )
        {
            if ( ! nearly_integer( H * crystal_structure_.space_group().symmetry_operator(i).translation(), 0.05 ) )
                return true;
        }
    }
    return false;
}

// ********************************************************************************

std::set< MillerIndices > PowderPatternCalculator::calculate_equivalent_reflections( const MillerIndices miller_indices ) const
{
    std::set< MillerIndices > result;
    for ( size_t i( 0 ); i != Laue_class_.nsymmetry_operators(); ++i )
        result.insert( miller_indices * Laue_class_.symmetry_operator( i ) );
    return result;
}

// ********************************************************************************

