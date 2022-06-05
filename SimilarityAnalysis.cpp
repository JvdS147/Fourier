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

#include "SimilarityAnalysis.h"

#include "CorrelationMatrix.h"
#include "CrystalStructure.h"
#include "FileList.h"
#include "PowderPattern.h"
#include "PowderPatternCalculator.h"
#include "ReadCif.h"
#include "Utilities.h"

#include <iostream>
#include <vector>

// ********************************************************************************

CorrelationMatrix calculate_correlation_matrix( const FileList & file_list )
{
    std::vector< PowderPattern > powder_patterns;
    powder_patterns.reserve( file_list.size() );
    Angle two_theta_start( 3.0, Angle::DEGREES );
    Angle two_theta_end(  35.0, Angle::DEGREES );
    Angle two_theta_step( 0.01, Angle::DEGREES );
    double FWHM( 0.1 );
    for ( size_t i( 0 ); i != file_list.size(); ++i )
    {
        CrystalStructure crystal_structure;
        std::cout << "Now reading cif... " + file_list.value( i ).full_name() << std::endl;
        read_cif( file_list.value( i ), crystal_structure );
        crystal_structure.apply_space_group_symmetry();
        std::cout << "Now calculating powder pattern... " + size_t2string( i, 4, '0' ) << std::endl;
        PowderPatternCalculator powder_pattern_calculator( crystal_structure );
        powder_pattern_calculator.set_two_theta_start( two_theta_start );
        powder_pattern_calculator.set_two_theta_end( two_theta_end );
        powder_pattern_calculator.set_two_theta_step( two_theta_step );
        powder_pattern_calculator.set_FWHM( FWHM );
        PowderPattern powder_pattern;
        powder_pattern_calculator.calculate( powder_pattern );
        ReflectionList reflection_list = powder_pattern_calculator.reflection_list();
        powder_patterns.push_back( powder_pattern );
    }
    CorrelationMatrix result( powder_patterns.size() );
    // To speed things up, for each powder pattern pre-calculate the weighted cross-correlation function
    std::vector< double > sqrt_weighted_cross_correlations;
    sqrt_weighted_cross_correlations.reserve( powder_patterns.size() );
    // When experimental patterns are involved, the default value is 3.0.
    Angle l = Angle( 1.0, Angle::DEGREES );
    std::cout << "Now starting the precalculations" << std::endl;
    for ( size_t i( 0 ); i != powder_patterns.size(); ++i )
        sqrt_weighted_cross_correlations.push_back( sqrt( weighted_cross_correlation( powder_patterns[i], powder_patterns[i], l ) ) );
    std::cout << "Precalculations done" << std::endl;
    size_t iTotal( 0 );
    for ( size_t i( 0 ); i != powder_patterns.size(); ++i )
    {
        for ( size_t j( i+1 ); j != powder_patterns.size(); ++j )
        {
            double value = weighted_cross_correlation( powder_patterns[i], powder_patterns[j], l ) / ( sqrt_weighted_cross_correlations[i] * sqrt_weighted_cross_correlations[j] );
            result.set_value( i, j, value );
            ++iTotal;
            if ( (iTotal % 100) == 0 )
                std::cout << size_t2string( iTotal ) + " comparisons done" << std::endl;
        }
    }
    return result;
}
    
// ********************************************************************************
