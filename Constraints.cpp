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

#include "Constraints.h"
#include "3DCalculations.h"
#include "BasicMathsFunctions.h"
#include "CrystalStructure.h"
#include "Fraction.h"
#include "PointGroup.h"
#include "Utilities.h"

#include <iostream>

namespace
{

bool only_one_column_is_non_zero( const Matrix3D & matrix, const size_t row, size_t & index )
{
    size_t number_of_non_zero_columns( 0 );
    for ( size_t i( 0 ); i != 3; ++i )
    {
        if ( ! nearly_zero( matrix.value( row, i ) ) )
        {
            ++number_of_non_zero_columns;
            index = i;
        }
    }
    return ( number_of_non_zero_columns == 1 );
}

std::string constraint_coordinate( const double x )
{
    std::string result = "= ";
    // I have no idea if 1/24 is the correct value
    Fraction fraction = double2fraction( x, Fraction( 1, 24 ) );
    // This should really be close to exact, because x has already been moved to be exactly on the special position
    if ( ! nearly_equal( fraction.to_double(), x ) )
        std::cout << "constraint_coordinate():Warning: conversion to fraction introduces large error." << std::endl;
    if ( ! fraction.is_pure_fraction() )
        result += int2string( fraction.integer_part() ) + ".0";
    if ( ( ! fraction.is_pure_fraction() ) &&
         ( ! fraction.is_integer() ) )
        result += " + ( ";
    if ( ! fraction.is_integer() )
        result += int2string( fraction.numerator() ) + ".0 / " + int2string( fraction.denominator() ) + ".0";
    if ( ( ! fraction.is_pure_fraction() ) &&
         ( ! fraction.is_integer() ) )
        result += " )";
    result += ";";
    return result;
}

} // namespace

// ********************************************************************************

std::string write_constraints( const CrystalStructure & crystal_structure, const Vector3D & point )
{
    std::string result;
    Vector3D exact_point( point );
    PointGroup point_group = crystal_structure.point_is_on_special_position( exact_point );
    std::string x = double2string( exact_point.x() );
    std::string y = double2string( exact_point.y() );
    std::string z = double2string( exact_point.z() );
    // First special case: general position.
    if ( point_group.nsymmetry_operators() == 1 )
    {
        result = "site C1 x ref_flag " + x + " y ref_flag " + y + " z ref_flag " + z + " occ C 1 beq = bnonh;";
        return result;
    }
    // Second special case: inversion.
    if ( point_group.has_inversion() )
    {
        result = "site C1 x " + constraint_coordinate( exact_point.x() ) + " y " + constraint_coordinate( exact_point.y() ) + " z " + constraint_coordinate( exact_point.z() ) + " occ C 1 beq = bnonh;";
        return result;
    } 
    // Third special case: three-fold rotation axis along the diagonal.
    // @@ There are four body diagonals, below we only cover one.
    // @@ Also, we assume that the body diaganal passes through the origin.
    std::vector< Matrix3D > symmetry_operators;
    symmetry_operators.push_back( Matrix3D() );
    symmetry_operators.push_back( Matrix3D( 0.0, 0.0, 1.0,
                                            1.0, 0.0, 0.0,
                                            0.0, 1.0, 0.0 ) );
    symmetry_operators.push_back( Matrix3D( 0.0, 1.0, 0.0,
                                            0.0, 0.0, 1.0,
                                            1.0, 0.0, 0.0 ) );
    PointGroup three_fold_on_diagonal( symmetry_operators );
    if ( same_symmetry_operators( point_group, three_fold_on_diagonal ) )
    {
        "prm coord_dash " + x;
        result = "site C1 x = coord_dash; y = coord_dash; z = coord_dash; occ C 1 beq = bnonh;";
        return result;
    }
    Matrix3D matrix = point_group.special_position_operator() - Matrix3D();
    // The matrix elements are not necessarily integers, so no need to try to keep them integers.
    Matrix3D T;
    matrix.convert_to_row_echelon_form( T );
    exact_point = T * exact_point;
    size_t rank = 3 - matrix.number_of_zero_rows();
    if ( rank == 0 )
    {
        // All rows are zeros, atom is in a general position.
        std::cout << "site rank == 0" << std::endl;
    }
    else if ( rank == 1 )
    {
        // If there is one column with a non-zero value, that variable is fixed.
        size_t index;
        if ( only_one_column_is_non_zero( matrix, 0, index ) )
        {
            if ( index == 0 )
                result = "site C1 x " + constraint_coordinate( exact_point.x() ) + " y ref_flag " + y + " z ref_flag " + z + " occ C 1 beq = bnonh;";
            else if ( index == 1 )
                result = "site C1 x ref_flag " + x + " y " + constraint_coordinate( exact_point.y() ) + " z ref_flag " + z + " occ C 1 beq = bnonh;";
            else
                result = "site C1 x ref_flag " + x + " y ref_flag " + y + " z " + constraint_coordinate( exact_point.z() ) + " occ C 1 beq = bnonh;";
        }
        else
        {
            std::cout << "Not yet implemented 1" << std::endl;
            std::cout << "point = " << point << std::endl;
                std::cout << "exact_point = " << exact_point << std::endl;
            std::cout << point_group << std::endl;
            std::cout << matrix << std::endl;
        }
    }
    else if ( rank == 2 )
    {
        // If there is one column with a non-zero value, that variable is fixed.
        size_t index_0;
        if ( only_one_column_is_non_zero( matrix, 0, index_0 ) )
        {
            size_t index_1;
            if ( only_one_column_is_non_zero( matrix, 1, index_1 ) )
            {
                if ( index_0 == 0 )
                {
                    if ( index_1 == 1 )
                    {
                        // z is not constrained.
                        result = "site C1 x " + constraint_coordinate( exact_point.x() ) + " y " + constraint_coordinate( exact_point.y() ) + " z ref_flag " + z + " occ C 1 beq = bnonh;";
                    }
                    else
                    {
                        // y is not constrained.
                        result = "site C1 x " + constraint_coordinate( exact_point.x() ) + " y ref_flag " + y + " z " + constraint_coordinate( exact_point.z() ) + " occ C 1 beq = bnonh;";
                    }
                }
                else
                {
                    // x is not constrained.
                    result = "site C1 x ref_flag " + x + " y " + constraint_coordinate( exact_point.y() ) + " z " + constraint_coordinate( exact_point.z() ) + " occ C 1 beq = bnonh;";
                }
            }
            else
            {
                std::cout << "Not yet implemented 2" << std::endl;
                std::cout << "point = " << point << std::endl;
                std::cout << "exact_point = " << exact_point << std::endl;
                std::cout << point_group << std::endl;
                std::cout << matrix << std::endl;
            }
        }
        else
        {
            // x, y and z are related. We have to figure out how they are related.
            // This can be a three-fold axis on a diagonal.
            std::cout << "Not yet implemented 3" << std::endl;
            std::cout << "point = " << point << std::endl;
                std::cout << "exact_point = " << exact_point << std::endl;
            std::cout << point_group << std::endl;
            std::cout << matrix << std::endl;

        }
    }
    else if ( rank == 3 )
    {
        // x, y and z all fixed. E.g. a rotoinversion.
        result = "site C1 x " + constraint_coordinate( exact_point.x() ) + " y " + constraint_coordinate( exact_point.y() ) + " z " + constraint_coordinate( exact_point.z() ) + " occ C 1 beq = bnonh;";
    }
    return result;
}

// ********************************************************************************

