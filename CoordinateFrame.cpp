/* *********************************************
Copyright (c) 2013-2024, Cornelis Jan (Jacco) van de Streek
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

#include "CoordinateFrame.h"
#include "3DCalculations.h"

#include <cmath>
#include <stdexcept>

// ********************************************************************************

CoordinateFrame::CoordinateFrame()
{
    *this = CoordinateFrame( NormalisedVector3D( 1.0, 0.0, 0.0 ), NormalisedVector3D( 0.0, 1.0, 0.0 ), NormalisedVector3D( 0.0, 0.0, 1.0 ) );
}

// ********************************************************************************

CoordinateFrame::CoordinateFrame( const NormalisedVector3D & x_axis )
{
    NormalisedVector3D best_attempt( 1.0, 0.0, 0.0 );
    double smallest_absolute_inner_product = std::abs( x_axis * best_attempt );
    if ( std::abs( x_axis * NormalisedVector3D( 0.0, 1.0, 0.0 ) ) < smallest_absolute_inner_product )
    {
        best_attempt = NormalisedVector3D( 0.0, 1.0, 0.0 );
        smallest_absolute_inner_product = x_axis * best_attempt;
    }
    if ( std::abs( x_axis * NormalisedVector3D( 0.0, 0.0, 1.0 ) ) < smallest_absolute_inner_product )
        best_attempt = NormalisedVector3D( 0.0, 0.0, 1.0 );
    NormalisedVector3D y_axis = orthogonalise( x_axis, best_attempt );
    *this = CoordinateFrame( x_axis, y_axis );
}

// ********************************************************************************

CoordinateFrame::CoordinateFrame( const NormalisedVector3D & x_axis, const NormalisedVector3D & y_axis ):
x_axis_(x_axis),
y_axis_(y_axis)
{
    if ( ! nearly_zero( x_axis * y_axis ) )
        throw std::runtime_error( "CoordinateFrame::CoordinateFrame( NormalisedVector3D, NormalisedVector3D ): error: x-axis and y-axis not perpendicular." );
    z_axis_ = NormalisedVector3D( x_axis.y() * y_axis.z() - x_axis.z() * y_axis.y(),
                                  x_axis.z() * y_axis.x() - x_axis.x() * y_axis.z(),
                                  x_axis.x() * y_axis.y() - x_axis.y() * y_axis.x() );
}

// ********************************************************************************

CoordinateFrame::CoordinateFrame( const NormalisedVector3D & x_axis, const NormalisedVector3D & y_axis, const NormalisedVector3D & z_axis ):
x_axis_(x_axis),
y_axis_(y_axis),
z_axis_(z_axis)
{
    // Check that the three axes are orthogonal.
    if ( ! nearly_zero( x_axis * y_axis ) )
        throw std::runtime_error( "CoordinateFrame::CoordinateFrame( NormalisedVector3D, NormalisedVector3D, NormalisedVector3D ): error: x-axis and y-axis not perpendicular." );
    if ( ! nearly_zero( y_axis * z_axis ) )
        throw std::runtime_error( "CoordinateFrame::CoordinateFrame( NormalisedVector3D, NormalisedVector3D, NormalisedVector3D ): error: y-axis and z-axis not perpendicular." );
    if ( ! nearly_zero( z_axis * x_axis ) )
        throw std::runtime_error( "CoordinateFrame::CoordinateFrame( NormalisedVector3D, NormalisedVector3D, NormalisedVector3D ): error: z-axis and x-axis not perpendicular." );
    // Check the coordinate frame is right-handed.
    
}

// ********************************************************************************

//void CoordinateFrame::print() const
//{
//    std::cout << "a = " << a() << ", " <<
//                 "b = " << b() << ", " <<
//                 "c = " << c() << ", " <<
//                 "al = " << alpha() << ", " <<
//                 "be = " << beta()  << ", " <<
//                 "ga = " << gamma() << std::endl;
//}
//
//// ********************************************************************************
//
//void CoordinateFrame::show() const
//{
//    std::cout << "a = " << a_vector_;
//    std::cout << "b = " << b_vector_;
//    std::cout << "c = " << c_vector_;
//    std::cout << "a* = " << a_star_vector_;
//    std::cout << "b* = " << b_star_vector_;
//    std::cout << "c* = " << c_star_vector_;
//}

// ********************************************************************************

//bool nearly_equal( const CoordinateFrame & lhs, const CoordinateFrame & rhs, double tolerance )
//{
//    return 
//}

// ********************************************************************************

