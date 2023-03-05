#ifndef RUNNINGAVERAGEANDESD_H
#define RUNNINGAVERAGEANDESD_H

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

#include "Angle.h"
#include "BasicMathsFunctions.h"
#include "Vector3D.h"

#include <cmath>
#include <stdexcept>
#include <vector>

#include <iostream> // For debugging

/*
  A running average and running estimated standard deviation, i.e. the average and ESD are calculated on the fly as numbers come in. This is <i>not</i> a running average with a fixed-size window.

  At all times do the member functions average() and estimated_standard_deviation() return the average and ESD of all numbers that have been supplied to the class up to that point.

  The numbers that are supplied to the class are <i>not</i> stored. The algorithm comes from Wikipedia and is meant to reduce rounding errors and overflow errors.

  So the purpose of this class is to calculate the average and ESD of a large number of numbers without having to store the individual numbers.

  I have made it a template class so that it works with Angle, I am not sure yet if that was smart.
  Awkward syntax, I had to add all sorts of extra functions for Angle and Vector3D with unexpected behaviour and debugging is much harder because a
  template is only compiled when it is used.
*/
template <class T>
class RunningAverageAndESD
{
public:

    RunningAverageAndESD() : A_(0.0), Q_(0.0), n_(0)
    {
    }

    RunningAverageAndESD( const T value ) : A_(0.0), Q_(0.0), n_(0)
    {
        add_value( value );
    }

    RunningAverageAndESD( const std::vector< T > & values ) : A_(0.0), Q_(0.0), n_(0)
    {
        add_values( values );
    }

    void add_value( const T value )
    {
        ++n_;
        Q_ = Q_ + ( static_cast<double>(n_-1) / n_ ) * square( value - A_ );
        A_ = A_ + ( ( value - A_ ) / n_ );
    }

    void add_values( const std::vector< T > & values )
    {
        for ( typename std::vector< T >::const_iterator it( values.begin() );
              it != values.end();
              ++it )
            add_value( *it );
    }

    // Throws std::runtime_error if no values have been added.
    T average() const
    {
        if ( n_ == 0 )
            throw std::runtime_error( "RunningAverageAndESD::average(): no values added yet." );
        return A_;
    }

    // Returns the estimated standard deviation of the sample (not the population)
    // i.e. the "n-1" ESD. Throws std::runtime_error if 0 or 1 values have been added.
    T estimated_standard_deviation() const
    {
        if ( n_ == 0 )
            throw std::runtime_error( "RunningAverageAndESD::estimated_standard_deviation(): no values added yet." );
        if ( n_ == 1 )
            throw std::runtime_error( "RunningAverageAndESD::estimated_standard_deviation(): only one value added." );
        return std::sqrt( Q_ / (n_-1) );
    }

    void clear() { A_ = 0.0; Q_ = 0.0; n_ = 0; }

    size_t nvalues() const { return n_; }

private:
    T A_;
    T Q_;
    size_t n_;

};

// Specialisation for Angle because square() and std::sqrt() did not exist for Angle.
template<>
class RunningAverageAndESD<Angle>
{
public:

    RunningAverageAndESD() : n_(0)
    {
    }

    RunningAverageAndESD( const Angle value ) : n_(0)
    {
        add_value( value );
    }

    void add_value( const Angle value )
    {
        ++n_;
        Q_ = Q_ + ( static_cast<double>(n_-1) / n_ ) * Angle::from_radians( ( value - A_ ).value_in_radians() * ( value - A_ ).value_in_radians() );
        A_ = A_ + ( ( value - A_ ) / n_ );
    }

    void add_values( const std::vector< Angle > & values )
    {
        for ( std::vector< Angle >::const_iterator it( values.begin() );
              it != values.end();
              ++it )
            add_value( *it );
    }

    // Throws std::runtime_error if no values have been added.
    Angle average() const
    {
        if ( n_ == 0 )
            throw std::runtime_error( "RunningAverageAndESD::average(): no values added yet." );
        return A_;
    }

    // Returns the estimated standard deviation of the sample (not the population)
    // i.e. the "n-1" ESD. Throws std::runtime_error if 0 or 1 values have been added.
    Angle estimated_standard_deviation() const
    {
        if ( n_ == 0 )
            throw std::runtime_error( "RunningAverageAndESD::estimated_standard_deviation(): no values added yet." );
        if ( n_ == 1 )
            throw std::runtime_error( "RunningAverageAndESD::estimated_standard_deviation(): only one value added." );
        return Angle::from_radians( std::sqrt( Q_.value_in_radians() / (n_-1) ) );
    }

    void clear() { A_ = Angle(); Q_ = Angle(); n_ = 0; }

    size_t nvalues() const { return n_; }

private:
    Angle A_;
    Angle Q_;
    size_t n_;

};

// Specialisation for Vector3D because square() and std::sqrt() did not exist for Vector3D.
template<>
class RunningAverageAndESD<Vector3D>
{
public:

    RunningAverageAndESD() : n_(0)
    {
    }

    RunningAverageAndESD( const Vector3D value ) : n_(0)
    {
        add_value( value );
    }

    void add_value( const Vector3D value )
    {
        ++n_;
        Q_ = Q_ + ( static_cast<double>(n_-1) / n_ ) * square( value - A_ );
        A_ = A_ + ( ( value - A_ ) / n_ );
    }

    void add_values( const std::vector< Vector3D > & values )
    {
        for ( std::vector< Vector3D >::const_iterator it( values.begin() );
              it != values.end();
              ++it )
            add_value( *it );
    }

    // Throws std::runtime_error if no values have been added.
    Vector3D average() const
    {
        if ( n_ == 0 )
            throw std::runtime_error( "RunningAverageAndESD::average(): no values added yet." );
        return A_;
    }

    // Returns the estimated standard deviation of the sample (not the population)
    // i.e. the "n-1" ESD. Throws std::runtime_error if 0 or 1 values have been added.
    Vector3D estimated_standard_deviation() const
    {
        if ( n_ == 0 )
            throw std::runtime_error( "RunningAverageAndESD::estimated_standard_deviation(): no values added yet." );
        if ( n_ == 1 )
            throw std::runtime_error( "RunningAverageAndESD::estimated_standard_deviation(): only one value added." );
        return sqrt( Q_ / (n_-1) );
    }

    void clear() { A_ = Vector3D(); Q_ = Vector3D(); n_ = 0; }

    size_t nvalues() const { return n_; }

private:
    Vector3D A_;
    Vector3D Q_;
    size_t n_;

};

#endif // RUNNINGAVERAGEANDESD_H

