#ifndef MATHSFUNCTIONS_H
#define MATHSFUNCTIONS_H

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

#include "BasicMathsFunctions.h"

#include <cmath>
#include <cstddef> // For definition of size_t.
#include <vector>

typedef double (*Function)( const double x );

// Returns the x-values and weights necessary for Gauss-Legendre quadrature.
// x1 is the lower limit for integration, x2 the upper limit. npoints is the number of points.
// x contains the x values, w contains the weights.
void Gauss_Legendre_quadrature( const double x1, const double x2, const size_t npoints, std::vector< double > & x, std::vector< double > & w );

double bisection( const Function f, const double target_y_value, const double initial_x_value, const double tolerance = TOLERANCE );

// Simple linear interpolation between sampling points.
double integral( const Function f, const double start, const double end, const size_t npoints );

// Assumes that the std::vector< double > represents f(x) values at a regular interval.
// Simple linear interpolation between sampling points.
double integral( const std::vector< double > & y_i, const double interval );

// Evaluates the integral by drawing random values for x.
// This has a minor advantage that you if you want to double the number of sampling points
// to improve the accuracy, you can simply call the function twice with the same number
// of points and average the results. If you do that with integral(), you would
// recalculate everything you have done before or you would need
// sophisticated bookkeeping of what you have and have not yet done.
// Some quick numerical tests with sin(x), exp(x) and f(x) = constant show that this is a really poor way
// to calculate an integral and the function integral() gives vastly superior results, especially
// if very few points (say, 10) are used.
double integral_Monte_Carlo( const Function f, const double start, const double end, const size_t npoints );

double integral_Gauss_Legendre_quadrature( const Function f, const double start, const double end, const size_t npoints );

double Legendre_polynomial( const size_t order, const double x );
void Legendre_polynomial_and_derivative( const size_t order, const double x, double & y, double & dydx );

// 0 <= m <= l, -1 <= x <= 1.
double associated_Legendre_polynomial( const size_t l, const size_t m, const double x );

// Must return size_t, because we use recursion.
// Simplistic algorithm, will overflow very quickly.
size_t factorial( const size_t n );

// double factorial is n!!.
// Must return size_t, because we use recursion.
// Simplistic algorithm, will overflow very quickly.
size_t double_factorial( const size_t n );

// Adds a list of doubles trying to avoid adding very small to very large numbers.
double add_doubles( const std::vector< double > & values );

// Adds a list of doubles trying to avoid adding very small to very large numbers.
double add_absolute_doubles( const std::vector< double > & values );

// Adds a list of doubles trying to avoid adding very small to very large numbers.
double add_squared_doubles( const std::vector< double > & values );

double calculate_average( const std::vector< double > & values );
double calculate_minimum( const std::vector< double > & values );
double calculate_maximum( const std::vector< double > & values );
size_t calculate_minimum( const std::vector< size_t > & values );
size_t calculate_maximum( const std::vector< size_t > & values );

// Calculates sqrt( x^2 + y^2 ) trying to prevent underflow or overflow.
double hypothenuse( const double x, const double y );

// Centered around 0.0, area normalised to 1.0.
// For pseudo-Voigt.
double Lorentzian( const double x, const double FWHM );

// Centered around 0.0, area normalised to 1.0.
// For pseudo-Voigt.
double Gaussian( const double x, const double FWHM );

// Centered around 0.0, area normalised to 1.0.
// Needs: FWHM (in degrees 2theta), eta, 2theta w.r.t. 0.0
double pseudo_Voigt( const double x, const double FWHM, const double eta = 0.9 );

// Returns a random number between 0.0 and 1.0, inclusive.
double uniform_distribution_1();

// Returns a random number between 0.0 and 1.0, exclusive.
double uniform_distribution_2();

double normal_distribution( const double mean = 0.0, const double sigma = 1.0 );

size_t Poisson_distribution( const double mean );

// For testing.
size_t Poisson_distribution_exact( const double mean );
size_t Poisson_distribution_Gaussian_0( const double mean );
size_t Poisson_distribution_Gaussian_1( const double mean );
size_t Poisson_distribution_Gaussian_2( const double mean );

double calculate_sample_RMSD( const std::vector< double > & lhs, const std::vector< double > & rhs );

double ln( const double x );

#endif // MATHSFUNCTIONS_H

