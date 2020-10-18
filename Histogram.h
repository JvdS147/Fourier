#ifndef HISTOGRAM_H
#define HISTOGRAM_H

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

#include <cstddef> // For definition of size_t
#include <vector>
/*
  A histogram.
  
  The histogram does not store the data it is given: only the histogram is stored.
  It is the responsibility of the caller to store the original data if it is anticipated
  that the data may need to be binned differently.
  
  There is quite a difference between a histogram for integers and a histogram for doubles.
  Maybe we should either make this class a template class or split this class into two classes.
  Hmmm... a template class would not solve the problem that the implementation must be *different*
  
  Example:

double start( 0.0 );
double finish( 100.0 );
size_t number_of_bins( 100 );
std::vector<double> data;
data.push_back(  );
Histogram histogram( start, finish, number_of_bins );
histogram.add_data( data );
std::cout << histogram.bin( 50 ) << std::endl;
std::cout << histogram.lower_than_start() << std::endl;

*/

// Limits are inclusive
class Histogram
{
public:

    Histogram( const double start, const double finish, const size_t number_of_bins );

    void add_data( const std::vector<double> data );
    
    void add_data( const double data );

    // The index is zero-based
    size_t bin( const size_t i ) const;
    
    size_t size() const { return number_of_bins_; }
    size_t number_of_bins() const { return number_of_bins_; }

// Sum of all bins + lower_than_start + greater_than_finish
//    size_t number_of_data_points() const;

    size_t lower_than_start() const { return lower_than_start_; }
    size_t greater_than_finish() const { return greater_than_finish_; }
    
private:
    double start_;
    double finish_;
    size_t number_of_bins_;
    std::vector<size_t> data_;
    size_t lower_than_start_;
    size_t greater_than_finish_;
};

#endif // HISTOGRAM_H

