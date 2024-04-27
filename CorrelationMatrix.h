#ifndef CORRELATIONMATRIX_H
#define CORRELATIONMATRIX_H

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

class FileName;

#include <cstddef> // For definition of size_t
#include <string>
#include <vector>

/*
  A correlation matrix (or similarity matrix).
  The matrix is symmetric.
  The similarity values must be normalised, so the diagonal contains 1.0.
  Access is boundary-checked, even for the diagonal.
  
  The value that is used for the diagonal can be set (but it must be the same for all entries on the diagonal).
*/
class CorrelationMatrix
{
public:

    explicit CorrelationMatrix( const size_t dimension );
    
    ~CorrelationMatrix();

    size_t size() const { return dimension_; }

    // In keeping with the silly C++ convention: zero-based.
    double value( size_t i, size_t j ) const;

    void set_value( size_t i, size_t j, const double value );

    double value_on_diagonal() const { return value_on_diagonal_; }
    
    void set_value_on_diagonal( const double value ) { value_on_diagonal_ = value; }

    // Diagonal is not included.
    double largest_value() const;

    // Diagonal is not included.
    double smallest_value() const;

    void save( const FileName & file_name ) const;
    
    // @@@not yet implemented>
    void save_lower_triangle( const FileName & file_name ) const;

    // Divides the entries into clusters, all entries more similar than threshold are put into a cluster.
    // If this leads to inconsistencies, e.g. because A = B, B = C, but A != C, then A = C.
    std::vector< std::vector< size_t > > clusters( const double threshold ) const;

    // Prints warning when similarity is higher than grey_area_threshold but lower than threshold.
    std::vector< std::vector< size_t > > clusters( const double grey_area_threshold, const double threshold ) const;

private:
    double* data_ptr_;
    size_t dimension_;
    double value_on_diagonal_;
};

#endif // CORRELATIONMATRIX_H

