#ifndef DOUBLESASTABLE_H
#define DOUBLESASTABLE_H

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

// Formats a list of doubles. Takes care of making all the same length, taking into account negative values, no separator after the last value

#include <string>
#include <vector>

std::vector< std::string > doubles_as_table( const std::vector< double > & doubles,
                                             const size_t ndecimals,
                                             const std::string & separator_begin,
                                             const std::string & separator_middle,
                                             const std::string & separator_end,
                                             const size_t target_length = 100,
                                             const bool force_minus_sign_padding = false,
                                             size_t ncolumns = 0 );

std::vector< std::string > doubles_as_table( const std::vector< double > & doubles,
                                             const size_t ndecimals,
                                             const std::string & separator,
                                             const size_t target_length = 100,
                                             const bool force_minus_sign_padding = false,
                                             size_t ncolumns = 0 );

std::vector< std::string > generate_Gauss_Legendre_quadrature_code( const double x1,
                                                                    const double x2,
                                                                    const size_t npoints,
                                                                    const size_t ndecimals = 6,
                                                                    const size_t target_length = 100,
                                                                    const bool force_minus_sign_padding = false,
                                                                    const size_t ncolumns = 0 );

#endif // DOUBLESASTABLE_H

