#ifndef DOUBLEWITHESD_H
#define DOUBLEWITHESD_H

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

#include <string>

/*
  A value with an ESD.
  
  The ESD can be 0.0.

*/
class DoubleWithESD
{
public:

    // Value 0.0, ESD 0.0
    DoubleWithESD(): value_(0.0), estimated_standard_deviation_(0.0), no_esd_(true) {}

    DoubleWithESD( const double value, const double estimated_standard_deviation ): value_(value), estimated_standard_deviation_(estimated_standard_deviation), no_esd_(false) {}

    // At the moment, can only cope with "0.12345(6)"
    explicit DoubleWithESD( std::string input );

    double value() const { return value_;}

    double estimated_standard_deviation() const { return estimated_standard_deviation_; }

    void set_value( const double value ) { value_ = value; }

    void set_estimated_standard_deviation( const double estimated_standard_deviation ) { estimated_standard_deviation_ = estimated_standard_deviation; }

    // 1.234(5), i.e. the value is 1.234 and the ESD is 0.0005. The ESD is rounded to one figure unless it is 1.
    std::string crystallographic_style() const;

private:
    double value_;
    double estimated_standard_deviation_;
    bool   no_esd_;

};

bool are_equal( const DoubleWithESD lhs, const DoubleWithESD rhs );

// Throws if estimated_standard_deviation is zero or negative.
std::string crystallographic_style( const double value, const double estimated_standard_deviation );

std::string crystallographic_style_pad_plus( const double value, const double estimated_standard_deviation );

#endif // DOUBLEWITHESD_H

