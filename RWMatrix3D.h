#ifndef RWMATRIX3D_H
#define RWMATRIX3D_H

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

class FileName;
class Matrix3D;

// Should these have been member functions? But then all class have to know about
// things like FileName and TextFileReader, which I do not really want.
// A powder diffraction pattern is different, because it is so natural to
// have a powder diffraction pattern in the form of a file,
// and the file format is dictated by external circumstances whereas for a matrix it
// should be clear that we are defining it ourselves.

// File extension: .mat
// Format:

//[ [ , , ], [ , , ], [ , , ] ]

// At the moment this is not what is written out by e.g. Matrix3D::show(), obviously, it should.

// Because we cannot overload on return type, at the moment it seems smart to pass the object
// to be read as an argument.
void read( const FileName & file_name, Matrix3D & matrix );

void write( const FileName & file_name, const Matrix3D & matrix );

#endif // RWMATRIX3D_H

