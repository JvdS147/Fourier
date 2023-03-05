#ifndef WRITECASTEPFILE_H
#define WRITECASTEPFILE_H

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

#include "CrystalStructure.h"

#include <string>

/*
*/
class WriteCASTEPFile
{
public:

    enum JobType { UNIT_CELL_FIXED, UNIT_CELL_FREE, H_ATOMS_ONLY, SS_NMR };

    WriteCASTEPFile();

    explicit WriteCASTEPFile( const CrystalStructure & crystal_structure,
                              const std::string & directory = "",
                              const std::string & base_name = "" );

    JobType job_type() const { return job_type_; }
    void set_job_type( const JobType job_type ) { job_type_ = job_type; }

    // Guaranteed to end in a backslash
    std::string directory() const { return directory_; }
    void set_directory( const std::string & directory );

    std::string base_name() const { return base_name_; }
    void set_base_name( const std::string & base_name ) { base_name_ = base_name; }

    CrystalStructure crystal_structure() const { return crystal_structure_; }
    void set_crystal_structure( const CrystalStructure crystal_structure ) { crystal_structure_ = crystal_structure; crystal_structure_.apply_space_group_symmetry(); }

    void write() const;

private:
    JobType job_type_;
    std::string directory_;
    std::string base_name_;
    CrystalStructure crystal_structure_;
    double cut_off_energy_;
};

#endif // WRITECASTEPFILE_H

