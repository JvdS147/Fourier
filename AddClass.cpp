/* *********************************************
Copyright (c) 2013-2025, Cornelis Jan (Jacco) van de Streek
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

#include "AddClass.h"
#include "FileName.h"
#include "StringFunctions.h"
#include "TextFileWriter.h"

#include <stdexcept>

// ********************************************************************************

void add_class( const std::string & class_name )
{
    std::string directory = "C:\\Users\\jacco\\Documents\\Cpp";
    if ( FileName( directory, class_name, "h" ).exists() )
        throw std::runtime_error( "add_class(): ERROR .h File exists." );
    if ( FileName( directory, class_name, "cpp" ).exists() )
        throw std::runtime_error( "add_class(): ERROR .cpp File exists." );
    if ( FileName( directory, "Test"+ class_name, "cpp" ).exists() )
        throw std::runtime_error( "add_class(): ERROR Test File exists." );
    TextFileWriter h_file_writer( FileName( directory, class_name, "h" ) );
    h_file_writer.write_line( "#ifndef " + to_upper( class_name ) + "_H" );
    h_file_writer.write_line( "#define " + to_upper( class_name ) + "_H" );
    h_file_writer.write_line();
    h_file_writer.write_line( "/* *********************************************" );
    h_file_writer.write_line( "Copyright (c) 2013-2025, Cornelis Jan (Jacco) van de Streek" );
    h_file_writer.write_line( "All rights reserved." );
    h_file_writer.write_line();
    h_file_writer.write_line( "Redistribution and use in source and binary forms, with or without" );
    h_file_writer.write_line( "modification, are permitted provided that the following conditions are met:" );
    h_file_writer.write_line( "    * Redistributions of source code must retain the above copyright" );
    h_file_writer.write_line( "      notice, this list of conditions and the following disclaimer." );
    h_file_writer.write_line( "    * Redistributions in binary form must reproduce the above copyright" );
    h_file_writer.write_line( "      notice, this list of conditions and the following disclaimer in the" );
    h_file_writer.write_line( "      documentation and/or other materials provided with the distribution." );
    h_file_writer.write_line( "    * Neither the name of my employers nor the" );
    h_file_writer.write_line( "      names of its contributors may be used to endorse or promote products" );
    h_file_writer.write_line( "      derived from this software without specific prior written permission." );
    h_file_writer.write_line();
    h_file_writer.write_line( "THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\" AND" );
    h_file_writer.write_line( "ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED" );
    h_file_writer.write_line( "WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE" );
    h_file_writer.write_line( "DISCLAIMED. IN NO EVENT SHALL CORNELIS JAN VAN DE STREEK BE LIABLE FOR ANY" );
    h_file_writer.write_line( "DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES" );
    h_file_writer.write_line( "(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;" );
    h_file_writer.write_line( "LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND" );
    h_file_writer.write_line( "ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT" );
    h_file_writer.write_line( "(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS" );
    h_file_writer.write_line( "SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE." );
    h_file_writer.write_line( "********************************************* */" );
    h_file_writer.write_line();
    h_file_writer.write_line( "/*" );
    h_file_writer.write_line();
    h_file_writer.write_line( "*/" );
    h_file_writer.write_line( "class " + class_name );
    h_file_writer.write_line( "{" );
    h_file_writer.write_line( "public:" );
    h_file_writer.write_line();
    h_file_writer.write_line( "    // Default constructor" );
    h_file_writer.write_line( "    " + class_name + "();" );
    h_file_writer.write_line();
    h_file_writer.write_line( "private:" );
    h_file_writer.write_line();
    h_file_writer.write_line( "};" );
    h_file_writer.write_line();
    h_file_writer.write_line( "#endif // " + to_upper( class_name ) + "_H" );
    h_file_writer.write_line();

    TextFileWriter cpp_file_writer( FileName( directory, class_name, "cpp" ) );
    cpp_file_writer.write_line( "/* *********************************************" );
    cpp_file_writer.write_line( "Copyright (c) 2013-2025, Cornelis Jan (Jacco) van de Streek" );
    cpp_file_writer.write_line( "All rights reserved." );
    cpp_file_writer.write_line();
    cpp_file_writer.write_line( "Redistribution and use in source and binary forms, with or without" );
    cpp_file_writer.write_line( "modification, are permitted provided that the following conditions are met:" );
    cpp_file_writer.write_line( "    * Redistributions of source code must retain the above copyright" );
    cpp_file_writer.write_line( "      notice, this list of conditions and the following disclaimer." );
    cpp_file_writer.write_line( "    * Redistributions in binary form must reproduce the above copyright" );
    cpp_file_writer.write_line( "      notice, this list of conditions and the following disclaimer in the" );
    cpp_file_writer.write_line( "      documentation and/or other materials provided with the distribution." );
    cpp_file_writer.write_line( "    * Neither the name of my employers nor the" );
    cpp_file_writer.write_line( "      names of its contributors may be used to endorse or promote products" );
    cpp_file_writer.write_line( "      derived from this software without specific prior written permission." );
    cpp_file_writer.write_line();
    cpp_file_writer.write_line( "THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\" AND" );
    cpp_file_writer.write_line( "ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED" );
    cpp_file_writer.write_line( "WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE" );
    cpp_file_writer.write_line( "DISCLAIMED. IN NO EVENT SHALL CORNELIS JAN VAN DE STREEK BE LIABLE FOR ANY" );
    cpp_file_writer.write_line( "DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES" );
    cpp_file_writer.write_line( "(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;" );
    cpp_file_writer.write_line( "LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND" );
    cpp_file_writer.write_line( "ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT" );
    cpp_file_writer.write_line( "(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS" );
    cpp_file_writer.write_line( "SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE." );
    cpp_file_writer.write_line( "********************************************* */" );
    cpp_file_writer.write_line();
    cpp_file_writer.write_line( "#include \"" + class_name + ".h\"");
    cpp_file_writer.write_line();
    cpp_file_writer.write_line( "// ********************************************************************************" );
    cpp_file_writer.write_line();
    cpp_file_writer.write_line( "// ********************************************************************************" );
    cpp_file_writer.write_line();

    TextFileWriter test_file_writer( FileName( directory, "Test"+ class_name, "cpp" ) );
    test_file_writer.write_line( "/* *********************************************" );
    test_file_writer.write_line( "Copyright (c) 2013-2025, Cornelis Jan (Jacco) van de Streek" );
    test_file_writer.write_line( "All rights reserved." );
    test_file_writer.write_line();
    test_file_writer.write_line( "Redistribution and use in source and binary forms, with or without" );
    test_file_writer.write_line( "modification, are permitted provided that the following conditions are met:" );
    test_file_writer.write_line( "    * Redistributions of source code must retain the above copyright" );
    test_file_writer.write_line( "      notice, this list of conditions and the following disclaimer." );
    test_file_writer.write_line( "    * Redistributions in binary form must reproduce the above copyright" );
    test_file_writer.write_line( "      notice, this list of conditions and the following disclaimer in the" );
    test_file_writer.write_line( "      documentation and/or other materials provided with the distribution." );
    test_file_writer.write_line( "    * Neither the name of my employers nor the" );
    test_file_writer.write_line( "      names of its contributors may be used to endorse or promote products" );
    test_file_writer.write_line( "      derived from this software without specific prior written permission." );
    test_file_writer.write_line();
    test_file_writer.write_line( "THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\" AND" );
    test_file_writer.write_line( "ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED" );
    test_file_writer.write_line( "WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE" );
    test_file_writer.write_line( "DISCLAIMED. IN NO EVENT SHALL CORNELIS JAN VAN DE STREEK BE LIABLE FOR ANY" );
    test_file_writer.write_line( "DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES" );
    test_file_writer.write_line( "(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;" );
    test_file_writer.write_line( "LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND" );
    test_file_writer.write_line( "ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT" );
    test_file_writer.write_line( "(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS" );
    test_file_writer.write_line( "SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE." );
    test_file_writer.write_line( "********************************************* */" );
    test_file_writer.write_line();
    test_file_writer.write_line( "#include \"" + class_name + ".h\"" );
    test_file_writer.write_line();
    test_file_writer.write_line( "#include \"TestSuite.h\"" );
    test_file_writer.write_line();
    test_file_writer.write_line( "#include <iostream>" );
    test_file_writer.write_line();
    test_file_writer.write_line( "void test_" + class_name + "( TestSuite & test_suite )" );
    test_file_writer.write_line( "{" );
    test_file_writer.write_line( "    std::cout << \"Now running tests for " + class_name + ".\" << std::endl;" );
    test_file_writer.write_line();
    test_file_writer.write_line( "    {" );
    test_file_writer.write_line( "    " + class_name + " dummy;" );
    test_file_writer.write_line( "    test_suite.test_equality( dummy, , \"" + class_name + "()\" );" );
    test_file_writer.write_line( "    }" );
    test_file_writer.write_line();
    test_file_writer.write_line( "}" );
    test_file_writer.write_line();

}

// ********************************************************************************

