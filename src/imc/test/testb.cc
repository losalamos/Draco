//----------------------------------*-C++-*----------------------------------//
// testb.cc
// Thomas M. Evans
// Wed Feb 18 10:45:23 1998
//---------------------------------------------------------------------------//
// @> test driver for OS_Builder
//---------------------------------------------------------------------------//

#include "imctest/test/OS_Parser.hh"
#include "imctest/test/OS_Builder.hh"
#include "imctest/test/OS_Mesh.hh"
#include "imctest/test/Opacity_Builder.hh"
#include "imctest/test/Mat_State.cc"
#include "ds++/SP.hh"
#include <iostream>
#include <string>

main()
{
    using IMC::OS_Parser;
    using IMC::OS_Builder;
    using IMC::OS_Mesh;
    using IMC::Mat_State;
    using IMC::Opacity_Builder;
    using namespace std;

    SP<OS_Mesh> mesh;
    SP< Mat_State<OS_Mesh> > mat_state;

  // scoping blocks
    {
	string infile;
	cout << "Name the input file" << endl;
	cin >> infile;
	
      // run the Parser
	SP<OS_Parser> parser = new OS_Parser(infile);
	parser->Parser();

      // initialize the mesh builder and build mesh
	OS_Builder build(parser);
	mesh = build.Build_Mesh();

      // initialize the Opacity builder and build state 
	Opacity_Builder<OS_Mesh> opacity_build(parser, mesh);
	mat_state = opacity_build.Build_Mat();
    }

    cout << "Coordinate System: " << mesh->Coord().Get_coord() << endl;
    cout << "Mesh Size: " << mesh->Num_cells() << endl;
    for (int cell = 1; cell <= mesh->Num_cells(); cell++)
	mesh->Print(cell);
    cout << endl;

    mat_state->Print(3);

}

//---------------------------------------------------------------------------//
//                              end of testb.cc
//---------------------------------------------------------------------------//
