//----------------------------------*-C++-*----------------------------------//
// testb.cc
// Thomas M. Evans
// Wed Feb 18 10:45:23 1998
//---------------------------------------------------------------------------//
// @> test driver for OS_Builder
//---------------------------------------------------------------------------//

#include "imctest/OS_Parser.hh"
#include "imctest/OS_Builder.hh"
#include "imctest/OS_Mesh.hh"
#include "imctest/Mat_State.hh"
#include "imctest/Opacity_Builder.hh"
#include "imctest/Opacity.hh"
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
    using IMC::Opacity;
    using namespace std;

    SP<OS_Mesh> mesh;
    SP< Mat_State<OS_Mesh> > mat_state;
    SP< Opacity<OS_Mesh> > opacity;

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
	opacity   = opacity_build.Build_Opacity();
    }

    cout << "Coordinate System: " << mesh->Coord().Get_coord() << endl;
    cout << "Mesh Size: " << mesh->Num_cells() << endl;
    for (int cell = 1; cell <= mesh->Num_cells(); cell++)
	mesh->Print(cell);
    cout << endl;

    cout << *mat_state;
    cout << *opacity;

    cout << "Mesh:      " << mesh->Num_cells() << endl;
    cout << "Opacity:   " << opacity->Num_cells() << endl;
    cout << "Mat_State: " << mat_state->Num_cells() << endl;
}

//---------------------------------------------------------------------------//
//                              end of testb.cc
//---------------------------------------------------------------------------//
