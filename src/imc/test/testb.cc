//----------------------------------*-C++-*----------------------------------//
// testb.cc
// Thomas M. Evans
// Wed Feb 18 10:45:23 1998
//---------------------------------------------------------------------------//
// @> test driver for OS_Builder
//---------------------------------------------------------------------------//

#include "../OS_Parser.hh"
#include "../OS_Builder.hh"
#include "../OS_Mesh.hh"
#include "SP.hh"
#include <iostream>
#include <string>

main()
{
    using IMC::OS_Parser;
    using IMC::OS_Builder;
    using IMC::OS_Mesh;
    using namespace std;

    SP<OS_Mesh> mesh;
    vector<int> zones;
    vector<double> temp;
    vector<int> matzone;

  // scoping blocks
    {
	string infile;
	cout << "Name the input file" << endl;
	cin >> infile;
	
      // run the Parser
	SP<OS_Parser> parser = new OS_Parser(infile);
	parser->Parser();
	zones   = parser->Zone();
	temp    = parser->Temperature();
	matzone = parser->Mat_zone();

      // initialize the builder
	OS_Builder build(parser);
	mesh = build.Build_Mesh();
    }

    cout << "Coordinate System: " << mesh->Coord().Get_coord() << endl;
    cout << "Mesh Size: " << mesh->Num_cells() << endl;
    for (int cell = 1; cell <= mesh->Num_cells(); cell++)
	mesh->Print(cell);
    cout << endl;

    for (int z = 1; z <= 4; z++)
    {
	cout << "Zone : " << z << endl;
	cout << "Mat  : " << matzone[z-1] << endl;
	cout << "Temp : " << temp[matzone[z-1]-1] << endl;
	for (int i = 1; i <= mesh->Num_cells(); i++)
	    if (zones[i-1] == z)
		cout << " " << i;
	cout << endl;
    }

    int cell;
    cout << "Pick a cell" << endl;
    cin >> cell;
    cout << "Cell is in zone " << zones[cell-1] << endl;
    cout << "Cell has temp   " << temp[matzone[zones[cell-1]-1]-1] << endl;
    cout << "Cell has material " << matzone[zones[cell-1]-1] << endl;
}

//---------------------------------------------------------------------------//
//                              end of testb.cc
//---------------------------------------------------------------------------//
