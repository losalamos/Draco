//----------------------------------*-C++-*----------------------------------//
// OS_Builder.hh
// Thomas M. Evans
// Mon Feb  9 16:16:07 1998
//---------------------------------------------------------------------------//
// @> OS_Builder class header file
//---------------------------------------------------------------------------//

#ifndef __imctest_OS_Builder_hh__
#define __imctest_OS_Builder_hh__

//===========================================================================//
// class OS_Builder - 
//
// Purpose      : builds an OS_Mesh object, parses the input file, 
//                redistributes the Mesh for parallel purposes
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

#include "Names.hh"
#include "Coord_sys.hh"
#include "Layout.hh"
#include "OS_Mesh.hh"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "SP.hh"

IMCSPACE

using std::vector;
using std::string;
using std::ifstream;

class OS_Builder
{
private:
  // data
    string input_file;
    string coord_system;
  // number of fine cells per coarse cell
    vector< vector<int> > fine_cells;
  // coarse edges
    vector< vector<double> > coarse_edge;
  // fine cell edges
    vector< vector<double> > fine_edge;
  // boundary conditions
    vector<string> bnd_cond;

  // Parser member functions
    void parser2D(ifstream &);
    void parser3D(ifstream &);
  // build Layout helper functions
    SP<Layout> buildLayout(const Coord_sys &coord);
    void Assign2D(Layout &);
    void Assign3D(Layout &);
  // build Coord_sys helper functions
    SP<Coord_sys> buildCoord();
    SP<Coord_sys> buildXY();
    SP<Coord_sys> buildXYZ();
  // build Mesh helper functions
    SP<OS_Mesh> build2DMesh(SP<Coord_sys>, const Layout &);
    SP<OS_Mesh> build3DMesh(SP<Coord_sys>, const Layout &);
public:
    explicit OS_Builder(const string &infile) 
	: input_file(infile), coord_system(""), fine_cells(0), 
	  coarse_edge(0), fine_edge(0), bnd_cond(0)
    {}
    SP<OS_Mesh> buildMesh();    
    void parser();
};

CSPACE

#endif                          // __imctest_OS_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/OS_Builder.hh
//---------------------------------------------------------------------------//
