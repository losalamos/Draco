//----------------------------------*-C++-*----------------------------------//
// OS_Parser.hh
// Thomas M. Evans
// Mon Feb 23 17:22:21 1998
//---------------------------------------------------------------------------//
// @> OS_Parser class header file
//---------------------------------------------------------------------------//

#ifndef __imctest_OS_Parser_hh__
#define __imctest_OS_Parser_hh__

//===========================================================================//
// class OS_Parser - 
//
// Purpose : parses a Orthogonal structured mesh input file
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

#include "Names.hh"
#include "OS_Mesh.hh"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

IMCSPACE

using std::string;
using std::vector;
using std::ifstream;

class OS_Parser 
{
private:
  // input file particulars
    string input_file;
    string coord_system;

  // data required for OS_Mesh generation

  // number of fine cells per coarse cell
    vector< vector<int> > fine_cells;
  // recursive total number of fine_cells per coarse cell
    vector< vector<int> > accum_cells;
  // coarse edges
    OS_Mesh::CCVF_a coarse_edge;
  // fine cell edges
    OS_Mesh::CCVF_a fine_edge;
  // boundary conditions
    vector<string> bnd_cond;

  // data required for Opacity<MT> generation

  // zone map
    vector<int> zone;

  // material zones
    vector<int> mat_zone;
    vector<double> density;
    vector<double> kappa;
    vector<double> temperature;

  // Parser member functions

  // OS_Mesh parser functions
    void Parser_mesh(ifstream &);
    void Parser2D(ifstream &);
    void Parser3D(ifstream &);

  // Source member functions
    void Parser_source(ifstream &);

  // Opacity parser functions
    void Parser_opacity(ifstream &);
    void Zone_mapper();
    void Cell_zone(int, int);
    void Cell_zone(int, int, int);
    void Zone_parser(ifstream &);
public:
  // constructor
    explicit OS_Parser(const string &infile)
	: input_file(infile), coord_system(""), fine_cells(0), 
	  accum_cells(0), coarse_edge(0), fine_edge(0), bnd_cond(0), 
	  zone(0), mat_zone(0), density(0), kappa(0), temperature(0)
    {}

  // public Parser member functions
    void Parser();

  // public copy functions for mesh
    string Coordinates() const { return coord_system; }
    const vector<string>& Boundaries() const { return bnd_cond; }
    const OS_Mesh::CCVF_a& Fine_edge() const { return fine_edge; }

  // public copy functions for Opacity<MT>
    const vector<int>& Zone() const { return zone; }
    const vector<int>& Mat_zone() const { return mat_zone; }
    const vector<double>& Density() const { return density; }
    const vector<double>& Kappa() const { return kappa; }
    const vector<double>& Temperature() const { return temperature; }
};

CSPACE

#endif                          // __imctest_OS_Parser_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/OS_Parser.hh
//---------------------------------------------------------------------------//
