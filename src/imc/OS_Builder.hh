//----------------------------------*-C++-*----------------------------------//
// OS_Builder.hh
// Thomas M. Evans
// Mon Feb  9 16:16:07 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __imctest_OS_Builder_hh__
#define __imctest_OS_Builder_hh__

//===========================================================================//
// class OS_Builder - 
//
// Date created : 2-9-98
// Purpose      : builds an OS_Mesh object, parses the input file, 
//                redistributes the Mesh for parallel purposes
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

#include "Names.hh"
#include <vector>
#include <string>
#include "SP.hh"

IMCSPACE

class OS_Builder
{
private:
  // data
    string input_file;
    string coord_system;
    vector<int> zones;
    vector<int> zone_cells;
    vector<int> cells;
    vector<double> zone_faces;
    vector<string> bnd_cond

  // Parser member functions
    void parser();
  // build Layout helper functions
    SP<Layout> buildLayout(const Coord_sys &coord)
    void Assign2D(Layout &);
    void Assign3D(Layout &);
  // build Coord_sys helper functions
    SP<Coord_sys> buildCoord();
    SP<Coord_sys> buildXY();
    SP<Coord_sys> buildXYZ();
  // build Mesh helper functions
    SP<OS_Mesh> build2dMesh(SP<Coord_sys>, const Layout &);
    SP<OS_Mesh> build3DMesh(SP<Coord_sys>, const Layout &);
public:
    SP<OS_Mesh> buildMesh();

CSPACE

#endif                          // __imctest_OS_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/OS_Builder.hh
//---------------------------------------------------------------------------//
