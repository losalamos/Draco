//----------------------------------*-C++-*----------------------------------//
// AMR_Builder.hh
// Thomas M. Evans
// Sat Jul 25 14:25:16 1998
//---------------------------------------------------------------------------//
// @> AMR_Builder class Header file
//---------------------------------------------------------------------------//

#ifndef __imc_AMR_Builder_hh__
#define __imc_AMR_Builder_hh__

//===========================================================================//
// class AMR_Builder - 
//
// Purpose : Build an AMR mesh
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

#include "Names.hh"
#include "AMR_Interface.hh"
#include "Coord_sys.hh"
#include "Layout.hh"
#include "OS_Mesh.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include <vector>
#include <string>

IMCSPACE

// stl components
using std::vector;
using std::string;

// draco components
using dsxx::SP;

class AMR_Builder
{
private:
  // data needed to make an AMR Mesh

  // coordinate system
    string coord_sys;
  // layout data
    vector<int> laydata;
  // vertex data
    vector<double> vertices;

  // member functions for building AMR_Meshes
    SP<Coord_sys> build_Coord();
    SP<Layout>    build_Layout(SP<Coord_sys>);
    SP<OS_Mesh>   build_Mesh(SP<Coord_sys>, SP<Layout>);

public:
  // explicit constructor
    explicit AMR_Builder(SP<AMR_Interface>);

  // build an AMR mesh
    SP<OS_Mesh> build_Mesh();
};

CSPACE

//---------------------------------------------------------------------------//
// inline functions for AMR_Interface
//---------------------------------------------------------------------------//

#endif                          // __imc_AMR_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of imc/AMR_Builder.hh
//---------------------------------------------------------------------------//
