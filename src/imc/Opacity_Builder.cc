//----------------------------------*-C++-*----------------------------------//
// Opacity_Builder.cc
// Thomas M. Evans
// Fri Mar  6 17:21:36 1998
//---------------------------------------------------------------------------//
// @> Opacity_Builder class implementation file
//---------------------------------------------------------------------------//

#include "imctest/Opacity_Builder.hh"
#include <cassert>

IMCSPACE

//---------------------------------------------------------------------------//
// constructors
//---------------------------------------------------------------------------//
template<class MT>
template<class PT>
Opacity_Builder<MT>::Opacity_Builder(SP<PT> parser, SP<MT> mesh)
    : rho(mesh), temp(mesh), opacity(mesh)
{
  // assign data members from the parser
    zone        = parser->Zone();
    mat_zone    = parser->Mat_zone();
    density     = parser->Density();
    kappa       = parser->Kappa();
    temperature = parser->Temperature();

  // check some asserts verifying size of mesh
    assert (mesh->Num_cells() == rho.Mesh().Num_cells());
    assert (mesh->Num_cells() == temp.Mesh().Num_cells());
    assert (zone.size() == mesh->Num_cells());
}

//---------------------------------------------------------------------------//
// build Mat_State object
//---------------------------------------------------------------------------//
template<class MT>
SP< Mat_State<MT> > Opacity_Builder<MT>::Build_Mat()
{
  // return Mat_State object
    SP< Mat_State<MT> > return_state;

  // number of cells
    int num_cells = rho.Mesh().Num_cells();

  // assign density and temperature to each cell
    for (int cell = 1; cell <= num_cells; cell++)
    {
	int den = density[mat_zone[zone[cell-1]-1]-1];
	int T   = temperature[mat_zone[zone[cell-1]-1]-1];
	rho(cell)  = den;
	temp(cell) = T;
    }
    
  // create Mat_State object
    return_state = new Mat_State<MT>(rho, temp);

  // return Mat_State SP
    return return_state;
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of Opacity_Builder.cc
//---------------------------------------------------------------------------//
