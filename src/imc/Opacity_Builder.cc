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
template<class IT>
Opacity_Builder<MT>::Opacity_Builder(SP<IT> interface, SP<MT> mesh)
    : rho(mesh), temp(mesh), opacity(mesh)
{
  // assign data members from the interface parser
    zone        = interface->Zone();
    mat_zone    = interface->Mat_zone();
    density     = interface->Density();
    kappa       = interface->Kappa();
    temperature = interface->Temperature();

  // check some asserts verifying size of mesh
    assert (mesh->Num_cells() == rho.Mesh().Num_cells());
    assert (mesh->Num_cells() == temp.Mesh().Num_cells());
    assert (mesh->Num_cells() == opacity.Mesh().Num_cells());
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
	double den = density[mat_zone[zone[cell-1]-1]-1];
	double T   = temperature[mat_zone[zone[cell-1]-1]-1];
	rho(cell)  = den;
	temp(cell) = T;
    }
    
  // create Mat_State object
    return_state = new Mat_State<MT>(rho, temp);

  // return Mat_State SP
    return return_state;
}

//---------------------------------------------------------------------------//
// build Opacity object
//---------------------------------------------------------------------------//
template<class MT>
SP< Opacity<MT> > Opacity_Builder<MT>::Build_Opacity()
{
  // return Opacity object
    SP< Opacity<MT> > return_opacity;

  // number of cells
    int num_cells = opacity.Mesh().Num_cells();

  // calculate and assign opacities to each cell
    for (int cell = 1; cell <= num_cells; cell++)
    {
	double k      = kappa[mat_zone[zone[cell-1]-1]-1];
	opacity(cell) = rho(cell) * k;
    }

  // create Opacity object
    return_opacity = new Opacity<MT>(opacity);

  // return Opacity SP
    return return_opacity;
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of Opacity_Builder.cc
//---------------------------------------------------------------------------//
