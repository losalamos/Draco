//----------------------------------*-C++-*----------------------------------//
// Opacity_Builder.cc
// Thomas M. Evans
// Fri Mar  6 17:21:36 1998
//---------------------------------------------------------------------------//
// @> Opacity_Builder class implementation file
//---------------------------------------------------------------------------//

#include "imc/Opacity_Builder.hh"
#include "imc/Global.hh"
#include "ds++/Assert.hh"
#include <cmath>

IMCSPACE

// stl components
using std::pow;

//---------------------------------------------------------------------------//
// constructors
//---------------------------------------------------------------------------//

template<class MT>
template<class IT>
Opacity_Builder<MT>::Opacity_Builder(SP<IT> interface)
{
    Require (interface);

  // assign data members from the interface parser
    zone             = interface->get_zone();
    mat_zone         = interface->get_mat_zone();
    density          = interface->get_density();
    kappa            = interface->get_kappa();
    temperature      = interface->get_temperature();
    implicitness     = interface->get_implicit();	
    specific_heat    = interface->get_specific_heat();
    delta_t          = interface->get_delta_t();
    analytic_opacity = interface->get_analytic_opacity();

  // some crucial checks about our data
    Check (implicitness >= 0 && implicitness <= 1);
    Check (delta_t > 0);
}

//---------------------------------------------------------------------------//
// build Mat_State object
//---------------------------------------------------------------------------//

template<class MT>
SP< Mat_State<MT> > Opacity_Builder<MT>::build_Mat(SP<MT> mesh)
{
    Require (mesh);

  // return Mat_State object
    SP< Mat_State<MT> > return_state;

  // number of cells
    int num_cells = mesh->num_cells();

  // make CCSF objects needed for Mat_State
    typename MT::CCSF_double rho(mesh);
    typename MT::CCSF_double temp(mesh);
    typename MT::CCSF_double Cv(mesh);

  // assign density and temperature to each cell
    for (int cell = 1; cell <= num_cells; cell++)
    {
	double den  = density[mat_zone[zone[cell-1]-1]-1];
	double T    = temperature[mat_zone[zone[cell-1]-1]-1];
	double heat = specific_heat[mat_zone[zone[cell-1]-1]-1];
	rho(cell)   = den;
	temp(cell)  = T;
	Cv(cell)    = heat * mesh->volume(cell) * den;;
    }
    
  // create Mat_State object
    return_state = new Mat_State<MT>(rho, temp, Cv);

  // return Mat_State SP
    return return_state;
}

//---------------------------------------------------------------------------//
// build Opacity object
//---------------------------------------------------------------------------//

template<class MT>
SP< Opacity<MT> > Opacity_Builder<MT>::build_Opacity(SP<MT> mesh,
						     SP<Mat_State<MT> > mat)
{
    Require (mesh);
    Require (mat);
    Require (mat->num_cells() == mesh->num_cells());

  // return Opacity object
    SP< Opacity<MT> > return_opacity;

  // number of cells
    int num_cells = mesh->num_cells();

  // instantiate objects needed for Opacity build
    typename MT::CCSF_double sigma_abs(mesh);
    typename MT::CCSF_double fleck(mesh);

  // calculate and assign opacities to each cell
    for (int cell = 1; cell <= num_cells; cell++)
    {
      // get updated temperature from mat_state
	double T     = mat->get_T(cell);

      // calculate real absorption opacities
	double k;
	if (analytic_opacity == "straight")
	    k = kappa[mat_zone[zone[cell-1]-1]-1];
	else if (analytic_opacity == "tcube")
	    k = kappa[mat_zone[zone[cell-1]-1]-1] / (T*T*T);
	else
	    Check (0);
	double den      = mat->get_rho(cell);
	sigma_abs(cell) = den * k;

      // calculate Fleck factor, for 1 group sigma_abs = planck
	double Cv    = mat->get_Cv(cell);
	Insist(Cv > 0, "The specific heat is <= 0!");
	double beta  = 4.0 * Global::a * T*T*T / Cv;
	double denom = implicitness * beta * Global::c * delta_t *
	    sigma_abs(cell);
	fleck(cell)  = 1.0 / (1.0 + denom);
    }

  // create Opacity object
    return_opacity = new Opacity<MT>(sigma_abs, sigma_abs, fleck);

  // return Opacity SP
    return return_opacity;
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of Opacity_Builder.cc
//---------------------------------------------------------------------------//
