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
    density          = interface->get_density();
    kappa            = interface->get_kappa();
    kappa_thomson    = interface->get_kappa_thomson();
    temperature      = interface->get_temperature();
    specific_heat    = interface->get_specific_heat();
    implicitness     = interface->get_implicit();	
    delta_t          = interface->get_delta_t();
    analytic_opacity = interface->get_analytic_opacity();
    analytic_sp_heat = interface->get_analytic_sp_heat();

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
    Require (mesh->num_cells() == density.size());
    Require (mesh->num_cells() == temperature.size());
    Require (mesh->num_cells() == specific_heat.size());

  // return Mat_State object
    SP< Mat_State<MT> > return_state;

  // number of cells
    int num_cells = mesh->num_cells();

  // make CCSF objects needed for Mat_State
    typename MT::CCSF_double rho(mesh);
    typename MT::CCSF_double temp(mesh);
    typename MT::CCSF_double dedt(mesh);
    typename MT::CCSF_double sp_heat(mesh);

  // assign density and temperature to each cell
    for (int cell = 1; cell <= num_cells; cell++)
    {
	double den    = density[cell-1];
	double T      = temperature[cell-1];
	double heat   = specific_heat[cell-1];
	rho(cell)     = den;
	temp(cell)    = T;
	sp_heat(cell) = heat;
        if (analytic_sp_heat == "straight")
          // specific heat units: [jks/g/keV]
	    dedt(cell) = heat * mesh->volume(cell) * den;
        else if (analytic_sp_heat == "tcube")
          // specific heat multiplier, units: [jks/cm^3/keV^4]
            dedt(cell) = heat * mesh->volume(cell) * T*T*T;
    }
    
  // create Mat_State object
    return_state = new Mat_State<MT>(rho, temp, dedt, sp_heat, 
				     analytic_sp_heat);

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
    Require (mesh->num_cells() == kappa.size());
    Require (mesh->num_cells() == kappa_thomson.size());

  // return Opacity object
    SP< Opacity<MT> > return_opacity;

  // number of cells
    int num_cells = mesh->num_cells();

  // instantiate objects needed for Opacity build
    typename MT::CCSF_double sigma_abs(mesh);
    typename MT::CCSF_double sigma_thomson(mesh);
    typename MT::CCSF_double fleck(mesh);

  // calculate and assign opacities to each cell
    for (int cell = 1; cell <= num_cells; cell++)
    {
      // get updated temperature from mat_state
	double T = mat->get_T(cell);

      // calculate real absorption opacities
	double k;
	if (analytic_opacity == "straight")
	    k = kappa[cell-1];
	else if (analytic_opacity == "tcube")
	    k = kappa[cell-1] / (T*T*T);
	else
	    Check (0);
	double den      = mat->get_rho(cell);
	sigma_abs(cell) = den * k;

      // calculate Thomson scattering opacity
	double kt = kappa_thomson[cell-1];
	if (analytic_opacity == "opacity")
	    sigma_thomson(cell) = kt;
	else
	    sigma_thomson(cell) = den * kt;

      // calculate Fleck factor, for 1 group sigma_abs = planck
	double dedt  = mat->get_dedt(cell);
	Insist(dedt > 0, "The specific heat is <= 0!");
	double beta  = 4.0 * Global::a * T*T*T * mesh->volume(cell) / dedt;
	double denom = implicitness * beta * Global::c * delta_t *
	    sigma_abs(cell);
	fleck(cell)  = 1.0 / (1.0 + denom);
    }

  // create Opacity object
    return_opacity = new Opacity<MT>(sigma_abs, sigma_thomson, sigma_abs, 
				     fleck);

  // return Opacity SP
    return return_opacity;
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of Opacity_Builder.cc
//---------------------------------------------------------------------------//
