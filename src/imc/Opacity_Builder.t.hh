//----------------------------------*-C++-*----------------------------------//
// Opacity_Builder.t.hh
// Thomas M. Evans
// Fri Mar  6 17:21:36 1998
//---------------------------------------------------------------------------//
// @> Opacity_Builder class implementation file
//---------------------------------------------------------------------------//

#include "Opacity_Builder.hh"
#include "Global.hh"
#include <cmath>

namespace rtt_imc 
{

// stl components
using std::pow;

// draco components
using rtt_dsxx::SP;

//---------------------------------------------------------------------------//
// constructors
//---------------------------------------------------------------------------//
// defined inline to avoid explicit instantiation

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

	Check (heat > 0);

        if (analytic_sp_heat == "straight")
	{
	    // specific heat units: [jks/g/keV]
	    dedt(cell)    = heat * mesh->volume(cell) * den;
	    sp_heat(cell) = heat;
	}
        else if (analytic_sp_heat == "tcube")
	{
	    // specific heat multiplier, units: [jks/cm^3/keV^4]
            dedt(cell)    = heat * mesh->volume(cell) * T*T*T;
	    sp_heat(cell) = heat;
	}
	else if (analytic_sp_heat == "dedt")
	{
	    // we are given de/dt, not Cv
	    dedt(cell)    = heat;
	    sp_heat(cell) = heat / (mesh->volume(cell) * den);
	}
    }
    
    // create Mat_State object
    return_state = new Mat_State<MT>(rho, temp, dedt, sp_heat);

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
    // draco components
    using rtt_mc::global::a;
    using rtt_mc::global::c;

    // DBC requirements
    Require (mesh);
    Require (mat);
    Require (mat->num_cells()  == mesh->num_cells());
    Require (mesh->num_cells() == kappa.size());
    Require (mesh->num_cells() == kappa_offset.size());
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
	// index into kappa's
	int cc = cell - 1;

	// get updated temperature from mat_state
	double T = mat->get_T(cell);

	// calculate real absorption opacities
	double den = mat->get_rho(cell);
	if (analytic_opacity == "straight")
	{
	    sigma_abs(cell)     = (kappa_offset[cc] + kappa[cc]) * den;  
	    sigma_thomson(cell) = kappa_thomson[cc] * den;
	}
	else if (analytic_opacity == "tcube")
	{
	    Check (T > 0.0);
	    sigma_abs(cell)     = (kappa_offset[cc] + kappa[cc]/(T*T*T)) * den;
	    sigma_thomson(cell) = kappa_thomson[cc] * den;
	}
	else if (analytic_opacity == "tlinear")
	{
	    Check (T > 0.0);
	    sigma_abs(cell)     = (kappa_offset[cc] + kappa[cc]/T) * den;
	    sigma_thomson(cell) = kappa_thomson[cc] * den;
	}
	else if (analytic_opacity == "opacity")
	{
	    sigma_abs(cell)     = kappa_offset[cc] + kappa[cc];
	    sigma_thomson(cell) = kappa_thomson[cc];
	}
	else
	    Check (0);

	// calculate Fleck factor, for 1 group sigma_abs = planck
	double dedt  = mat->get_dedt(cell);
	Insist(dedt > 0, "The specific heat is <= 0!");
	double beta  = 4.0 * a * T*T*T * mesh->volume(cell) / dedt;
	double denom = implicitness * beta * c * delta_t * sigma_abs(cell);
	fleck(cell)  = 1.0 / (1.0 + denom);
    }

    // create Opacity object
    return_opacity = new Opacity<MT>(sigma_abs, sigma_thomson, sigma_abs, 
				     fleck);

    // return Opacity SP
    return return_opacity;
}

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of Opacity_Builder.t.hh
//---------------------------------------------------------------------------//
