//----------------------------------*-C++-*----------------------------------//
// Source_Init.cc
// Thomas M. Evans
// Fri Mar 20 13:13:54 1998
//---------------------------------------------------------------------------//
// @> Source_Init class implementation file
//---------------------------------------------------------------------------//

#include "imctest/Source_Init.hh"
#include "imctest/Constants.hh"
#include "imctest/Math.hh"
#include "imctest/Particle.hh"
#include "imctest/Particle_Buffer.hh"
#include "rng/Sprng.hh"
#include "ds++/Assert.hh"
#include <cmath>
#include <iostream>
#include <fstream>

IMCSPACE

// draco components
using RNG::Sprng;
using Global::min;

// STL components
using std::pow;
using std::ofstream;

//---------------------------------------------------------------------------//
// constructor
//---------------------------------------------------------------------------//
// initialization of member data in general constructor

template<class MT>
template<class IT>
Source_Init<MT>::Source_Init(SP<IT> interface, SP<MT> mesh)
    : evol(mesh), evoltot(0), ess(mesh), fss(mesh), esstot(0), erad(mesh), 
      eradtot(0), ncen(mesh), ncentot(0), nvol(mesh), nss(mesh),
      nvoltot(0), nsstot(0), eloss_vol(0), eloss_ss(0), ew_vol(mesh),
      ew_ss(mesh)
{
  // get values from interface
    evol_ext = interface->get_evol_ext();
    ss_pos   = interface->get_ss_pos();
    ss_temp  = interface->get_ss_temp();
    rad_temp = interface->get_rad_temp();
    delta_t  = interface->get_delta_t();
    npmax    = interface->get_npmax();
    npwant   = interface->get_npnom();
    dnpdt    = interface->get_dnpdt();
    capacity = interface->get_capacity();
    
  // do some assertions to check that all is well
    int num_cells = mesh->num_cells();
    Check (evol_ext.size() == num_cells);
    Check (rad_temp.size() == num_cells);
    Check (ss_pos.size() == ss_temp.size());

  // temporary assertions
    Check (evol.get_Mesh()   == *mesh);
    Check (ess.get_Mesh()    == *mesh);
    Check (fss.get_Mesh()    == *mesh);
    Check (erad.get_Mesh()   == *mesh);
    Check (ncen.get_Mesh()   == *mesh);
    Check (nvol.get_Mesh()   == *mesh);
    Check (nss.get_Mesh()    == *mesh);
    Check (ew_vol.get_Mesh() == *mesh);
    Check (ew_ss.get_Mesh()  == *mesh);
}

//---------------------------------------------------------------------------//
// public member functions
//---------------------------------------------------------------------------//
// source initialyzer -- this is the main guy

template<class MT>
void Source_Init<MT>::initialize(SP<MT> mesh, SP<Opacity<MT> > opacity, 
				 SP<Mat_State<MT> > state, 
				 SP<Rnd_Control> rcontrol, int cycle)
{
  // check to make sure objects exist
    Check (mesh);
    Check (opacity);
    Check (state);
    Check (rcontrol);

  // calculate number of particles this cycle
    npwant = min(npmax, static_cast<int>(npwant + dnpdt * delta_t));
    Check (npwant != 0);

  // on first pass do initial census, on all cycles calc source energies 
    if (cycle == 1)
	calc_initial_census(*mesh, *opacity, *state, *rcontrol);
    else
	calc_source_energies(*opacity, *state);
	
  // calculate source numbers
    calc_source_numbers();
}

//---------------------------------------------------------------------------//
// private member functions used in Initialize
//---------------------------------------------------------------------------//

template<class MT>
void Source_Init<MT>::calc_initial_census(const MT &mesh,
					  const Opacity<MT> &opacity,
					  const Mat_State<MT> &state,
					  Rnd_Control &rcontrol)
{
  // calculate and write the initial census source

  // calc volume emission and surface source energies
    calc_source_energies(opacity, state);
    
  // calc radiation energy
    calc_erad();

  // calc initial number of census particles
    calc_ncen_init();

  // write out the initial census
    write_initial_census(mesh, rcontrol);  
}

//---------------------------------------------------------------------------//

template<class MT>
void Source_Init<MT>::calc_source_energies(const Opacity<MT> &opacity, 
					   const Mat_State<MT> &state)
{
  // calc volume emission energy per cell, total
    calc_evol(opacity, state);

  // calc surface source energy per cell, total
    calc_ess();
}

//---------------------------------------------------------------------------//

template<class MT>
void Source_Init<MT>::calc_source_numbers()
{
  // for ss and volume emission, calculate numbers of particles, energy
  // weight, and energy loss due to inadequate sampling

    int nsource       = npwant - ncentot;
    double esource    = evoltot + esstot;
    double part_per_e;
    if (esource > 0)
	part_per_e = nsource / esource;
    else 
	part_per_e = 0.0;

  // calculate volume source number and surface source numbers

  // initialize totals each cycle
    nvoltot = 0;
    nsstot  = 0;
    eloss_vol = 0.0;
    eloss_ss  = 0.0;

  // loop over cells to calculate number of vol and ss particles
    for (int cell = 1; cell <= nvol.get_Mesh().num_cells(); cell++)
    {
      // calculate volume source info
	nvol(cell) = static_cast<int>(evol(cell) * part_per_e + .5);
	if (nvol(cell) > 0)
	    ew_vol(cell) = evol(cell) / nvol(cell);
	else
	    ew_vol(cell) = 0.0;
	nvoltot += nvol(cell);
	eloss_vol += evol(cell) - nvol(cell) * ew_vol(cell);

      // calculate surface source info
	nss(cell) = static_cast<int>(ess(cell) * part_per_e + .5);
	if (nss(cell) > 0)
	    ew_ss(cell) = ess(cell) / nss(cell);
	else
	    ew_ss(cell) = 0.0;
	nsstot += nss(cell);
	eloss_ss += ess(cell) - nss(cell) * ew_ss(cell);
    }	
}

//---------------------------------------------------------------------------//
// private member functions which calculate source parameters
//---------------------------------------------------------------------------//

template<class MT>
void Source_Init<MT>::calc_evol(const Opacity<MT> &opacity,
				const Mat_State<MT> &state)
{
  // reset evoltot
    evoltot = 0.0;

  // calc volume source and tot volume source
    for (int cell = 1; cell <= evol.get_Mesh().num_cells(); cell++)
    {
      // calc cell centered volume source
	evol(cell) = opacity.fplanck(cell) * Global::a * Global::c *
	    pow(state.get_T(cell), 4) * evol.get_Mesh().volume(cell) * 
	    delta_t + (1.0 - opacity.get_fleck(cell)) * evol_ext[cell-1] *
	    evol.get_Mesh().volume(cell) * delta_t;

      // accumulate evoltot
	evoltot += evol(cell);
    }
}

//---------------------------------------------------------------------------//
// caculate the total surface source and the surface source in each cell
    
template<class MT>
void Source_Init<MT>::calc_ess()
{
  // reset esstot
    esstot = 0.0;

  // loop over surface sources in problem
    for (int ss = 0; ss < ss_pos.size(); ss++)
    {
	vector<int> surcells = ess.get_Mesh().get_surcells(ss_pos[ss]);
	for (int sc = 0; sc < surcells.size(); sc++)
	{      
	  // make sure this cell doesn't already have a surface source
	    Check (fss(surcells[sc]) == 0);

	  // assign energy to surface source cell
	    ess(surcells[sc]) = ss_temp[ss];

	  // assign source face to surface source cell
	    fss(surcells[sc]) = fss.get_Mesh().
		get_bndface(ss_pos[ss], surcells[sc]);

	  // accumulate esstot
	    esstot += ess(surcells[sc]);
	}
    }
}  

//---------------------------------------------------------------------------//
// calculate radiation energy in each cell and total radiation energy

template<class MT>
void Source_Init<MT>::calc_erad()
{
      // reset eradtot
    eradtot = 0.0;

  // calc radiation energy in each cell and accumulate total radiation energy 
    for (int cell = 1; cell <= erad.get_Mesh().num_cells(); cell++)
    {
      // calc cell centered radiation energy
	erad(cell) = Global::a * erad.get_Mesh().volume(cell) *
	    pow(rad_temp[cell-1], 4);

      // accumulate evoltot
	eradtot += erad(cell);
    }
}

//---------------------------------------------------------------------------//
// calculate initial census particles per cell and total

template<class MT>
void Source_Init<MT>::calc_ncen_init()
{
  // first guess at census particles per cell
    Insist ((evoltot+esstot+eradtot) != 0, "You must specify some source!");
    int ncenguess = static_cast<int>(eradtot) / (evoltot + esstot + eradtot) 
	* npwant;

  // particles per unit energy
    double part_per_e;
    if (eradtot > 0)
	part_per_e = ncenguess / eradtot;
    else
	part_per_e = 0.0;

  // attempt to make all census particles have the same energy weight,
  // iterate on number of initial census particles
    bool retry = true;
    while (retry)
    {
      // calculate census particles per cell
	ncentot   =   0;
	eloss_cen = 0.0;
	double ew;
	for (int cell = 1; cell <= ncen.get_Mesh().num_cells(); cell++)
	{
	    ncen(cell) = static_cast<int>(erad(cell) * part_per_e + 0.5);
	    ncentot += ncen(cell);
	    if (ncen(cell) > 0)
		ew = erad(cell) / ncen(cell);
	    else
		ew = 0.0;
	    eloss_cen += erad(cell) - ew * ncen(cell);
	}

      // check to see we haven't exceeded total particles for this cycle
	if (ncentot <= npwant)
	    retry = false;
	else
	    part_per_e = (ncenguess - (ncentot-npwant)) / eradtot;
    }
}

//---------------------------------------------------------------------------//
// write the initial census
	
template<class MT>
void Source_Init<MT>::write_initial_census(const MT &mesh, Rnd_Control &rcon)
{
  // open census file and make Particle Buffer
    ofstream cen_file("census");
    Particle_Buffer<Particle<MT> > buffer(mesh, rcon);

  // loop over cells
    for (int cell = 1; cell <= mesh.num_cells(); cell++)
	for (int i = 1; i <= ncen(cell); i++)
	{
	  // make a new random number for delivery to Particle
	    Sprng random = rcon.get_rn();
	    
	  // sample particle location
	    vector<double> r = mesh.sample_pos("uniform", cell, random);

	  // sample particle direction
	    vector<double> omega = mesh.get_Coord().
		sample_dir("isotropic", random);
	    
	  // sample frequency (not now, 1 group)

	  // calculate energy weight
	    double ew = erad(cell) / ncen(cell);

	  // create Particle
	    Particle<MT> particle(r, omega, ew, cell, random);
	    
	  // write particle
	    buffer.write_census(cen_file, particle);
	}
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of Source_Init.cc
//---------------------------------------------------------------------------//
