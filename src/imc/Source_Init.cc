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
#include "rng/Sprng.hh"
#include <cmath>
#include <iostream>
#include <fstream>

IMCSPACE

using RNG::Sprng;
using Global::min;
using std::pow;
using std::ofstream;

//---------------------------------------------------------------------------//
// public member functions
//---------------------------------------------------------------------------//
// source initialyzer -- this is the main guy

template<class MT>
void Source_Init<MT>::initialize(const MT &mesh, const Opacity<MT> &opacity,
				 const Mat_State<MT> &state, 
				 const Rnd_Control &rcontrol, int cycle)
{
  // calculate number of particles
    calc_num_part();

  // on first pass do initial source census, on all cycles calc source
  // energies 
    if (cycle == 1)
	calc_initial_census(mesh, opacity, state, rcontrol);
    else
	calc_source_energies(opacity, state);
	
  // calculate source numbers
    calc_source_numbers();
}

//---------------------------------------------------------------------------//
// private member functions used in Initialize
//---------------------------------------------------------------------------//

template<class MT>
void Source_Init<MT>::calc_num_part()
{
  // calculate the number of particles used in the timestep
    npwant = min(npmax, npwant + dnpdt * delta_t);
}

//---------------------------------------------------------------------------//

template<class MT>
void Source_Init<MT>::calc_initial_census(const MT &mesh,
					  const Opacity<MT> &opacity,
					  const Mat_State<MT> &state,
					  const Rnd_Control &rcontrol)
{
  // calculate and write the initial census source

  // calc volume emission
    calc_evol(opacity, state);

  // calc surface source emission
    calc_ess();
    
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
  // calc vol emission
    calc_evol(opacity, state);

  // calc surface source
    calc_ess();
}

//---------------------------------------------------------------------------//

template<class MT>
void calc_source_numbers()
{
  // calculate numbers of different quantities necessary for the source

    int nsource       = npwant - ncentot;
    double esource    = evoltot + esstot;
    double part_per_e = nsource / esource;

  // calculate volume source number and surface source numbers

  // re-initialize totals
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
	    pow(state.get_temperature(cell), 4) *
	    evol.get_Mesh().volume(cell) * delta_t + 
	    (1.0 - opacity.get_fleck(cell)) * evol_ext[cell-1] *
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
	    assert (fss(surcells[sc]) == 0);

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
    int ncenguess = eradtot / (evoltot + esstot + eradtot) * npwant;

  // particles per unit energy
    double part_per_e = ncenguess / eradtot;

  // attempt to make all census particles have the same energy weight,
  // iterate
    bool retry = true;
    while (retry)
    {
      // calculate census particles per cell
	ncentot = 0;
	for (int cell = 1; cell <= ncen.get_Mesh().num_cells(); cell++)
	{
	    ncen(cell) = static_cast<int>(erad(cell)*part_per_e+0.5);
	    ncentot += ncen(cell);
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
void Source_Init<MT>::write_initial_census(const MT &mesh, 
					   const Rnd_Control &rcon)
{
  // open census file
    ofstream cen_file("census");

  // loop over cells
    for (int cell = 1; cell <= mesh.num_cells(); cell++)
	for (int i = 1; i <= ncen(cell); i++)
	{
	  // make a new random number for delivery to Particle
	    Sprng random = rcon.get_rn();
	    
	  // sample particle location
	    vector<double> r = sample_pos("uniform", cell, random);

	  // sample particle direction
	    vector<double> omega = sample_dir("isotropic", random);
	    
	  // sample frequency (not now, 1 group)

	  // calculate energy weight
	    double ew = erad(cell)/ncen(cell);

	  // create Particle
	    Particle<MT> particle(r, omega, ew, cell, random);
	    
	  // write particle
	    cen_file << particle;
	}
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of Source_Init.cc
//---------------------------------------------------------------------------//
