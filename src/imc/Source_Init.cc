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
#include <cmath>

IMCSPACE

using Global::Min;
using std::pow;

//---------------------------------------------------------------------------//
// public member functions
//---------------------------------------------------------------------------//
// source initialyzer -- this is the main guy
template<class MT>
void Source_Init<MT>::initialize(const Opacity<MT> &opacity,
				 const Mat_State<MT> &state)
{
  // calculate number of particles
    calc_num_part();

  // on first pass do initial source census
    if (cycle == 1)
	calc_initial_census(opacity, state);

  // calculate source energies
    calc_source_energies();

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
    npwant = min(npmax, npwant + dupdt * delta_t);
}

//---------------------------------------------------------------------------//

template<class MT>
void Source_Init<MT>::calc_initial_census(const Opacity<MT> &opacity,
					  const Mat_State<MT> &state)
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
    write_initial_census();
}

//---------------------------------------------------------------------------//

template<class MT>
void Source_Init<MT>::calc_source_energies()
{
  // calc vol emission
    calc_evol();

  // calc surface source
    calc_ess();
}

//---------------------------------------------------------------------------//

template<class MT>
void calc_source_numbers()
{
  // calculate numbers of different quantities necessary for the source

    nsource    = npwant - ncensus;
    esource    = evoltot - esstot;
    part_per_e = nsource / esource;

  // calculate volume source number
    calc_nvol();
    
  // calculate surface source number
    calc_nss();
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

template<class MT>
void Source_Init<MT>::calc_ess()
{
  // caculate the total surface source and the surface source in each cell
    
  // reset esstot

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

CSPACE

//---------------------------------------------------------------------------//
//                              end of Source_Init.cc
//---------------------------------------------------------------------------//
