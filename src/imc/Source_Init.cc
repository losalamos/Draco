//----------------------------------*-C++-*----------------------------------//
// Source_Init.cc
// Thomas M. Evans
// Fri Mar 20 13:13:54 1998
//---------------------------------------------------------------------------//
// @> Source_Init class implementation file
//---------------------------------------------------------------------------//

#include "imc/Source_Init.hh"
#include "imc/Global.hh"
#include "imc/Particle.hh"
#include "imc/Particle_Buffer.hh"
#include "rng/Sprng.hh"
#include "ds++/Assert.hh"
#include <cmath>
#include <iomanip>
#include <fstream>

IMCSPACE

// draco components
using RNG::Sprng;
using Global::min;

// STL components
using std::pow;
using std::ofstream;
using std::ios;
using std::setw;
using std::setiosflags;
using std::endl;

//---------------------------------------------------------------------------//
// constructor
//---------------------------------------------------------------------------//
// initialization of member data in general constructor

template<class MT>
template<class IT>
Source_Init<MT>::Source_Init(SP<IT> interface, SP<MT> mesh)
    : evol(mesh), evol_net(mesh), evoltot(0), ess(mesh), fss(mesh), 
      esstot(0), erad(mesh), eradtot(0), ncen(mesh), ncentot(0), nvol(mesh),
      nss(mesh), nvoltot(0), nsstot(0), eloss_vol(0), eloss_ss(0), 
      ew_vol(mesh), ew_ss(mesh), t4_slope(mesh)
{
    Require (interface);
    Require (mesh);

  // get values from interface
    evol_ext = interface->get_evol_ext();
    ss_pos   = interface->get_ss_pos();
    ss_temp  = interface->get_ss_temp();
    rad_temp = interface->get_rad_temp();
    delta_t  = interface->get_delta_t();
    npmax    = interface->get_npmax();
    npnom    = interface->get_npnom();
    dnpdt    = interface->get_dnpdt();
    capacity = interface->get_capacity();
    ss_dist  = interface->get_ss_dist();
    
  // do some assertions to check that all is well
    int num_cells = mesh->num_cells();
    Check (evol_ext.size() == num_cells);
    Check (rad_temp.size() == num_cells);
    Check (ss_pos.size()   == ss_temp.size());

  // temporary assertions
    Check (evol.get_Mesh()     == *mesh);
    Check (ess.get_Mesh()      == *mesh);
    Check (fss.get_Mesh()      == *mesh);
    Check (erad.get_Mesh()     == *mesh);
    Check (ncen.get_Mesh()     == *mesh);
    Check (nvol.get_Mesh()     == *mesh);
    Check (nss.get_Mesh()      == *mesh);
    Check (ew_vol.get_Mesh()   == *mesh);
    Check (ew_ss.get_Mesh()    == *mesh);
    Check (evol_net.get_Mesh() == *mesh);
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
    npwant = min(npmax, static_cast<int>(npnom + dnpdt * delta_t));
    Check (npwant != 0);

  // on first pass do initial census, on all cycles calc source energies 
    if (cycle == 1)
	calc_initial_census(*mesh, *opacity, *state, *rcontrol);
    else
	calc_source_energies(*opacity, *state);
	
  // calculate source numbers
    calc_source_numbers(*opacity);

  // calculate the slopes of T_electron^4
    calc_t4_slope(*mesh, *state);
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
    if (ncentot > 0)
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
void Source_Init<MT>::calc_source_numbers(const Opacity<MT> &opacity)
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
    double eloss_cell;

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

      // calculate energy loss due to sampling and adjust evol
	eloss_cell = evol(cell) - nvol(cell) * ew_vol(cell);
	eloss_vol += eloss_cell;
	evol(cell) = nvol(cell) * ew_vol(cell);
	evol_net(cell) = evol(cell) - (1.0 - opacity.get_fleck(cell)) * 
	    evol_ext[cell-1] * evol.get_Mesh().volume(cell) * delta_t;

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
	    ess(surcells[sc]) = Global::a * Global::c * 0.25 * 
		pow(ss_temp[ss],4) * delta_t;

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

  // we should not have made any Random numbers yet
    Require (RNG::rn_stream == 0);

  // loop over cells
    for (int cell = 1; cell <= mesh.num_cells(); cell++)
	for (int i = 1; i <= ncen(cell); i++)
	{
	  // make a new random number for delivery to Particle
	    Sprng random = rcon.get_rn();
	    
	  // sample particle location
	    vector<double> r = mesh.sample_pos(cell, random);

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

  // update the rn_stream constant
    RNG::rn_stream = rcon.get_num();

  // a final assertion
    Ensure (RNG::rn_stream == ncentot);
}

//---------------------------------------------------------------------------//
// calculate slope of T_electron^4 using temporarily calc'd  edge t^4's.

template<class MT>
void Source_Init<MT>::calc_t4_slope(const MT &mesh, 
				    const Mat_State<MT> &state)
{
    double t4_low;
    double t4_high;
    double delta_r;

    for (int cell = 1; cell <= mesh.num_cells(); cell++)
    {
	double t4 = pow(state.get_T(cell), 4);
	for ( int coord = 1; coord <= mesh.get_Coord().get_dim(); coord++)
	{
	    int face_low  = 2*coord - 1;
	    int face_high = 2*coord;
	    int cell_low  = mesh.next_cell(cell, face_low);
	    int cell_high = mesh.next_cell(cell, face_high);

	  // set slope to zero if either side is radiatively reflecting
	    if (cell_low == cell || cell_high == cell)
		t4_slope(coord, cell) = 0.0;

	  // set slope to zero if both sides are radiatively vacuum
	    else if (cell_low == 0 && cell_high == 0)
		t4_slope(coord, cell) = 0.0;

	  // if low side is vacuum, use only two t^4's
	    else if (cell_low == 0)
	    {
		t4_high = pow(state.get_T(cell_high), 4);
		delta_r = 0.5 * (mesh.dim(coord, cell) + 
				 mesh.dim(coord, cell_high));

		t4_slope(coord, cell) = (t4_high - t4) / delta_r;

	      // make sure slope isn't too large so as to give a negative
	      // t4_low.  If so, limit slope so t4_low is zero.
		t4_low = t4 - t4_slope(coord, cell) * 0.5 * 
		    mesh.dim(coord, cell); 
		if (t4_low < 0.0)
		    t4_slope(coord, cell) = 2.0 * t4 / mesh.dim(coord, cell); 
	    }

	  // if high side is vacuum, use only two t^4's
	    else if (cell_high == 0)
	    {
		t4_low = pow(state.get_T(cell_low), 4);
		delta_r = 0.5 * (mesh.dim(coord, cell) + 
				 mesh.dim(coord, cell_low));
		t4_slope(coord, cell) = (t4 - t4_low) / delta_r;

	      // make sure slope isn't too large so as to give a negative
	      // t4_high.  If so, limit slope so t4_high is zero.
		t4_high = t4 + t4_slope(coord,cell) * 0.5 * 
		    mesh.dim(coord, cell); 
		if (t4_high < 0.0)
		    t4_slope(coord, cell) = 2.0 * t4 / mesh.dim(coord, cell);
	    }

	  // no conditions on calculating slope; just do it
	    else
	    {
		t4_low  = pow(state.get_T(cell_low),  4);
		t4_high = pow(state.get_T(cell_high), 4);

		double low_slope = (t4 - t4_low) /
		    (0.5 * (mesh.dim(coord, cell_low) +
			    mesh.dim(coord, cell)) );

		double high_slope = (t4_high - t4) /
		    (0.5 * (mesh.dim(coord, cell) +
			    mesh.dim(coord, cell_high)) );

		double t4_lo_edge = t4 - low_slope  * 0.5 * 
		    mesh.dim(coord, cell);
		double t4_hi_edge = t4 + high_slope * 0.5 *
		    mesh.dim(coord, cell);

		t4_slope(coord, cell) = (t4_hi_edge - t4_lo_edge) / 
		                         mesh.dim(coord, cell);
	    }
	}
    }
}

//---------------------------------------------------------------------------//
// diagnostic functions for Source_Init
//---------------------------------------------------------------------------//
// print out the Source Initialization 

template<class MT>
void Source_Init<MT>::print(ostream &out) const
{
    out << "*** SOURCE INITIALIZATION ***" << endl;
    out << "-----------------------------" << endl;

  // give them the particulars of the source init
    out << setw(35) << setiosflags(ios::right) 
	<< "Number of particles requested: " << setw(10) << npnom << endl;
    out << setw(35) << setiosflags(ios::right)
	<< "Total number calculated: " << setw(10) << npwant << endl;
    out << " ** Breakdown ** " << endl;
    out << setw(20) << "Census Particles: " << setw(10)
	<< ncentot << endl;
    out << setw(20) << "Volume Particles: " << setw(10)
	<< nvoltot << endl;
    out << setw(20) << "Surface Particles: " << setw(10)
	<< nsstot << endl;

    out << endl << " ** Source Energies ** " << endl;
    out.precision(3);
    out << setiosflags(ios::fixed);
    out << setw(10) << setiosflags(ios::right) << "Cell"
        << setw(15) << setiosflags(ios::right) << "Volume ew"
        << setw(15) << setiosflags(ios::right) << "Surface ew" << endl;
    for (int i = 1; i <= ew_vol.get_Mesh().num_cells(); i++)
        out << setw(10) << i << setw(15) << ew_vol(i) << setw(15)
            << ew_ss(i) << endl;	

    out << "-----------------------------" << endl;
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of Source_Init.cc
//---------------------------------------------------------------------------//
