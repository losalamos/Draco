//----------------------------------*-C++-*----------------------------------//
// Source_Init.cc
// Thomas M. Evans
// Fri Mar 20 13:13:54 1998
//---------------------------------------------------------------------------//
// @> Source_Init class implementation file
//---------------------------------------------------------------------------//

#include "imc/Source_Init.hh"
#include "imc/Global.hh"
#include "rng/Sprng.hh"
#include "ds++/Assert.hh"
#include <cmath>
#include <iomanip>
#include <fstream>

IMCSPACE

// draco components
using RNG::Sprng;
using Global::min;
using Global::max;

// STL components
using std::pow;
using std::ofstream;
using std::ios;
using std::setw;
using std::setiosflags;
using std::endl;
using std::fabs;

//---------------------------------------------------------------------------//
// constructor
//---------------------------------------------------------------------------//
// initialization of member data in general constructor

template<class MT, class PT>
template<class IT>
Source_Init<MT,PT>::Source_Init(SP<IT> interface, SP<MT> mesh)
    : evol(mesh), evol_net(mesh), evoltot(0),  ess(mesh), fss(mesh), 
      esstot(0), ecen(mesh), ecentot(0), ncen(mesh), ncentot(0), nvol(mesh), 
      nss(mesh), nvoltot(0), nsstot(0), eloss_vol(0), eloss_ss(0), 
      eloss_cen(0), ew_vol(mesh), ew_ss(mesh), ew_cen(mesh), t4_slope(mesh)
{
    Require (interface);
    Require (mesh);

  // get values from interface
    evol_ext   = interface->get_evol_ext();
    rad_source = interface->get_rad_source();
    rad_s_tend = interface->get_rad_s_tend();
    ss_pos     = interface->get_ss_pos();
    ss_temp    = interface->get_ss_temp();
    rad_temp   = interface->get_rad_temp();
    delta_t    = interface->get_delta_t();
    npmax      = interface->get_npmax();
    npnom      = interface->get_npnom();
    dnpdt      = interface->get_dnpdt();
    capacity   = interface->get_capacity();
    ss_dist    = interface->get_ss_dist();
    
  // do some assertions to check that all is well
    int num_cells = mesh->num_cells();
    Check (evol_ext.size()   == num_cells);
    Check (rad_source.size() == num_cells);
    Check (rad_temp.size()   == num_cells);
    Check (ss_pos.size()     == ss_temp.size());

  // temporary assertions
    Check (evol.get_Mesh()     == *mesh);
    Check (ess.get_Mesh()      == *mesh);
    Check (fss.get_Mesh()      == *mesh);
    Check (ecen.get_Mesh()     == *mesh);
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

template<class MT, class PT>
void Source_Init<MT,PT>::initialize(SP<MT> mesh, SP<Opacity<MT> > opacity, 
				    SP<Mat_State<MT> > state, 
				    SP<Rnd_Control> rcontrol, int cycle)
{
  // check to make sure objects exist
    Require (mesh);
    Require (opacity);
    Require (state);
    Require (rcontrol);

  // calculate number of particles this cycle
    npwant = min(npmax, static_cast<int>(npnom + dnpdt*(cycle-1)*delta_t));
    Check (npwant != 0);

  // on first pass do initial census, on all cycles calc source energies 
    if (cycle == 1)
	calc_initial_census(*mesh, *opacity, *state, *rcontrol, cycle);
    else
	calc_source_energies(*opacity, *state, cycle);
	
  // calculate source numbers
    calc_source_numbers(*opacity, cycle);

  // comb the census
    if (census->size() > 0)
	comb_census(*mesh, *rcontrol); 

  // calculate the slopes of T_electron^4
    calc_t4_slope(*mesh, *state);

  // make sure a Census has been created
    Ensure (census);
    Ensure (ncentot == census->size());
}

//---------------------------------------------------------------------------//
// host_init : this is a TEMPORARY fixup to do source initialization on the
// host code

template<class MT, class PT>
void Source_Init<MT,PT>::host_init(SP<MT> mesh, SP<Opacity<MT> > opacity, 
				   SP<Mat_State<MT> > state, 
				   SP<Rnd_Control> rcontrol, int cycle)
{
  // check to make sure objects exist
    Require (mesh);
    Require (opacity);
    Require (state);
    Require (rcontrol);

  // calculate number of particles this cycle
    npwant = min(npmax, static_cast<int>(npnom + dnpdt * delta_t));
    Check (npwant != 0);

  // on first pass do initial census, on all cycles calc source energies 
    if (cycle == 0)
	calc_initial_census(*mesh, *opacity, *state, *rcontrol, cycle);
    else
    {
      // calculate source energies
	calc_source_energies(*opacity, *state, cycle);
	
      // calculate source numbers
	calc_source_numbers(*opacity, cycle);
	
      // comb the census
	if (census->size() > 0)
	    comb_census(*mesh, *rcontrol); 
	
      // calculate the slopes of T_electron^4
	calc_t4_slope(*mesh, *state);
    }

  // make sure a Census has been created
    Ensure (census);
    Ensure (ncentot == census->size());
}

//---------------------------------------------------------------------------//
// private member functions used in Initialize
//---------------------------------------------------------------------------//

template<class MT, class PT>
void Source_Init<MT,PT>::calc_initial_census(const MT &mesh,
					     const Opacity<MT> &opacity,
					     const Mat_State<MT> &state,
					     Rnd_Control &rcontrol, 
					     int cycle)
{
  // calculate and write the initial census source
    Require (!census);

  // make the Census 
    census = new Particle_Buffer<PT>::Census();

  // calc volume emission and surface source energies
    calc_source_energies(opacity, state, cycle);
    
  // calc radiation energy for census
    calc_ecen();

  // calc initial number of census particles
    calc_ncen_init();

  // write out the initial census
    if (ncentot > 0)
	write_initial_census(mesh, rcontrol);  
}

//---------------------------------------------------------------------------//

template<class MT, class PT>
void Source_Init<MT,PT>::calc_source_energies(const Opacity<MT> &opacity, 
					      const Mat_State<MT> &state,
					      int cycle)
{
  // calc volume emission energy per cell, total
    calc_evol(opacity, state, cycle);

  // calc surface source energy per cell, total
    calc_ess();
}

//---------------------------------------------------------------------------//

template<class MT, class PT>
void Source_Init<MT,PT>::calc_source_numbers(const Opacity<MT> &opacity, 
					     const int cycle)
{
  // iterate on number of census, surface source, and volume emission
  // particles so that all particles have nearly the same ew (i.e., 
  // variance reduction, in its most basic form).  The actual census
  // particles will be combed to give the number and weight determined
  // here.

    Insist ((evoltot+esstot+ecentot) != 0, "You must specify some source!");

    int  nptryfor = npwant;
    bool retry    = true;
    int  ntry     = 0;

    double part_per_e;
    double d_ncen;
    double d_nvol;
    double d_nss;
    int    numtot;
  // check size of existing census list
    int    censize = census->size();

    while (retry)
    {
	ntry++;
	numtot = 0;

	if (nptryfor < 1)
	    nptryfor = 1;

	part_per_e = nptryfor / (evoltot+esstot+ecentot);
	for (int cell = 1; cell <= nvol.get_Mesh().num_cells(); cell++)
	{
	  // census
	    if (ecen(cell) > 0.0  &&  censize > 0)
	    {
		d_ncen = ecen(cell) * part_per_e;
	      // * bias(cell)
		ncen(cell) = static_cast<int>(d_ncen + 0.5);
	      // try our darnedest to get at least one particle
		if (ncen(cell) == 0) 
		    ncen(cell) = static_cast<int>(d_ncen + 0.9999);
		numtot += ncen(cell);
	    }
	    else
		ncen(cell) = 0;

	  // volume emission
	    if (evol(cell) > 0)
	    {
		d_nvol = evol(cell) * part_per_e;
	      // * bias(cell)
		nvol(cell) = static_cast<int>(d_nvol + 0.5);
	      // try our darnedest to get at least one particle
		if (nvol(cell) == 0) 
		    nvol(cell) = static_cast<int>(d_nvol + 0.9999);
		numtot += nvol(cell);
	    }
	    else
		nvol(cell) = 0;

	  // surface source
	    if (ess(cell) > 0)
	    {
		d_nss = ess(cell) * part_per_e;
	      // * bias(cell)
		nss(cell) = static_cast<int>(d_nss + 0.5);
	      // try our darnedest to get at least one particle
		if (nss(cell) == 0) 
		    nss(cell) = static_cast<int>(d_nss + 0.9999);
		numtot += nss(cell);
	    }
	    else
		nss(cell) = 0;
	}

	if (numtot > npwant  &&  ntry < 100  &&  nptryfor > 1)
	    nptryfor -= (numtot - npwant);
	else
	    retry = false;
    }

  // with numbers per cell calculated, calculate ew and eloss
    ncentot = 0;
    nvoltot = 0;
    nsstot  = 0;
    if (cycle > 1) eloss_cen = 0.0;   // include eloss from initial census
    eloss_vol = 0.0;
    eloss_ss  = 0.0;

    for (int cell = 1; cell <= nvol.get_Mesh().num_cells(); cell++)
    {
      // census
	if (ncen(cell) > 0)
	{
	    ew_cen(cell) = ecen(cell) / ncen(cell);
	    ncentot += ncen(cell);
	}
	else
	    ew_cen(cell) = 0.0;

	eloss_cen += ecen(cell) - ncen(cell) * ew_cen(cell);

      // volume emission (evol is adjusted for accurate temperature update)
	if (nvol(cell) > 0)
	{
	    ew_vol(cell) = evol(cell) / nvol(cell);
	    nvoltot += nvol(cell);
	}
	else
	    ew_vol(cell) = 0.0;

	eloss_vol += evol(cell) - nvol(cell) * ew_vol(cell);
	evol(cell) = nvol(cell) * ew_vol(cell);
	if (nvol(cell) == 0)
	    evol_net(cell) = 0.0;
      // we assume here that all sampling e_loss comes from actual 
      // volume emission, not the external volume source
      //	evol_net(cell) = evol(cell) - (1.0 - opacity.get_fleck(cell)) * 
      //	    evol_ext[cell-1] * evol.get_Mesh().volume(cell) * delta_t;

      // surface source 
	if (nss(cell) > 0)
	{
	    ew_ss(cell) = ess(cell) / nss(cell);
	    nsstot += nss(cell);
	}
	else
	    ew_ss(cell) = 0.0;

	eloss_ss += ess(cell) - nss(cell) * ew_ss(cell);
    }	
}

//---------------------------------------------------------------------------//
// private member functions which calculate source parameters
//---------------------------------------------------------------------------//

template<class MT, class PT>
void Source_Init<MT,PT>::calc_evol(const Opacity<MT> &opacity,
				   const Mat_State<MT> &state, int cycle)
{
  // reset evoltot
    evoltot = 0.0;

  // calc volume source and tot volume source
  // evol_net needed for temperature update
    for (int cell = 1; cell <= evol.get_Mesh().num_cells(); cell++)
    {
      // calc cell centered volume source
	evol_net(cell) = opacity.fplanck(cell) * Global::a * Global::c *
	    pow(state.get_T(cell), 4) * evol.get_Mesh().volume(cell) * 
	    delta_t;
	evol(cell) = evol_net(cell) + 
	    evol_ext[cell-1] * (1.0 - opacity.get_fleck(cell)) *  
	    evol.get_Mesh().volume(cell) * delta_t;

      // accumulate evoltot
	evoltot += evol(cell);
    }

  // calculate evol due to external radiation source
    if (rad_s_tend > 0.0)
    {
      // calculate time duration [sh]
	double duration;
	double t_remain = rad_s_tend - (cycle - 1) * delta_t;
	if (t_remain > delta_t)
	    duration = delta_t;
	else 
	    duration = t_remain;

      // calculate radiation source energy and add to evol
	if (duration > 0.0)
	{
	  // normalize volume integral of rad_source to 1
	    double rad_volume = 0.0;
	    for (int cell = 1; cell <= evol.get_Mesh().num_cells(); cell++)
		rad_volume += rad_source[cell-1] * evol.get_Mesh().volume(cell);
	    Ensure (rad_volume > 0.0);
	    double evol_add;
	    for (int cell = 1; cell <= evol.get_Mesh().num_cells(); cell++)
	    {
		evol_add = rad_source[cell-1] / rad_volume * 
		    evol.get_Mesh().volume(cell) * duration; 
		evol(cell) += evol_add;
		evoltot    += evol_add;
	    }
	}
    }
}

//---------------------------------------------------------------------------//
// caculate the total surface source and the surface source in each cell
    
template<class MT, class PT>
void Source_Init<MT,PT>::calc_ess()
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

	  // assign source face to surface source cell
	    fss(surcells[sc]) = fss.get_Mesh().
		get_bndface(ss_pos[ss], surcells[sc]);

	  // assign energy to surface source cell
	    ess(surcells[sc]) = Global::a * Global::c * 0.25 *
		ess.get_Mesh().face_area(surcells[sc], fss(surcells[sc])) *
		pow(ss_temp[ss],4) * delta_t;

	  // accumulate esstot
	    esstot += ess(surcells[sc]);
	}
    }
}  

//---------------------------------------------------------------------------//
// calculate radiation energy in each cell and total radiation energy

template<class MT, class PT>
void Source_Init<MT,PT>::calc_ecen()
{
      // reset ecentot
    ecentot = 0.0;

  // calc census radiation energy in each cell and accumulate
    for (int cell = 1; cell <= ecen.get_Mesh().num_cells(); cell++)
    {
      // calc cell centered census radiation energy
	ecen(cell) = Global::a * ecen.get_Mesh().volume(cell) *
	    pow(rad_temp[cell-1], 4);

      // accumulate evoltot
	ecentot += ecen(cell);
    }
}

//---------------------------------------------------------------------------//
// calculate initial census particles per cell and total

template<class MT, class PT>
void Source_Init<MT,PT>::calc_ncen_init()
{
  // first guess at census particles per cell
    Insist ((evoltot+esstot+ecentot) != 0, "You must specify some source!");

    int ncenwant = static_cast<int>((ecentot) / (evoltot + esstot + ecentot) 
	* npwant);

  // particles per unit energy
    double part_per_e;

  // attempt to make all census particles have the same energy weight,
  // iterate on number of initial census particles
    bool   retry = true;
    int    ntry = 0;
    int    ncenguess = ncenwant;
    double ew;

    while (retry)
    {
      // calculate census particles per cell
	ntry++;
	ncentot   =   0;
	eloss_cen = 0.0;
	if (ncenguess < 1)
	    ncenguess = 1;
	if (ecentot > 0)
	    part_per_e = ncenguess / ecentot;
	else
	    part_per_e = 0.0;

	for (int cell = 1; cell <= ncen.get_Mesh().num_cells(); cell++)
	{
	    if (ecen(cell) > 0.0)
	    {
		double d_ncen = ecen(cell) * part_per_e + 0.5;
		ncen(cell) = static_cast<int>(d_ncen);
	      // try our darnedest to get at least one particle
		if (ncen(cell) == 0)
		    ncen(cell) = static_cast<int>(d_ncen + 0.9999);
	    }
	    else
		ncen(cell) = 0;

	    if (ncen(cell) > 0)
	    {
		ew = ecen(cell) / ncen(cell);
		ncentot += ncen(cell);
	    }
	    else
		ew = 0.0;

	    eloss_cen += ecen(cell) - ew * ncen(cell);
	}

      // check to see we haven't exceeded total particles for this cycle
	if (ncentot > ncenwant  &&  ntry < 100  &&  ncenguess > 1)
	    ncenguess -= (ncentot - ncenwant);
	else
	    retry = false;
    }
}

//---------------------------------------------------------------------------//
// write the initial census
	
template<class MT, class PT>
void Source_Init<MT,PT>::write_initial_census(const MT &mesh, 
					      Rnd_Control &rcon) 
{
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
	    double ew = ecen(cell) / ncen(cell);

	  // create Particle
	    SP<PT> particle = new PT(r, omega, ew, cell, random);

	  // write particle to census
	    census->push(particle);
	}

  // update the rn_stream constant
    RNG::rn_stream = rcon.get_num();

  // a final assertion
    Ensure (RNG::rn_stream == ncentot);
    Ensure (census->size() == ncentot);
}

//---------------------------------------------------------------------------//
// comb the census -- old way, doesn't reproduce with multiple processors
// because the combing is order-dependent.
	
template<class MT, class PT>
void Source_Init<MT,PT>::old_comb_census(const MT &mesh, Rnd_Control &rcon) 
{
  // make double sure we have census particles to comb
    Require (census->size() > 0);

  // obtain a random number stream for the comb
    Sprng random = rcon.get_rn();

  // update the rn_stream constant
    RNG::rn_stream = rcon.get_num();

  // declare and initialize the comb in each cell
    vector<double> comb(mesh.num_cells());
    for (int cell = 1; cell <= mesh.num_cells(); cell++)
    {
	comb[cell-1] = random.ran();
	comb[cell-1] = max(0.0001, min(0.9999, comb[cell-1]));
    }

  // initialize number of (new, combed) census particles per cell
    ncentot = 0;
    for (int cell = 1; cell <= mesh.num_cells(); cell++)
	ncen(cell) = 0;

  // read census particle and comb it
    double dbl_cen_part;
    double prev_comb;
    int    numcomb;
    double ecencheck;
    int    cencell;
    double cenew;

  // make new census bank to hold combed census particles
    SP<Particle_Buffer<PT>::Census> comb_census = new
	Particle_Buffer<PT>::Census();

    while (census->size())
    {
      // read census particle and get cencell, ew
	SP<PT> particle = census->top();
	census->pop();
	cencell = particle->get_cell();
	cenew   = particle->get_ew();
		
	if (ew_cen(cencell) > 0)
	{
	    dbl_cen_part = cenew / ew_cen(cencell);
	    prev_comb    = comb[cencell-1];
	    comb[cencell-1] += dbl_cen_part;
	    numcomb = static_cast<int>
		(dbl_cen_part + (prev_comb - static_cast<int>(prev_comb)));
	    ecencheck += numcomb * ew_cen(cencell);

	  // create newly combed census particles
	    if (numcomb > 0)
	    {
		particle->set_ew(ew_cen(cencell));
		comb_census->push(particle);

		if (numcomb > 1)
		    for (int nc = 1; nc <= numcomb-1; nc++)
		    {
			SP<PT> another = particle;
			Sprng nran = rcon.spawn(particle->get_random());
			another->set_random(nran);
			comb_census->push(another);
		    }
		   
	      // add up newly combed census particles
		ncen(cencell) += numcomb;
		ncentot       += numcomb;
	    }
	}
    }

    Check (census->size() == 0);
    Check (comb_census->size() == ncentot);

  // assign the newly combed census to census
    census = comb_census;

    Require (fabs(ecencheck + eloss_cen - ecentot) <= 1.0e-6 * ecentot);
}
//---------------------------------------------------------------------------//
// comb the census -- different way, numcomb depends on particle's own
// random number (thus, it is not order dependent).  perform a post-comb ew
// adjustment if there is no unsampled census energy.
	
template<class MT, class PT>
void Source_Init<MT,PT>::comb_census(const MT &mesh, Rnd_Control &rcon) 
{
  // make double sure we have census particles to comb
    Require (census->size() > 0);

  // initialize number of (new, combed) census particles per cell
    ncentot = 0;
    for (int cell = 1; cell <= mesh.num_cells(); cell++)
	ncen(cell) = 0;

    double dbl_cen_part = 0.0;
    int    numcomb = 0;
    double ecencheck = 0.0;
    int    cencell = 0;
    double cenew = 0.0;

  // add total existing census energy to eloss_cen, then take away as 
  // new census particles are combed.
    eloss_cen += (ecentot - eloss_cen);

  // make new census bank to hold combed census particles
    SP<Particle_Buffer<PT>::Census> comb_census = new
	Particle_Buffer<PT>::Census();

    while (census->size())
    {
      // read census particle and get cencell, ew
	SP<PT> particle = census->top();
	census->pop();
	cencell = particle->get_cell();
	cenew   = particle->get_ew();
	Sprng random = particle->get_random();
		
	if (ew_cen(cencell) > 0)
	{
	    dbl_cen_part = (cenew / ew_cen(cencell)) + random.ran();
	    numcomb = static_cast<int>(dbl_cen_part);


	  // create newly combed census particles
	    if (numcomb > 0)
	    {
		particle->set_ew(ew_cen(cencell));
		comb_census->push(particle);

		if (numcomb > 1)
		    for (int nc = 1; nc <= numcomb-1; nc++)
		    {
			SP<PT> another = particle;
			Sprng nran = rcon.spawn(particle->get_random());
			another->set_random(nran);
			comb_census->push(another);
		    }
		   
	      // add up newly combed census particles
		ncen(cencell) += numcomb;
		ncentot       += numcomb;

	      // subtract newly combed particles from eloss_cen
		eloss_cen -= numcomb * ew_cen(cencell);
		ecencheck += numcomb * ew_cen(cencell);
	    }
	}
    }

  // Combing is a variance reduction technique. 
  // 
  // If there is imminent loss of census energy (unsampled), we must settle 
  // for the statistical conservation of energy (to adjust the ew's with 
  // imminent loss would bias the energy).  Unfortunately, energy loss
  // propagates.  If there is no imminent loss, we may adjust the ew's for 
  // exact conservation of energy and some degradation of variance reduction.
  // The check for imminent loss may be loosened by checking on some minimal
  // energy loss instead any nonzero loss.
    bool imminent_loss = false;
    for (int cell = 1; cell <= mesh.num_cells(); cell++)
    {
	if (ncen(cell) == 0 && ecen(cell) > 0)
	    imminent_loss = true;
    }

    Check (census->size() == 0);
    Check (comb_census->size() == ncentot);

    if (imminent_loss)
    {
      // assign newly combed census to census
	census = comb_census;
    }
    else
    {
      // post-comb ew adjustment: 
      // read from comb_census, modify ew, push to census
	for (int cell = 1; cell <= mesh.num_cells(); cell++)
	{
	    if (ncen(cell) > 0)
		ew_cen(cell) = ecen(cell) / ncen(cell);
	    else
		ew_cen(cell) = 0.0;
	}

	eloss_cen += ecencheck;
	ecencheck = 0.0;

	while (comb_census->size() > 0)
	{
	    SP<PT> particle = comb_census->top();
	    comb_census->pop();
	    cencell = particle->get_cell();	
	    particle->set_ew(ew_cen(cencell));
	    census->push(particle);	  
	    ecencheck += ew_cen(cencell);
	    eloss_cen -= ew_cen(cencell);
	}

	Check (census->size() == ncentot);
	Check (comb_census->size() == 0);
    }

    Require (fabs(ecencheck + eloss_cen - ecentot) <= 1.0e-6 * ecentot);

}
//---------------------------------------------------------------------------//
// calculate slope of T_electron^4 using temporarily calc'd  edge t^4's.

template<class MT, class PT>
void Source_Init<MT,PT>::calc_t4_slope(const MT &mesh, 
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

template<class MT, class PT>
void Source_Init<MT,PT>::print(ostream &out) const
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
