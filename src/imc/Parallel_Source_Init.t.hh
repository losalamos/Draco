//----------------------------------*-C++-*----------------------------------//
// Parallel_Source_Init.t.hh
// Todd J. Urbatsch
// Mon Aug  3 09:31:56 1998
//---------------------------------------------------------------------------//
// @> Parallel_Source_Init implementation file
//---------------------------------------------------------------------------//

#include "Parallel_Source_Init.hh"
#include "Global.hh"
#include "c4/global.hh"
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cmath>

IMCSPACE

// draco components
using C4::node;
using C4::nodes;
using C4::Send;
using C4::Recv;
using Global::min;

// std components
using std::fill;
using std::endl;
using std::setw;
using std::ios;
using std::fabs;
using std::setiosflags;
using std::ofstream;

//---------------------------------------------------------------------------//
// constructors
//---------------------------------------------------------------------------//
// master node constructor; requires global source quantities

template<class MT, class PT>
template<class IT>
Parallel_Source_Init<MT,PT>::Parallel_Source_Init(SP<IT> interface, 
						  SP<MT> mesh)
    : evol(mesh), evol_net(mesh), evoltot(0), ess(mesh), fss(mesh),
      esstot(0), ecen(mesh), ecentot(0), ncen(mesh), ncentot(0), nvol(mesh),
      nss(mesh), nvoltot(0), nsstot(0), eloss_vol(0), eloss_ss(0),
      eloss_cen(0), ew_vol(mesh), ew_ss(mesh), ew_cen(mesh), volrn(mesh),
      ssrn(mesh), cenrn(mesh), t4_slope(mesh, interface->get_t4_slope()), 
      npwant(0)
{
    Require (interface);
    Require (mesh);
    Require (!census);

  // get data from MT_Interface
    evol_ext         = interface->get_evol_ext();
    rad_source       = interface->get_rad_source();
    rad_s_tend       = interface->get_rad_s_tend();
    ss_pos           = interface->get_ss_pos();
    ss_temp          = interface->get_ss_temp();
    rad_temp         = interface->get_rad_temp();
    delta_t          = interface->get_delta_t();
    elapsed_t        = interface->get_elapsed_t();
    npmax            = interface->get_npmax();
    npnom            = interface->get_npnom();
    dnpdt            = interface->get_dnpdt();
    ss_dist          = interface->get_ss_dist();
    num_global_cells = interface->get_num_global_cells();
    cycle            = interface->get_cycle();
    global_cell      = interface->get_global_cells();

  // assign the census
    census = IT::get_census();

  // some assertions
    Check (num_global_cells >= mesh->num_cells());
    Check (num_global_cells > 0);
    Check (nodes() > 0);
    Check (global_cell.size() == mesh->num_cells());
    Check (rad_temp.size() == mesh->num_cells());

  // get number of particles to transport
    npwant = min(npmax, static_cast<int>(npnom + dnpdt*elapsed_t));
    Check (npwant > 0);

  // on master processor initialize global vectors
    if (!node())
    {
      // size the cells_on_proc array
	cells_on_proc.resize(nodes());

      // resize all global source vectors
	global_ecen.resize(num_global_cells);
	global_evol.resize(num_global_cells);
	global_ess.resize(num_global_cells);
	global_ncen.resize(num_global_cells);
	global_nvol.resize(num_global_cells);
	global_nss.resize(num_global_cells);
	global_ew_cen.resize(num_global_cells);
	global_ew_vol.resize(num_global_cells);
	global_ew_ss.resize(num_global_cells);
	global_cenrn.resize(num_global_cells);
	global_volrn.resize(num_global_cells);
	global_ssrn.resize(num_global_cells);

      // initialize global data
	fill(global_ecen.begin(), global_ecen.end(), 0.0);
	fill(global_evol.begin(), global_evol.end(), 0.0);
	fill(global_ess.begin(), global_ess.end(), 0.0);
	fill(global_ncen.begin(), global_ncen.end(), 0);
	fill(global_nvol.begin(), global_nvol.end(), 0);
	fill(global_nss.begin(), global_nss.end(), 0);
	fill(global_ew_cen.begin(), global_ew_cen.end(), 0.0);
	fill(global_ew_vol.begin(), global_ew_vol.end(), 0.0);
	fill(global_ew_ss.begin(), global_ew_ss.end(), 0.0);
	fill(global_cenrn.begin(), global_cenrn.end(), 0);
	fill(global_volrn.begin(), global_volrn.end(), 0);
	fill(global_ssrn.begin(), global_ssrn.end(), 0);
	
      // initialize global variables
	global_ecentot   = 0;
	global_evoltot   = 0;
	global_esstot    = 0;
	global_eloss_cen = 0;
	global_eloss_vol = 0;
	global_eloss_ss  = 0;
	global_ncentot   = 0;
	global_nvoltot   = 0;
	global_nsstot    = 0;
    }
}

//---------------------------------------------------------------------------//
// initial census and source functions
//---------------------------------------------------------------------------//
// calculate the initial census

template<class MT, class PT> SP<typename Particle_Buffer<PT>::Census>
Parallel_Source_Init<MT,PT>::calc_initial_census(SP<MT> mesh, 
						 SP<Opacity<MT> > opacity, 
						 SP<Mat_State<MT> > state,  
						 SP<Rnd_Control> rcontrol)
{
  // calculate and write the initial census source
    Require (!census);

  // make the Census 
    census = new Particle_Buffer<PT>::Census();

  // calc volume emission and surface source energies
    calc_source_energies(*opacity, *state);
    
  // calc radiation energy for census
    calc_init_ecen();

  // Collapse energies to global vectors on the master node
    if (!node())
	recv_source_energies(*mesh);
    else if (node())
	send_source_energies(*mesh);
    else
	Check (0);

  // master node determines census numbers and ew's
    if (!node())
	calc_ncen_init();

  // Send census numbers from master node to IMC-nodes
    if (!node())
	send_census_numbers(*mesh);
    else if (node())
	recv_census_numbers(*mesh);
    else
	Check (0);

  // write out the initial census on this processor
    if (ncentot > 0)
	write_initial_census(*mesh, *rcontrol);

  // return the census
    Ensure (census);
    return census;
}

//---------------------------------------------------------------------------//
// parallel source initializer

template<class MT, class PT> SP<Source<MT> > 
Parallel_Source_Init<MT,PT>::initialize(SP<MT> mesh, 
					SP<Opacity<MT> > opacity, 
					SP<Mat_State<MT> > state, 
					SP<Rnd_Control> rcontrol,
					const Particle_Buffer<PT> &buffer)
{
  // check to make sure objects exist on each processor
    Require (mesh);
    Require (opacity);
    Require (state);
    Require (rcontrol);
    Require (census);

  // each processor calculates its own source energies (evol and ess)
    calc_source_energies(*opacity, *state);

  // each processor reads and sums up its own census
    sum_up_ecen(*mesh);
	
  // Collapse energies to global vectors on the master node
    if (!node())
	recv_source_energies(*mesh);
    else if (node())
	send_source_energies(*mesh);
    else
	Check (0);

  // master node determines source numbers with global information
    if (!node())
	calc_source_numbers(*opacity);

  // Send source numbers from master node to IMC-nodes
    if (!node())
	send_source_numbers(*mesh);
    else if (node())
	recv_source_numbers(*mesh);
    else
	Check (0);

  // comb the census
    if (census->size() > 0)
	comb_census(*mesh, *rcontrol); 
    Check(ncentot == census->size());

    // calculate T4_slope (until host works)
    // calc_t4_slope(*mesh, *state);

  // build the source
    SP<Source<MT> > source;
    source = new Source<MT>(volrn, nvol, ew_vol, t4_slope, ssrn, nss, fss,
			    ew_ss, *census, ss_dist, nvoltot, nsstot,
			    rcontrol, buffer, state);

  // return source
    return source;
}

//---------------------------------------------------------------------------//
// census operations
//---------------------------------------------------------------------------//
// write the initial census
	
template<class MT, class PT>
void Parallel_Source_Init<MT,PT>::write_initial_census(const MT &mesh, 
						       Rnd_Control &rcon) 
{
  // loop over cells
    for (int cell = 1; cell <= mesh.num_cells(); cell++)
	for (int i = 1; i <= ncen(cell); i++)
	{
	  // make a new random number for delivery to Particle
	    Sprng random = rcon.get_rn(global_cenrn[cell-1] + i - 1);
	    
	  // sample particle location
	    vector<double> r = mesh.sample_pos(cell, random);

	  // sample particle direction
	    vector<double> omega = mesh.get_Coord().
		sample_dir("isotropic", random);
	    
	  // sample frequency (not now; 1 group)
	  // ew was calculated in Parallel_Source_Init

	  // create Particle
	    SP<PT> particle = new PT(r, omega, ew_cen(cell), cell, random);

	  // write particle to census
	    census->push(particle);
	}

  // a final assertion
    Ensure (census->size() == ncentot);
}

//---------------------------------------------------------------------------//
// comb the census (reproducible)-- numcomb depends on particle's own
// random number (thus, it is not order dependent).  perform a post-comb ew
// adjustment if there is no unsampled census energy.
	
template<class MT, class PT>
void Parallel_Source_Init<MT,PT>::comb_census(const MT &mesh, 
					      Rnd_Control &rcon) 
{
  // make double sure we have census particles to comb
    Require (census->size() > 0);
    
  // initialize number of (new, combed) census particles per cell
    ncentot = 0;
    for (int cell = 1; cell <= mesh.num_cells(); cell++)
	ncen(cell) = 0;

    double dbl_cen_part = 0.0;
    int    numcomb      = 0;
    double ecencheck    = 0.0;
    int    cencell      = 0;
    double cenew        = 0.0;

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
	cencell      = particle->get_cell();
	cenew        = particle->get_ew();
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
		      // COPY a new particle and spawn a new RN state
		      	SP<PT> another = new PT(*particle);
			Sprng nran     = rcon.spawn(particle->get_random());
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
	ecencheck  = 0.0;
	
      // make a temporary census to put ew-adjusted particles into
	SP<Particle_Buffer<PT>::Census> adjusted_census = 
	    new Particle_Buffer<PT>::Census();

	while (comb_census->size() > 0)
	{
	    SP<PT> particle = comb_census->top();
	    comb_census->pop();
	    cencell = particle->get_cell();	
	    particle->set_ew(ew_cen(cencell));
	    adjusted_census->push(particle);	  
	    ecencheck += ew_cen(cencell);
	    eloss_cen -= ew_cen(cencell);
	}
	
	census = adjusted_census;
	Check (census->size() == ncentot);
	Check (comb_census->size() == 0);
    }

    Ensure (fabs(ecencheck + eloss_cen - ecentot) <= 1.0e-6 * ecentot);
}

//---------------------------------------------------------------------------//
// general source functions
//---------------------------------------------------------------------------//
// calculate number of source particles for each type of source

template<class MT, class PT> void 
Parallel_Source_Init<MT,PT>::calc_source_numbers(const Opacity<MT> &opacity)
{
  // iterate on global numbers of census, surface source, and volume emission
  // particles so that all particles have nearly the same ew (i.e., 
  // variance reduction, in its most basic form).  The actual census
  // particles will be combed to give approximately the number and weight 
  // determined here.

  // make sure we're on the master only
    Require (!node());

  // calculate total source energy
    double global_etot = global_evoltot + global_esstot + global_ecentot;
    Insist (global_etot != 0, "You must specify some source!");

    int  nptryfor = npwant;
    bool retry    = true;
    int  ntry     = 0;

    double part_per_e;
    double d_ncen;
    double d_nvol;
    double d_nss;
    int    numtot;

  // do iteration step
    while (retry)
    {
	ntry++;
	numtot = 0;

	if (nptryfor < 1)
	    nptryfor = 1;

	part_per_e = nptryfor / global_etot;
	for (int cell = 0; cell < num_global_cells; cell++)
	{
	  // census
	    if (global_ecen[cell] > 0.0  &&  census->size() > 0)
	    {
		d_ncen = global_ecen[cell] * part_per_e;
	      // * global_bias[cell]
		global_ncen[cell] = static_cast<int>(d_ncen + 0.5);
	      // try our darnedest to get at least one particle
		if (global_ncen[cell] == 0) 
		    global_ncen[cell] = static_cast<int>(d_ncen + 0.9999);
		numtot += global_ncen[cell];
	    }
	    else
		global_ncen[cell] = 0;

	  // volume emission
	    if (global_evol[cell] > 0)
	    {
		d_nvol = global_evol[cell] * part_per_e;
	      // * global_bias[cell]
		global_nvol[cell] = static_cast<int>(d_nvol + 0.5);
	      // try our darnedest to get at least one particle
		if (global_nvol[cell] == 0) 
		    global_nvol[cell] = static_cast<int>(d_nvol + 0.9999);
		numtot += global_nvol[cell];
	    }
	    else
		global_nvol[cell] = 0;

	  // surface source
	    if (global_ess[cell] > 0)
	    {
		d_nss = global_ess[cell] * part_per_e;
	      // * global_bias[cell]
		global_nss[cell] = static_cast<int>(d_nss + 0.5);
	      // try our darnedest to get at least one particle
		if (global_nss[cell] == 0) 
		    global_nss[cell] = static_cast<int>(d_nss + 0.9999);
		numtot += global_nss[cell];
	    }
	    else
		global_nss[cell] = 0;
	}

      // convergence condition
	if (numtot > npwant  &&  ntry < 100  &&  nptryfor > 1)
	    nptryfor -= (numtot - npwant);
	else
	    retry = false;
    }

  // with numbers per cell calculated, calculate ew and eloss.
    global_ncentot   = 0;
    global_nvoltot   = 0;
    global_nsstot    = 0;
    global_eloss_cen = 0.0; 
    global_eloss_vol = 0.0;
    global_eloss_ss  = 0.0;

    for (int cell = 0; cell < num_global_cells; cell++)
    {
      // census
	if (global_ncen[cell] > 0)
	{
	    global_ew_cen[cell] = global_ecen[cell] / global_ncen[cell];
	    global_ncentot     += global_ncen[cell];
	}
	else
	    global_ew_cen[cell] = 0.0;

	global_eloss_cen += global_ecen[cell] - 
	    global_ncen[cell] * global_ew_cen[cell];

      // volume emission (evol is adjusted for accurate temperature update)
	if (global_nvol[cell] > 0)
	{
	    global_ew_vol[cell] = global_evol[cell] / global_nvol[cell];
	    global_nvoltot     += global_nvol[cell];
	}
	else
	    global_ew_vol[cell] = 0.0;

	global_eloss_vol += global_evol[cell] - 
	    global_nvol[cell] * global_ew_vol[cell];
	global_evol[cell] = global_nvol[cell] * global_ew_vol[cell];

      // surface source 
	if (global_nss[cell] > 0)
	{
	    global_ew_ss[cell] = global_ess[cell] / global_nss[cell];
	    global_nsstot     += global_nss[cell];
	}
	else
	    global_ew_ss[cell] = 0.0;

	global_eloss_ss += global_ess[cell] - 
	    global_nss[cell] * global_ew_ss[cell];
    }

  // set random number stream numbers for vol, first, and ss, second.
  // update the global stream number.
    int rn_count = RNG::rn_stream;
    for (int cell = 0; cell < num_global_cells; cell++)
    {
	global_volrn[cell] = rn_count;
	rn_count          += global_nvol[cell];
    }

    RNG::rn_stream += global_nvoltot;
    Check (rn_count == RNG::rn_stream);

    for (int cell = 0; cell < num_global_cells; cell++)
    {
	global_ssrn[cell]  = rn_count;
	rn_count          += global_nss[cell];
    }

    RNG::rn_stream += global_nsstot;
    Ensure (rn_count == RNG::rn_stream);
}

//---------------------------------------------------------------------------//
// calculate initial census particles per cell, ew, and random number stream

template<class MT, class PT>
void Parallel_Source_Init<MT,PT>::calc_ncen_init()
{
  // first guess at census particles per cell
  // done only on master node on zeroth cycle
    Check (!node());
    double global_etot = global_evoltot + global_esstot + global_ecentot;
    Insist (global_etot != 0, "You must specify some source!");

    int ncenwant = static_cast<int>((global_ecentot / global_etot) * npwant);

  // particles per unit energy
    double part_per_e;

  // attempt to make all census particles have the same energy weight,
  // iterate on number of initial census particles
    bool   retry = true;
    int    ntry  = 0;
    int    ncenguess = ncenwant;
    double ew;

    while (retry)
    {
      // calculate census particles per cell
	ntry++;
	global_ncentot   =   0;
	global_eloss_cen = 0.0;
	if (ncenguess < 1)
	    ncenguess = 1;
	if (global_ecentot > 0)
	    part_per_e = ncenguess / global_ecentot;
	else
	    part_per_e = 0.0;

	for (int cell = 0; cell < num_global_cells; cell++)
	{
	    if (global_ecen[cell] > 0.0)
	    {
		double d_ncen     = global_ecen[cell] * part_per_e + 0.5;
		global_ncen[cell] = static_cast<int>(d_ncen);
	      // try our darnedest to get at least one particle
		if (global_ncen[cell] == 0)
		    global_ncen[cell] = static_cast<int>(d_ncen + 0.9999);
	    }
	    else
		global_ncen[cell] = 0;

	    if (global_ncen[cell] > 0)
	    {
		ew = global_ecen[cell] / global_ncen[cell];
		global_ncentot += global_ncen[cell];
	    }
	    else
		ew = 0.0;

	    global_eloss_cen += global_ecen[cell] - ew * global_ncen[cell];
	}

      // check to see we haven't exceeded total particles for this cycle
	if (global_ncentot > ncenwant  &&  ntry < 100  &&  ncenguess > 1)
	    ncenguess -= (global_ncentot - ncenwant);
	else
	    retry = false;
    }

  // set census ew and random number stream number per global cell
    Require (RNG::rn_stream == 0);
    int rn_count = RNG::rn_stream;

    for (int cell = 0; cell < num_global_cells; cell++)
    {
	if (global_ncen[cell] > 0.0)
	    global_ew_cen[cell] = global_ecen[cell] / global_ncen[cell];
	else
	    global_ew_cen[cell] = 0.0;

	global_cenrn[cell] = rn_count;
	rn_count          += global_ncen[cell];
    }

  // update global random number stream number
    Ensure (rn_count == global_ncentot);
    RNG::rn_stream += global_ncentot;
    Ensure (RNG::rn_stream == global_ncentot);
}

//---------------------------------------------------------------------------//
// calculate volume and surface source energies

template<class MT, class PT> void 
Parallel_Source_Init<MT,PT>::calc_source_energies(const Opacity<MT> &opacity, 
						  const Mat_State<MT> &state)
{
  // calc volume emission energy per cell, total
    calc_evol(opacity, state);

  // calc surface source energy per cell, total
    calc_ess();
}

//---------------------------------------------------------------------------//
// private member functions which calculate source parameters
//---------------------------------------------------------------------------//
// calculate total volume emission, which is the sum of the fleck-explicit
// volume emission, the fleck-implicit emission due to a material energy
// source density, and "emission" due to a radiation energy source (which, 
// on input, MUST be normalized to unity over space for Su/Olson benchmark).

template<class MT, class PT>
void Parallel_Source_Init<MT,PT>::calc_evol(const Opacity<MT> &opacity,
					    const Mat_State<MT> &state)
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
	double t_remain = rad_s_tend - elapsed_t;
	if (t_remain > delta_t)
	    duration = delta_t;
	else 
	    duration = t_remain;

      // calculate radiation source energy and add to evol
	if (duration > 0.0)
	{
	  // ASSUME rad_source normalized to 1 over space
	    double evol_add;
	    for (int cell = 1; cell <= evol.get_Mesh().num_cells(); cell++)
	    {
		evol_add = rad_source[cell-1] * 
		    evol.get_Mesh().volume(cell) * duration; 
		evol(cell) += evol_add;
		evoltot    += evol_add;
	    }
	}
    }
}

//---------------------------------------------------------------------------//
// calculate the total surface source and the surface source in each cell
    
template<class MT, class PT>
void Parallel_Source_Init<MT,PT>::calc_ess()
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
// calculate the initial census energy on each processor

template<class MT, class PT>
void Parallel_Source_Init<MT,PT>::calc_init_ecen()
{
  // reset ecentot on this processor
    ecentot = 0.0;

  // calc census radiation energy in each cell and accumulate
    for (int cell = 1; cell <= ecen.get_Mesh().num_cells(); cell++)
    {
      // calc cell centered census radiation energy
	ecen(cell) = Global::a * ecen.get_Mesh().volume(cell) *
	    pow(rad_temp[cell-1], 4);

      // accumulate ecentot
	ecentot += ecen(cell);
    }
}

//---------------------------------------------------------------------------//
// accumulate the census energy on each processor

template<class MT, class PT>
void Parallel_Source_Init<MT,PT>::sum_up_ecen(const MT &mesh)
{
    Require (census);
    ncentot = census->size();

  // initialize census energy per cell and total
    for (int cell = 1; cell <= mesh.num_cells(); cell++)
	ecen(cell) = 0.0;

    ecentot = 0.0;

  // read each census particle, get its cell (proc-local cell index), 
  // and accumulate its ew to ecen(local_cell).
    int cencell;
    double cenew;
    Particle_Buffer<PT>::Census &census_ref = *census;
    for (int i = 0; i < ncentot; i++)
    {
	SP<PT> particle = census_ref[i]; 
	cencell         = particle->get_cell();
	cenew           = particle->get_ew();
	ecen(cencell)  += cenew;
	ecentot        += cenew;
    }
}

//---------------------------------------------------------------------------//
// send source energies from IMC-nodes

template<class MT, class PT>
void Parallel_Source_Init<MT,PT>::send_source_energies(const MT &mesh)
{
  // performed by IMC-nodes, make sure we're not on master node
    Check(node());

  // number of cells on this processor
    int num_cells = mesh.num_cells();

  // define c-style arrays for sending source energies and cell info
    int *global_cell_send = new int[num_cells];
    double *ecen_send     = new double[num_cells];
    double *evol_send     = new double[num_cells];
    double *ess_send      = new double[num_cells];

  // assign arrays to temporary sending arrays
    for(int cell = 1; cell <= num_cells; cell++)
    {
	global_cell_send[cell-1] = global_cell[cell-1];
	ecen_send[cell-1]        = ecen(cell);
	evol_send[cell-1]        = evol(cell);
	ess_send[cell-1]         = ess(cell);
    }

  // Send global_cell info and source_energies to master node
    Send (num_cells, 0, 400);
    Send (global_cell_send, num_cells, 0, 401);
    Send (ecen_send, num_cells, 0, 402);
    Send (evol_send, num_cells, 0, 403);
    Send (ess_send, num_cells, 0, 404);

  // reclaim storage
    delete [] global_cell_send;
    delete [] ecen_send;
    delete [] evol_send;
    delete [] ess_send;
}

//---------------------------------------------------------------------------//
// receive source energies at master
 
template<class MT, class PT>
void Parallel_Source_Init<MT,PT>::recv_source_energies(const MT &mesh)
{
  // performed by master node, make sure we're on master node
    Check(!node());

  // number of cells on master node, ncells
    int num_cells = mesh.num_cells();
    int ncells_on_proc;

  // Receive number of cells on sending-processor, build global
  // source energy vectors and topology map
    for (int p_recv = 1; p_recv < nodes(); p_recv++)
    {
      // get number of cells on proc p_recv and their global indices 
	Recv (ncells_on_proc, p_recv, 400);
	int *global_cell_recv = new int[ncells_on_proc];
	Recv (global_cell_recv, ncells_on_proc, p_recv, 401);
	cells_on_proc[p_recv].resize(ncells_on_proc);

      // receive the source energies from proc p_recv
	double *ecen_recv = new double[ncells_on_proc];
	double *evol_recv = new double[ncells_on_proc];
	double *ess_recv  = new double[ncells_on_proc];
	Recv (ecen_recv, ncells_on_proc, p_recv, 402);
	Recv (evol_recv, ncells_on_proc, p_recv, 403);
	Recv (ess_recv,  ncells_on_proc, p_recv, 404);

      // map processor p_recv's cells and energies to global vectors.
      // for replicated cells, ess and evol are not accumulated, and 
      // ecen is accumulated except on zeroth cycle.
	for (int nc = 0; nc < ncells_on_proc; nc++)
	{
	  // assign received energies to global vectors on master node
	    int gcell                 = global_cell_recv[nc];
	    cells_on_proc[p_recv][nc] = gcell;

	    if (cycle > 0)
		global_ecen[gcell-1] += ecen_recv[nc];
	    else
	    {
		if (global_ecen[gcell-1] == 0)
		    global_ecen[gcell-1] = ecen_recv[nc];
	    }

	    if (global_evol[gcell-1] == 0.0)
		global_evol[gcell-1] = evol_recv[nc];

	    if (global_ess[gcell-1]  == 0.0)
		global_ess[gcell-1]  = ess_recv[nc];
	}

      // reclaim memory
	delete [] global_cell_recv;
	delete [] ecen_recv;
	delete [] evol_recv;
	delete [] ess_recv;
    }

  // Add contributions from the master's cells
    cells_on_proc[0].resize(num_cells);
    for (int nc = 1; nc <= num_cells; nc++)
    {
	int gcell = global_cell[nc-1];
	cells_on_proc[0][nc-1] = gcell;
	if (cycle > 0)
	    global_ecen[gcell-1] += ecen(nc);
	else
	    if (global_ecen[gcell-1] == 0.0)
		global_ecen[gcell-1] = ecen(nc);
	if (global_evol[gcell-1] == 0.0)
	    global_evol[gcell-1] = evol(nc);
	if (global_ess[gcell-1] == 0.0)
	    global_ess[gcell-1] = ess(nc);
    }

  // initialize global totals
    global_ecentot = 0.0;
    global_evoltot = 0.0;
    global_esstot  = 0.0;

  // accumulate global energy totals
    for (int gc = 0; gc < num_global_cells; gc++)
    {
	global_ecentot += global_ecen[gc];
	global_evoltot += global_evol[gc];
	global_esstot  += global_ess[gc];
    }
    
    Ensure (global_ecentot + global_evoltot + global_esstot);
}

//---------------------------------------------------------------------------//
// send source numbers from master 

template<class MT, class PT>
void Parallel_Source_Init<MT,PT>::send_source_numbers(const MT &mesh)
{
  // performed by master node; make sure we're on master node
    Check(!node());

  // loop over other processors
    for (int p_send = 1; p_send < nodes(); p_send++)
    {
	int ncells_on_proc = cells_on_proc[p_send].size();

      // define c-style arrays for sending source numbers
	int *ncen_send      = new int[ncells_on_proc];
	int *nvol_send      = new int[ncells_on_proc];
	int *nss_send       = new int[ncells_on_proc];
	double *ew_cen_send = new double[ncells_on_proc];
	double *ew_vol_send = new double[ncells_on_proc];
	double *ew_ss_send  = new double[ncells_on_proc];
	int *ssrn_send      = new int[ncells_on_proc];
	int *volrn_send     = new int[ncells_on_proc];

      // assign values to c-style arrays for sending
      // ASSUMES full domain decomposition 
	for (int nc = 0; nc < ncells_on_proc; nc++)
	{
	    int gcell = cells_on_proc[p_send][nc];
	    ncen_send[nc]   = global_ncen[gcell-1];
	    nvol_send[nc]   = global_nvol[gcell-1];
	    nss_send[nc]    = global_nss[gcell-1];
	    ew_cen_send[nc] = global_ew_cen[gcell-1];
	    ew_vol_send[nc] = global_ew_vol[gcell-1];
	    ew_ss_send[nc]  = global_ew_ss[gcell-1];
	    ssrn_send[nc]   = global_ssrn[gcell-1];
	    volrn_send[nc]  = global_volrn[gcell-1];
	}

      // send source number info to IMC-processors
	Send (ncen_send,   ncells_on_proc, p_send, 440);
	Send (nvol_send,   ncells_on_proc, p_send, 441);
	Send (nss_send,    ncells_on_proc, p_send, 442);
	Send (ew_cen_send, ncells_on_proc, p_send, 443);
	Send (ew_vol_send, ncells_on_proc, p_send, 444);
	Send (ew_ss_send,  ncells_on_proc, p_send, 445);
	Send (ssrn_send,   ncells_on_proc, p_send, 446);
	Send (volrn_send,  ncells_on_proc, p_send, 447);

      // reclaim memory
	delete [] ncen_send;
	delete [] nvol_send;
	delete [] nss_send;
	delete [] ew_cen_send;
	delete [] ew_vol_send;
	delete [] ew_ss_send;
	delete [] ssrn_send;
	delete [] volrn_send;
    }

  // assign data on master processor and accumulate totals on master
    ncentot = 0;
    nvoltot = 0;
    nsstot  = 0;
    int ncells_on_proc = mesh.num_cells();
    for (int nc = 1; nc <= ncells_on_proc; nc++)
    {
      // assign global_cell data to local data
	int gcell  = cells_on_proc[0][nc-1];
	ncen(nc)   = global_ncen[gcell-1];
	nvol(nc)   = global_nvol[gcell-1];
	nss(nc)    = global_nss[gcell-1];
	ew_cen(nc) = global_ew_cen[gcell-1];
	ew_vol(nc) = global_ew_vol[gcell-1];
	ew_ss(nc)  = global_ew_ss[gcell-1];
	ssrn(nc)   = global_ssrn[gcell-1];
	volrn(nc)  = global_volrn[gcell-1];

      // accumulate totals on master
	ncentot += ncen(nc);
	nvoltot += nvol(nc);
	nsstot  += nss(nc);
    }
}

//---------------------------------------------------------------------------//
// receive source numbers at IMC-nodes

template<class MT, class PT>
void Parallel_Source_Init<MT,PT>::recv_source_numbers(const MT &mesh)
{
  // performed by IMC nodes; make sure we're not on master node
    Check (node());
      
  // number of cells on this processor
    int num_cells = mesh.num_cells();

  // define c-style arrays for receiving source numbers
    int *ncen_recv      = new int[num_cells];
    int *nvol_recv      = new int[num_cells];
    int *nss_recv       = new int[num_cells];
    double *ew_cen_recv = new double[num_cells];
    double *ew_vol_recv = new double[num_cells];
    double *ew_ss_recv  = new double[num_cells];
    int *ssrn_recv      = new int[num_cells];
    int *volrn_recv     = new int[num_cells];

  // receive source number info from master
    Recv (ncen_recv,   num_cells, 0, 440);
    Recv (nvol_recv,   num_cells, 0, 441);
    Recv (nss_recv,    num_cells, 0, 442);
    Recv (ew_cen_recv, num_cells, 0, 443);
    Recv (ew_vol_recv, num_cells, 0, 444);
    Recv (ew_ss_recv,  num_cells, 0, 445);
    Recv (ssrn_recv,   num_cells, 0, 446);
    Recv (volrn_recv,  num_cells, 0, 447);

  // assign to processor's cell-centered scalar fields
    for (int cell = 1; cell <= num_cells; cell++)
    {
      // assign received data
	ncen(cell)   = ncen_recv[cell-1];
	nvol(cell)   = nvol_recv[cell-1];
	nss(cell)    = nss_recv[cell-1];
	ew_cen(cell) = ew_cen_recv[cell-1];
	ew_vol(cell) = ew_vol_recv[cell-1];
	ew_ss(cell)  = ew_ss_recv[cell-1];
	ssrn(cell)   = ssrn_recv[cell-1];
	volrn(cell)  = volrn_recv[cell-1];

      // zero out evol_net if there are no volume emission particles in the
      // cell, this is correct ONLY for full DD
	if (nvol(cell) == 0)
	    evol_net(cell) = 0.0;
    }

  // reclaim memory
    delete [] ncen_recv;
    delete [] nvol_recv;
    delete [] nss_recv;
    delete [] ew_cen_recv;
    delete [] ew_vol_recv;
    delete [] ew_ss_recv;
    delete [] ssrn_recv;
    delete [] volrn_recv;

  // accumulate totals on this processor
    ncentot = 0;
    nvoltot = 0;
    nsstot  = 0;
    for (int cell = 1; cell <= num_cells; cell++)
    {
	ncentot += ncen(cell);
	nvoltot += nvol(cell);
	nsstot  += nss(cell);
    }
}

//---------------------------------------------------------------------------//
// receive census numbers at IMC-nodes

template<class MT, class PT>
void Parallel_Source_Init<MT,PT>::recv_census_numbers(const MT &mesh)
{
  // performed by IMC nodes; make sure we're not on master node
    Check (node());
      
  // number of cells on this processor
    int num_cells = mesh.num_cells();

  // define c-style arrays for receiving source numbers
    int *ncen_recv      = new int[num_cells];
    double *ew_cen_recv = new double[num_cells];
    int *cenrn_recv     = new int[num_cells];

  // receive source number info from master
    Recv (ncen_recv,   num_cells, 0, 460);
    Recv (ew_cen_recv, num_cells, 0, 461);
    Recv (cenrn_recv,  num_cells, 0, 462);

  // assign to processor's cell-centered scalar fields
    for (int cell = 1; cell <= num_cells; cell++)
    {
	ncen(cell)   = ncen_recv[cell-1];
	ew_cen(cell) = ew_cen_recv[cell-1];
	cenrn(cell)  = cenrn_recv[cell-1];
    }

  // reclaim memory
    delete [] ncen_recv;
    delete [] ew_cen_recv;
    delete [] cenrn_recv;

  // accumulate totals on this processor
    ncentot = 0;
    for (int cell = 1; cell <= num_cells; cell++)
	ncentot += ncen(cell);
}

//---------------------------------------------------------------------------//
// send census numbers from master 

template<class MT, class PT>
void Parallel_Source_Init<MT,PT>::send_census_numbers(const MT &mesh)
{
  // performed by master node; make sure we're on master node
    Check(!node());

  // loop over other processors
    for (int p_send = 1; p_send < nodes(); p_send++)
    {
	int ncells_on_proc = cells_on_proc[p_send].size();

      // define c-style arrays for sending source numbers
	int *ncen_send      = new int[ncells_on_proc];
	double *ew_cen_send = new double[ncells_on_proc];
	int *cenrn_send     = new int[ncells_on_proc];

      // assign values to c-style arrays for sending
      // ASSUMES full domain decomposition 
	for (int nc = 0; nc < ncells_on_proc; nc++)
	{
	    int gcell       = cells_on_proc[p_send][nc];
	    ncen_send[nc]   = global_ncen[gcell-1];
	    ew_cen_send[nc] = global_ew_cen[gcell-1];
	    cenrn_send[nc]  = global_cenrn[gcell-1];
	}

      // send source number info to IMC-processors
	Send (ncen_send,   ncells_on_proc, p_send, 440);
	Send (ew_cen_send, ncells_on_proc, p_send, 443);
	Send (cenrn_send,  ncells_on_proc, p_send, 447);

      // reclaim memory
	delete [] ncen_send;
	delete [] ew_cen_send;
	delete [] cenrn_send;
    }

  // assign data on master processor
    ncentot = 0;
    int ncells_on_proc = mesh.num_cells();
    Check (cells_on_proc[0].size() == ncells_on_proc);
    for (int nc = 1; nc <= ncells_on_proc; nc++)
    {
      // map global census data to master processor cells
	int gcell  = cells_on_proc[0][nc-1];
	ncen(nc)   = global_ncen[gcell-1];
	ew_cen(nc) = global_ew_cen[gcell-1];
	cenrn(nc)  = global_cenrn[gcell-1];

      // accumulate census totals
	ncentot += ncen(nc);
    }
}
//---------------------------------------------------------------------------//
// calculate slope of T_electron^4 using temporarily calc'd  edge t^4's.

template<class MT, class PT>
void Parallel_Source_Init<MT,PT>::calc_t4_slope(const MT &mesh, 
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
// diagnostic functions for Parallel_Source_Init
//---------------------------------------------------------------------------//
// print out the Source_Initialization GLOBAL data

template<class MT, class PT>
void Parallel_Source_Init<MT,PT>::print(ostream &out) const
{
    out << ">>> PARALLEL SOURCE INITIALIZATION <<<" << endl;
    out << "======================================" << endl;

  // give them the particulars of the source init
    out << setw(35) << setiosflags(ios::right) 
	<< "Number of particles requested: " << setw(10) << npnom << endl;
    out << setw(35) << setiosflags(ios::right)
	<< "Total number calculated: " << setw(10) 
	<< global_ncentot + global_nvoltot + global_nsstot << endl;
    out << " ** Breakdown ** " << endl;
    out << setw(28) << "Census Particles (est): " << setw(10)
	<< global_ncentot << endl;
    out << setw(28) << "Volume Particles: " << setw(10)
	<< global_nvoltot << endl;
    out << setw(28) << "Surface Particles: " << setw(10)
	<< global_nsstot << endl;

    out << endl << " ** Source Energies ** " << endl;
    out.precision(4);
    out.setf(ios::scientific, ios::floatfield);
    out << setw(10) << setiosflags(ios::right) << "Cell"
        << setw(15) << setiosflags(ios::right) << "Volume ew"
        << setw(15) << setiosflags(ios::right) << "Surface ew"  
	<< setw(15) << setiosflags(ios::right) << "Census ew" << endl;
    for (int i = 1; i <= num_global_cells; i++)
        out << setw(10) << i << setw(15) << global_ew_vol[i-1] << setw(15)
            << global_ew_ss[i-1] << setw(15) << global_ew_cen[i-1] << endl;	

    out << "======================================" << endl;
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of Parallel_Source_Init.t.hh
//---------------------------------------------------------------------------//
