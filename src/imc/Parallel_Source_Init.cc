//----------------------------------*-C++-*----------------------------------//
// Parallel_Source_Init.cc
// Todd J. Urbatsch
// Mon Aug  3 09:31:56 1998
//---------------------------------------------------------------------------//
// @> Parallel_Source_Init implementation file
//---------------------------------------------------------------------------//

#include "imc/Parallel_Source_Init.hh"
#include "imc/Global.hh"
#include "c4/global.hh"
#include <iostream>
#include <algorithm>
#include <cmath>

IMCSPACE

// draco components
using C4::node;
using C4::nodes;
using C4::Send;
using C4::Recv;

// std components
using std::fill;
using std::endl;
using std::setw;
using std::ios;
using std::fabs;

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
      eloss_cen(0), ew_vol(mesh), ew_ss(mesh), ew_cen(mesh)
{
    Require (interface);
    Require (mesh);

  // get data from MT_Interface
    evol_ext         = interface->get_evol_ext();
    rad_source       = interface->get_rad_source();
    rad_s_tend       = interface->get_rad_s_tend();
    ss_pos           = interface->get_ss_pos();
    ss_temp          = interface->get_ss_temp();
    delta_t          = interface->get_delta_t();
    elapsed_t        = interface->get_elapsed_t();
    npmax            = interface->get_npmax();
    npnom            = interface->get_npnom();
    dnpdt            = interface->get_dnpdt();
    ss_dist          = interface->get_ss_dist();
    num_global_cells = interface->get_num_global_cells();

  // some assertions
    Check (num_global_cells >= mesh->num_cells());
    Check (num_global_cells > 0);
    Check (nodes() > 0);

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
	
      // initialize global variables
	global_ecentot   = 0;
	global_evoltot   = 0;
	global_esstot    = 0;
	global_eloss_cen = 0;
	global_eloss_vol = 0;
	global_eloss_ss  = 0;
    }
}

//---------------------------------------------------------------------------//
initial census and source functions
//---------------------------------------------------------------------------//
calculate the initial census

template<class MT, class PT> void 
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
    calc_source_energies(opacity, state);
    
  // calc radiation energy for census
    calc_initial_ecen();

  // Collapse energies to global vectors on the master node
    if (!node())
	recv_source_energies();
    else if (node())
	send_source_energies();
    else
	Check (0);

  // master node determines census numbers and ew's
    if (!node())
	calc_ncen_init();

  // Send census numbers from master node to IMC-nodes
    if (!node())
	send_census_numbers();
    else if (node())
	recv_census_numbers();
    else
	Check (0);

  // write out the initial census on this processor
    if (ncentot > 0)
	write_initial_census(mesh, rcontrol);  
}

// //---------------------------------------------------------------------------//
// // parallel source initializer

// template<class MT, class PT>
// void Parallel_Source_Init<MT,PT>::initialize(SP<MT> mesh, 
// 					     SP<Opacity<MT> > opacity, 
// 					     SP<Mat_State<MT> > state, 
// 					     SP<Rnd_Control> rcontrol, 
// 					     int cycle, double elapsed_t)
// {
//   // check to make sure objects exist on each processor
//     Require (mesh);
//     Require (opacity);
//     Require (state);
//     Require (rcontrol);
//     Require (census);

//   // each processor calculates its own source energies (evol and ess)
//     calc_source_energies(*opacity, *state, elapsed_t);

//   // each processor reads and sums up its own census
//     sum_up_census();
	
//   // Collapse energies to global vectors on the master node
//     if (!node())
// 	recv_source_energies();
//     else if (node())
// 	send_source_energies();
//     else
// 	Check (0);

//   // master node determines source numbers with global information
//     if (!node())
// 	calc_source_numbers(*opacity, cycle, elapsed_t);

//   // Send source numbers from master node to IMC-nodes
//     if (!node())
// 	send_source_numbers();
//     else if (node())
// 	recv_source_numbers();
//     else
// 	Check (0);

//   // comb the census
//     if (census->size() > 0)
// 	comb_census(*mesh, *rcontrol); 

//     Ensure (ncentot == census->size());
// }

// //---------------------------------------------------------------------------//
// // write the initial census on each processor
	
// template<class MT, class PT>
// void Parallel_Source_Init<MT,PT>::write_initial_census(const MT &mesh, 
// 						       Rnd_Control &rcon) 
// {
//   // loop over cells
//     for (int cell = 1; cell <= mesh.num_cells(); cell++)
// 	for (int i = 1; i <= ncen(cell); i++)
// 	{
// 	  // make a new random number for delivery to Particle
// 	    Sprng random = rcon.get_rn();
// 	    rcon->set_num(global_cenrn[cell-1] + i - 1);
	    
// 	  // sample particle location
// 	    vector<double> r = mesh.sample_pos(cell, random);

// 	  // sample particle direction
// 	    vector<double> omega = mesh.get_Coord().
// 		sample_dir("isotropic", random);
	    
// 	  // sample frequency (not now; 1 group)
// 	  // ew was calculated in Parallel_Source_Init

// 	  // create Particle
// 	    SP<PT> particle = new PT(r, omega, ew_cen(cell), cell, random);

// 	  // write particle to census
// 	    census->push(particle);

//   // a final assertion
//     Ensure (census->size() == ncentot);
// }
// //---------------------------------------------------------------------------//
// // comb the census (reproducible)-- numcomb depends on particle's own
// // random number (thus, it is not order dependent).  perform a post-comb ew
// // adjustment if there is no unsampled census energy.
	
// template<class MT, class PT>
// void Parallel_Source_Init<MT,PT>::comb_census(const MT &mesh, 
// 					      Rnd_Control &rcon) 
// {
//   // make double sure we have census particles to comb
//     Require (census->size() > 0);

//   // initialize number of (new, combed) census particles per cell
//     ncentot = 0;
//     for (int cell = 1; cell <= mesh.num_cells(); cell++)
// 	ncen(cell) = 0;

//     double dbl_cen_part = 0.0;
//     int    numcomb      = 0;
//     double ecencheck    = 0.0;
//     int    cencell      = 0;
//     double cenew        = 0.0;

//   // add total existing census energy to eloss_cen, then take away as 
//   // new census particles are combed.
//     eloss_cen += (ecentot - eloss_cen);

//   // make new census bank to hold combed census particles
//     SP<Particle_Buffer<PT>::Census> comb_census = new
// 	Particle_Buffer<PT>::Census();

//     while (census->size())
//     {
//       // read census particle and get cencell, ew
// 	SP<PT> particle = census->top();
// 	census->pop();
// 	cencell      = particle->get_cell();
// 	cenew        = particle->get_ew();
// 	Sprng random = particle->get_random();
		
// 	if (ew_cen(cencell) > 0)
// 	{
// 	    dbl_cen_part = (cenew / ew_cen(cencell)) + random.ran();
// 	    numcomb = static_cast<int>(dbl_cen_part);


// 	  // create newly combed census particles
// 	    if (numcomb > 0)
// 	    {
// 		particle->set_ew(ew_cen(cencell));
// 		comb_census->push(particle);

// 		if (numcomb > 1)
// 		    for (int nc = 1; nc <= numcomb-1; nc++)
// 		    {
// 			SP<PT> another = particle;
// 			Sprng nran     = rcon.spawn(particle->get_random());
// 			another->set_random(nran);
// 			comb_census->push(another);
// 		    }
		   
// 	      // add up newly combed census particles
// 		ncen(cencell) += numcomb;
// 		ncentot       += numcomb;

// 	      // subtract newly combed particles from eloss_cen
// 		eloss_cen -= numcomb * ew_cen(cencell);
// 		ecencheck += numcomb * ew_cen(cencell);
// 	    }
// 	}
//     }

//   // Combing is a variance reduction technique. 
//   // 
//   // If there is imminent loss of census energy (unsampled), we must settle 
//   // for the statistical conservation of energy (to adjust the ew's with 
//   // imminent loss would bias the energy).  Unfortunately, energy loss
//   // propagates.  If there is no imminent loss, we may adjust the ew's for 
//   // exact conservation of energy and some degradation of variance reduction.
//   // The check for imminent loss may be loosened by checking on some minimal
//   // energy loss instead any nonzero loss.
//     bool imminent_loss = false;
//     for (int cell = 1; cell <= mesh.num_cells(); cell++)
//     {
// 	if (ncen(cell) == 0 && ecen(cell) > 0)
// 	    imminent_loss = true;
//     }

//     Check (census->size() == 0);
//     Check (comb_census->size() == ncentot);

//     if (imminent_loss)
//     {
//       // assign newly combed census to census
// 	census = comb_census;
//     }
//     else
//     {
//       // post-comb ew adjustment: 
//       // read from comb_census, modify ew, push to census
// 	for (int cell = 1; cell <= mesh.num_cells(); cell++)
// 	{
// 	    if (ncen(cell) > 0)
// 		ew_cen(cell) = ecen(cell) / ncen(cell);
// 	    else
// 		ew_cen(cell) = 0.0;
// 	}

// 	eloss_cen += ecencheck;
// 	ecencheck  = 0.0;

// 	while (comb_census->size() > 0)
// 	{
// 	    SP<PT> particle = comb_census->top();
// 	    comb_census->pop();
// 	    cencell = particle->get_cell();	
// 	    particle->set_ew(ew_cen(cencell));
// 	    census->push(particle);	  
// 	    ecencheck += ew_cen(cencell);
// 	    eloss_cen -= ew_cen(cencell);
// 	}

// 	Check (census->size() == ncentot);
// 	Check (comb_census->size() == 0);
//     }

//     Require (fabs(ecencheck + eloss_cen - ecentot) <= 1.0e-6 * ecentot);
// }

// //---------------------------------------------------------------------------//

// template<class MT, class PT>
// void Parallel_Source_Init<MT,PT>::calc_source_numbers(const Opacity<MT> 
// 						         &opacity, 
// 						      const int cycle, 
// 						      const double elapsed_t)
// {
//   // iterate on global numbers of census, surface source, and volume emission
//   // particles so that all particles have nearly the same ew (i.e., 
//   // variance reduction, in its most basic form).  The actual census
//   // particles will be combed to give approximately the number and weight 
//   // determined here.

//   // make sure we're on the master only
//     Check(!node());

//   // calculate total source energy
//     double global_etot = global_evoltot + global_esstot + global_ecentot;
//     Insist (global_etot != 0, "You must specify some source!");

//   // calculate number of particles this cycle
//     npwant = min(npmax, static_cast<int>(npnom + dnpdt*elapsed_t));
//     Check (npwant != 0);

//     int  nptryfor = npwant;
//     bool retry    = true;
//     int  ntry     = 0;

//     double part_per_e;
//     double d_ncen;
//     double d_nvol;
//     double d_nss;
//     int    numtot;
//     int global_num_cells = nvol.get_Mesh().get_num_global_cells();


//     while (retry)
//     {
// 	ntry++;
// 	numtot = 0;

// 	if (nptryfor < 1)
// 	    nptryfor = 1;

// 	part_per_e = nptryfor / global_etot;
// 	for (int cell = 0; cell < global_num_cells ; cell++)
// 	{
// 	  // census
// 	    if (global_ecen[cell] > 0.0  &&  censize > 0)
// 	    {
// 		d_ncen = global_ecen[cell] * part_per_e;
// 	      // * global_bias[cell]
// 		global_ncen[cell] = static_cast<int>(d_ncen + 0.5);
// 	      // try our darnedest to get at least one particle
// 		if (global_ncen[cell] == 0) 
// 		    global_ncen[cell] = static_cast<int>(d_ncen + 0.9999);
// 		numtot += global_ncen[cell];
// 	    }
// 	    else
// 		global_ncen[cell] = 0;

// 	  // volume emission
// 	    if (global_evol[cell] > 0)
// 	    {
// 		d_nvol = global_evol[cell] * part_per_e;
// 	      // * global_bias[cell]
// 		global_nvol[cell] = static_cast<int>(d_nvol + 0.5);
// 	      // try our darnedest to get at least one particle
// 		if (global_nvol[cell] == 0) 
// 		    global_nvol[cell] = static_cast<int>(d_nvol + 0.9999);
// 		numtot += global_nvol[cell];
// 	    }
// 	    else
// 		global_nvol[cell] = 0;

// 	  // surface source
// 	    if (global_ess[cell] > 0)
// 	    {
// 		d_nss = global_ess[cell] * part_per_e;
// 	      // * global_bias[cell]
// 		global_nss[cell] = static_cast<int>(d_nss + 0.5);
// 	      // try our darnedest to get at least one particle
// 		if (global_nss[cell] == 0) 
// 		    global_nss[cell] = static_cast<int>(d_nss + 0.9999);
// 		numtot +=global_nss[cell];
// 	    }
// 	    else
// 		global_nss[cell] = 0;
// 	}

// 	if (numtot > npwant  &&  ntry < 100  &&  nptryfor > 1)
// 	    nptryfor -= (numtot - npwant);
// 	else
// 	    retry = false;
//     }

//   // with numbers per cell calculated, calculate ew and eloss.
//     global_ncentot   = 0;
//     global_nvoltot   = 0;
//     global_nsstot    = 0;
//     global_eloss_cen = 0.0; 
//     global_eloss_vol = 0.0;
//     global_eloss_ss  = 0.0;

//     for (int cell = 0; cell < global_num_cells; cell++)
//     {
//       // census
// 	if (global_ncen[cell] > 0)
// 	{
// 	    global_ew_cen[cell] = global_ecen[cell] / global_ncen[cell];
// 	    global_ncentot     += global_ncen[cell];
// 	}
// 	else
// 	    global_ew_cen[cell] = 0.0;

// 	global_eloss_cen += global_ecen[cell] - 
// 	    global_ncen[cell] * global_ew_cen[cell];

//       // volume emission (evol is adjusted for accurate temperature update)
// 	if (global_nvol[cell] > 0)
// 	{
// 	    global_ew_vol[cell] = global_evol[cell] / global_nvol[cell];
// 	    global_nvoltot     += global_nvol[cell];
// 	}
// 	else
// 	    global_ew_vol[cell] = 0.0;

// 	global_eloss_vol += global_evol[cell] - 
// 	    global_nvol[cell] * global_ew_vol[cell];
// 	global_evol[cell] = global_nvol[cell] * global_ew_vol[cell];
// 	if (global_nvol[cell] == 0)
// 	    global_evol_net[cell] = 0.0;

//       // surface source 
// 	if (global_nss[cell] > 0)
// 	{
// 	    global_ew_ss[cell] = global_ess[cell] / global_nss[cell];
// 	    global_nsstot     += global_nss[cell];
// 	}
// 	else
// 	    global_ew_ss[cell] = 0.0;

// 	global_eloss_ss += global_ess[cell] - 
// 	    global_nss[cell] * global_ew_ss[cell];
//     }

//   // set random number stream numbers for vol, first, and ss, second.
//   // update the global stream number.
//     int rn_count = RNG::rn_stream;
//     for (int cell = 0; cell < global_num_cells; cell++)
//     {
// 	global_volrn[cell] = rn_count;
// 	rn_count          += global_nvol[cell];
//     }

//     RNG::rn_stream += global_nvoltot;
//     Check (rn_count == RNG::rn_stream);

//     for (int cell = 0; cell < global_num_cells; cell++)
//     {
// 	global_ssrn[cell]  = rn_count;
// 	rn_count          += global_nss[cell];
//     }

//     RNG::rn_stream += global_nsstot;
//     Check (rn_count == RNG::rn_stream);
// }

// //---------------------------------------------------------------------------//
// // calculate initial census particles per cell, ew, and random number stream

// template<class MT, class PT>
// void Parallel_Source_Init<MT,PT>::calc_ncen_init()
// {
//   // first guess at census particles per cell
//   // done only on master node on zeroth cycle
//     double global_etot = global_evoltot + global_esstot + global_ecentot;
//     Insist (global_etot != 0, "You must specify some source!");

//     int ncenwant = static_cast<int>((global_ecentot / global_etot) * npwant);

//   // particles per unit energy
//     double part_per_e;

//   // attempt to make all census particles have the same energy weight,
//   // iterate on number of initial census particles
//     bool   retry = true;
//     int    ntry  = 0;
//     int    ncenguess = ncenwant;
//     double ew;
//     int global_num_cells = ncen.get_Mesh().get_num_global_cells();

//     while (retry)
//     {
//       // calculate census particles per cell
// 	ntry++;
// 	global_ncentot   =   0;
// 	global_eloss_cen = 0.0;
// 	if (ncenguess < 1)
// 	    ncenguess = 1;
// 	if (global_ecentot > 0)
// 	    part_per_e = ncenguess / global_ecentot;
// 	else
// 	    part_per_e = 0.0;

// 	for (int cell = 0; cell < global_num_cells; cell++)
// 	{
// 	    if (global_ecen[cell] > 0.0)
// 	    {
// 		double d_ncen     = global_ecen[cell] * part_per_e + 0.5;
// 		global_ncen[cell] = static_cast<int>(d_ncen);
// 	      // try our darnedest to get at least one particle
// 		if (global_ncen[cell] == 0)
// 		    global_ncen[cell] = static_cast<int>(d_ncen + 0.9999);
// 	    }
// 	    else
// 		global_ncen[cell] = 0;

// 	    if (global_ncen[cell] > 0)
// 	    {
// 		ew = global_ecen[cell] / global_ncen[cell];
// 		global_ncentot += global_ncen[cell];
// 	    }
// 	    else
// 		ew = 0.0;

// 	    global_eloss_cen += global_ecen[cell] - ew * global_ncen[cell];
// 	}

//       // check to see we haven't exceeded total particles for this cycle
// 	if (global_ncentot > ncenwant  &&  ntry < 100  &&  ncenguess > 1)
// 	    ncenguess -= (ncentot - ncenwant);
// 	else
// 	    retry = false;
//     }

//   // set census ew and random number stream number per global cell
//     Require (RNG::rn_stream == 0);
//     int rn_count = RNG::rn_stream;

//     for (int cell = 0; cell < global_num_cells; cell++)
//     {
// 	if (global_ecen[cell] > 0.0)
// 	    global_ew_cen[cell] = global_ecen[cell] / global_ncen[cell];
// 	else
// 	    global_ew_cen[cell] = 0.0;

// 	global_cenrn[cell] = rn_count;
// 	rn_count          += global_ncen[cell];
//     }

//   // update global random number stream number
//     Ensure (rn_count == global_ncentot);
//     RNG::rn_stream += global_ncentot;
//     Ensure (RNG::rn_stream == global_ncentot);
// }

// //---------------------------------------------------------------------------//

// template<class MT, class PT>
// void Parallel_Source_Init<MT,PT>::calc_source_energies(const Opacity<MT> 
// 						         &opacity, 
// 						       const Mat_State<MT> 
// 						         &state,
// 						       double elapsed_t)
// {
//   // calc volume emission energy per cell, total
//     calc_evol(opacity, state, elapsed_t);

//   // calc surface source energy per cell, total
//     calc_ess();
// }

// //---------------------------------------------------------------------------//
// // private member functions which calculate source parameters
// //---------------------------------------------------------------------------//

// // calculate total volume emission, which is the sum of the fleck-explicit
// // volume emission, the fleck-implicit emission due to a material energy
// // source density, and "emission" due to a radiation energy source (which, 
// // on input, MUST be normalized to unity over space for Su/Olson benchmark).

// template<class MT, class PT>
// void Parallel_Source_Init<MT,PT>::calc_evol(const Opacity<MT> &opacity,
// 					    const Mat_State<MT> &state, 
// 					    double elapsed_t)
// {
//   // reset evoltot
//     evoltot = 0.0;

//   // calc volume source and tot volume source
//   // evol_net needed for temperature update 
//     for (int cell = 1; cell <= evol.get_Mesh().num_cells(); cell++)
//     {
//       // calc cell centered volume source
// 	evol_net(cell) = opacity.fplanck(cell) * Global::a * Global::c *
// 	    pow(state.get_T(cell), 4) * evol.get_Mesh().volume(cell) * 
// 	    delta_t;
// 	evol(cell) = evol_net(cell) + 
// 	    evol_ext[cell-1] * (1.0 - opacity.get_fleck(cell)) *  
// 	    evol.get_Mesh().volume(cell) * delta_t;

//       // accumulate evoltot
// 	evoltot += evol(cell);
//     }

//   // calculate evol due to external radiation source
//     if (rad_s_tend > 0.0)
//     {
//       // calculate time duration [sh]
// 	double duration;
// 	double t_remain = rad_s_tend - elapsed_t;
// 	if (t_remain > delta_t)
// 	    duration = delta_t;
// 	else 
// 	    duration = t_remain;

//       // calculate radiation source energy and add to evol
// 	if (duration > 0.0)
// 	{
// 	  // ASSUME rad_source normalized to 1 over space
// 	    double evol_add;
// 	    for (int cell = 1; cell <= evol.get_Mesh().num_cells(); cell++)
// 	    {
// 		evol_add = rad_source[cell-1] * 
// 		    evol.get_Mesh().volume(cell) * duration; 
// 		evol(cell) += evol_add;
// 		evoltot    += evol_add;
// 	    }
// 	}
//     }
// }

// //---------------------------------------------------------------------------//
// // caculate the total surface source and the surface source in each cell
    
// template<class MT, class PT>
// void Parallel_Source_Init<MT,PT>::calc_ess()
// {
//   // reset esstot
//     esstot = 0.0;

//   // loop over surface sources in problem
//     for (int ss = 0; ss < ss_pos.size(); ss++)
//     {
// 	vector<int> surcells = ess.get_Mesh().get_surcells(ss_pos[ss]);
// 	for (int sc = 0; sc < surcells.size(); sc++)
// 	{      
// 	  // make sure this cell doesn't already have a surface source
// 	    Check (fss(surcells[sc]) == 0);

// 	  // assign source face to surface source cell
// 	    fss(surcells[sc]) = fss.get_Mesh().
// 		get_bndface(ss_pos[ss], surcells[sc]);

// 	  // assign energy to surface source cell
// 	    ess(surcells[sc]) = Global::a * Global::c * 0.25 *
// 		ess.get_Mesh().face_area(surcells[sc], fss(surcells[sc])) *
// 		pow(ss_temp[ss],4) * delta_t;

// 	  // accumulate esstot
// 	    esstot += ess(surcells[sc]);
// 	}
//     }
// }  

// //---------------------------------------------------------------------------//
// // calculate the initial census energy on each processor

// template<class MT>
// void Parallel_Source_Init<MT>::calc_init_ecen()
// {
//   // reset ecentot on this processor
//     ecentot = 0.0;

//   // calc census radiation energy in each cell and accumulate
//     for (int cell = 1; cell <= ecen.get_Mesh().num_cells(); cell++)
//     {
//       // calc cell centered census radiation energy
// 	ecen(cell) = Global::a * ecen.get_Mesh().volume(cell) *
// 	    pow(rad_temp[cell-1], 4);

//       // accumulate evoltot
// 	ecentot += ecen(cell);
//     }
// }

// //---------------------------------------------------------------------------//
// // accumulate the census energy on each processor

// template<class MT>
// void Parallel_Source_Init<MT>::sum_up_ecen()
// {
//     Require (census);
//     ncentot = census->size();

//   // initialize census energy per cell and total
//     for (int cell = 1; cell <= mesh.num_cells(); cell++)
// 	ecen(cell) = 0.0;

//     ecentot = 0.0;

//   // read each census particle, get its cell (proc-local cell index), 
//   // and accumulate its ew to ecen(local_cell).
//     for (int cell = 1; cell <= mesh.num_cells(); cell++)
//     {
// 	SP<PT> particle = census->top();
// 	census->pop();
// 	cencell = particle->get_cell();
// 	cenew   = particle->get_ew();
// 	ecen(cencell) += cenew;
//     }
// }


// //---------------------------------------------------------------------------//
// // send source energies from IMC-nodes

// template<class MT>
// void Parallel_Source_Init<MT>::send_source_energies(const MT &mesh)
// {
//   // performed by IMC-nodes, make sure we're not on master node
//     Check(node());

//   // number of cells on this processor
//     int num_cells = mesh.get_num_cells();

//   // define c-style arrays for sending source energies and cell info
//     int *global_cell_send = new int[num_cells];
//     double *ecen_send     = new double[num_cells];
//     double *evol_send     = new double[num_cells];
//     double *ess_send      = new double[num_cells];

//   // assign arrays to temporary sending arrays
//     for(int cell = 1; cell <= num_cells; cell++)
//     {
// 	global_cell_send[cell-1] = mesh.get_global_cell(cell);
// 	ecen_send[cell-1]        = ecen(cell);
// 	evol_send[cell-1]        = evol(cell);
// 	ess_send[cell-1]         = ess(cell);
//     }

//   // Send global_cell info and source_energies to master node
//     Send (num_cells, 0, 400);
//     Send (global_cell_send, num_cells, 0, 401);
//     Send (ecen_send, num_cells, 0, 402);
//     Send (evol_send, num_cells, 0, 403);
//     Send (ess_send, num_cells, 0, 404);

//   // reclaim storage
//     delete [] global_cell_send;
//     delete [] ecen_send;
//     delete [] evol_send;
//     delete [] ess_send;
// }

// //---------------------------------------------------------------------------//
// // receive source energies at master
 
// template<class MT>
// void Parallel_Source_Init<MT>::recv_source_energies(const MT &mesh, 
// 						    const int cycle)
// {
//   // performed by master node, make sure we're on master node
//     Check(!node());

//   // number of cells: on master processor, total global
//     int num_cells        = mesh.get_num_cells();
//     int global_num_cells = mesh.get_global_num_cells();
//     int ncells_on_proc;

//   // Receive number of cells on sending-processor, build global
//   // source energy vectors and topology map
//     for (int p_recv = 1; p_recv < nodes(); p_recv++)
//     {
//       // get number of cells on proc p_recv and their global indices 
// 	Recv (ncells_on_proc, p_recv, 400);
// 	int *global_cell_recv = new int[ncells_on_proc];
// 	Recv (global_cell_recv, ncells_on_proc, p_recv, 401);
// 	cells_on_proc[p_recv].resize(ncells_on_proc);

//       // receive the source energies from proc p_recv
// 	double *ecen_recv = new double[ncells_on_proc];
// 	double *evol_recv = new double[ncells_on_proc];
// 	double *ess_recv  = new double[ncells_on_proc];
// 	Recv (ecen_recv, ncells_on_proc, p_recv, 402);
// 	Recv (evol_recv, ncells_on_proc, p_recv, 403);
// 	Recv (ess_recv,  ncells_on_proc, p_recv, 404);

//       // map processor p_recv's cells and energies to global vectors.
//       // for replicated cells, ess and evol are not accumulated, and 
//       // ecen is accumulated except on zeroth cycle.
// 	for (int nc = 0; nc < ncells_on_proc; nc++)
// 	{
// 	  // assign received energies to global vectors on master node
// 	    int gcell                 = global_cell_recv[nc];
// 	    cells_on_proc[p_recv][nc] = gcell;

// 	    if (cycle > 0)
// 		global_ecen[gcell-1] += ecen_recv[nc];
// 	    else
// 	    {
// 		if (global_ecen[gcell-1] == 0)
// 		    global_ecen[gcell-1] = ecen_recv[nc];
// 	    }

// 	    if (global_evol[gcell-1] == 0.0)
// 		global_evol[gcell-1] = evol_recv[nc];

// 	    if (global_ess[gcell-1]  == 0.0)
// 		global_ess[gcell-1]  = ess_recv[nc];
// 	}

//       // reclaim memory
// 	delete [] global_cell_recv;
// 	delete [] ecen_recv;
// 	delete [] evol_recv;
// 	delete [] ess_recv;
//     }

//   // Add contributions from the master's cells
//     for (int nc = 1; nc <= num_cells; nc++)
//     {
// 	int gcell = mesh.get_global_cell(nc);
// 	cells_on_proc[0][nc-1] = gcell;
// 	if (cycle > 0)
// 	    global_ecen[gcell-1] += ecen(nc);
// 	else
// 	    if (global_ecen[gcell-1] == 0.0)
// 		global_ecen[gcell-1] = ecen(nc);
// 	if (global_evol[gcell-1] == 0.0)
// 	    global_evol[gcell-1] = evol(nc);
// 	if (global_ess[gcell-1] == 0.0)
// 	    global_ess[gcell-1] = ess(nc);
//     }

//   // initialize global totals
//     global_ecentot = 0.0;
//     global_evoltot = 0.0;
//     global_esstot  = 0.0;

//   // accumulate global energy totals
//     for (int gc = 0; gc < global_num_cells; gc++)
//     {
// 	global_ecentot += global_ecen[gc];
// 	global_evoltot += global_evol[gc];
// 	global_esstot  += global_ess[gc];
//     }
// }

// //---------------------------------------------------------------------------//
// // send source numbers from master 

// template<class MT>
// void Parallel_Source_Init<MT>::send_source_numbers(const MT &mesh)
// {
//   // performed by master node; make sure we're on master node
//     Check(!node());

//   // loop over other processors
//     for (int p_send = 1; p_send <= nodes(); p_send)
//     {
// 	int ncells_on_proc = cells_on_proc[p_send].size();

//       // define c-style arrays for sending source numbers
// 	int *ncen_send      = new[ncells_on_proc];
// 	int *nvol_send      = new[ncells_on_proc];
// 	int *nss_send       = new[ncells_on_proc];
// 	double *ew_cen_send = new[ncells_on_proc];
// 	double *ew_vol_send = new[ncells_on_proc];
// 	double *ew_ss_send  = new[ncells_on_proc];
// 	int *ssrn_send      = new[ncells_on_proc];
// 	int *volrn_send     = new[ncells_on_proc];

//       // assign values to c-style arrays for sending
//       // ASSUMES full domain decomposition 
// 	for (int nc = 0; nc < ncells_on_proc; nc++)
// 	{
// 	    int gcell = cells_on_proc[p_send][nc];
// 	    ncen_send[nc]   = global_ncen[gcell];
// 	    nvol_send[nc]   = global_nvol[gcell];
// 	    nss_send[nc]    = global_nss[gcell];
// 	    ew_cen_send[nc] = global_ew_cen[gcell];
// 	    ew_vol_send[nc] = global_ew_cen[gcell];
// 	    ew_ss_send[nc]  = global_ew_ss[gcell];
// 	    ssrn_send[nc]   = global_ssrn[gcell];
// 	    volrn_send[nc]  = global_volrn[gcell];
// 	}

//       // send source number info to IMC-processors
// 	Send (ncen_send,   ncells_on_proc, p_send, 440);
// 	Send (nvol_send,   ncells_on_proc, p_send, 441);
// 	Send (nss_send,    ncells_on_proc, p_send, 442);
// 	Send (ew_cen_send, ncells_on_proc, p_send, 443);
// 	Send (ew_vol_send, ncells_on_proc, p_send, 444);
// 	Send (ew_ss_send,  ncells_on_proc, p_send, 445);
// 	Send (ssrn_send,   ncells_on_proc, p_send, 446);
// 	Send (volrn_send,  ncells_on_proc, p_send, 447);

//       // reclaim memory
// 	delete [] ncen_send;
// 	delete [] nvol_send;
// 	delete [] nss_send;
// 	delete [] ew_cen_send;
// 	delete [] ew_vol_send;
// 	delete [] ew_ss_send;
// 	delete [] ssrn_send;
// 	delete [] volrn_send;
//     }

//   // assign data on master processor
//     int ncells_on_proc = mesh.get_num_cells();
//     for (int nc = 1; nc <= ncells_on_proc; nc++)
//     {
// 	int gcell  = cells_on_proc[0][nc-1];
// 	ncen(nc)   = global_ncen[gcell];
// 	nvol(nc)   = global_nvol[gcell];
// 	nss(nc)    = global_nss[gcell];
// 	ew_cen(nc) = global_ew_cen[gcell];
// 	ew_vol(nc) = global_ew_cen[gcell];
// 	ew_ss(nc)  = global_ew_ss[gcell];
// 	ssrn(nc)   = global_ssrn[gcell];
// 	volrn(nc)  = global_volrn[gcell];
//     }
// }
// //---------------------------------------------------------------------------//
// // receive source numbers at IMC-nodes

// template<class MT>
// void Parallel_Source_Init<MT>::recv_source_numbers(const MT &mesh)
// {
//   // performed by IMC nodes; make sure we're not on master node
//     Check (node());
      
//   // number of cells on this processor
//     int num_cells = mesh.get_num_cells();

//   // define c-style arrays for receiving source numbers
//     int *ncen_recv      = new[num_cells];
//     int *nvol_recv      = new[num_cells];
//     int *nss_recv       = new[num_cells];
//     double *ew_cen_recv = new[num_cells];
//     double *ew_vol_recv = new[num_cells];
//     double *ew_ss_recv  = new[num_cells];
//     int *ssrn_recv      = new[num_cells];
//     int *volrn_send     = new[num_cells];

//   // receive source number info from master
//     Recv (ncen_recv,   num_cells, 0, 440);
//     Recv (nvol_recv,   num_cells, 0, 441);
//     Recv (nss_recv,    num_cells, 0, 442);
//     Recv (ew_cen_recv, num_cells, 0, 443);
//     Recv (ew_vol_recv, num_cells, 0, 444);
//     Recv (ew_ss_recv,  num_cells, 0, 445);
//     Recv (ssrn_recv,   num_cells, 0, 446);
//     Recv (volrn_recv,  num_cells, 0, 447);

//   // assign to processor's cell-centered scalar fields
//     for (int cell = 1; cell <= num_cells; cell++)
//     {
// 	ncen(cell)   = ncen_recv[cell-1];
// 	nvol(cell)   = nvol_recv[cell-1];
// 	nss(cell)    = nss_recv[cell-1];
// 	ew_cen(cell) = ew_cen_recv[cell-1];
// 	ew_vol(cell) = ew_vol_recv[cell-1];
// 	ew_ss(cell)  = ew_ss_recv[cell-1];
// 	ssrn(cell)   = ssrn_recv[cell-1];
// 	volrn(cell)  = volrn_recv[cell-1];
//     }

//   // reclaim memory
//     delete [] ncen_recv;
//     delete [] nvol_recv;
//     delete [] nss_recv;
//     delete [] ew_cen_recv;
//     delete [] ew_vol_recv;
//     delete [] ew_ss_recv;
//     delete [] ssrn_recv;
//     delete [] volrn_recv;

//   // accumulate totals on this processor
//     ncentot = 0;
//     nvoltot = 0;
//     nsstot  = 0;
//     for (int cell = 1; cell <= num_cells; cell++)
//     {
// 	ncentot += ncen(cell);
// 	nvoltot += nvol(cell);
// 	nsstot  += nss(cell);
//     }
// }

// //---------------------------------------------------------------------------//
// // receive census numbers at IMC-nodes

// template<class MT>
// void Parallel_Source_Init<MT>::recv_census_numbers(const MT &mesh)
// {
//   // performed by IMC nodes; make sure we're not on master node
//     Check (node());
      
//   // number of cells on this processor
//     int num_cells = mesh.get_num_cells();

//   // define c-style arrays for receiving source numbers
//     int *ncen_recv      = new[num_cells];
//     double *ew_cen_recv = new[num_cells];
//     int *cenrn_send     = new[num_cells];

//   // receive source number info from master
//     Recv (ncen_recv,   num_cells, 0, 460);
//     Recv (ew_cen_recv, num_cells, 0, 461);
//     Recv (cenrn_recv,  num_cells, 0, 462);

//   // assign to processor's cell-centered scalar fields
//     for (int cell = 1; cell <= num_cells; cell++)
//     {
// 	ncen(cell)   = ncen_recv[cell-1];
// 	ew_cen(cell) = ew_cen_recv[cell-1];
// 	cenrn(cell)  = cenrn_recv[cell-1];
//     }

//   // reclaim memory
//     delete [] ncen_recv;
//     delete [] ew_cen_recv;
//     delete [] cenrn_recv;
// }
// //---------------------------------------------------------------------------//
// // send census numbers from master 

// template<class MT>
// void Parallel_Source_Init<MT>::send_census_numbers(const MT &mesh)
// {
//   // performed by master node; make sure we're on master node
//     Check(!node());

//   // loop over other processors
//     for (int p_send = 1; p_send <= nodes(); p_send)
//     {
// 	int ncells_on_proc = cells_on_proc[p_send].size();

//       // define c-style arrays for sending source numbers
// 	int *ncen_send      = new[ncells_on_proc];
// 	double *ew_cen_send = new[ncells_on_proc];
// 	int *cenrn_send     = new[ncells_on_proc];

//       // assign values to c-style arrays for sending
//       // ASSUMES full domain decomposition 
// 	for (int nc = 0; nc < ncells_on_proc; nc++)
// 	{
// 	    int gcell       = cells_on_proc[p_send][nc];
// 	    ncen_send[nc]   = global_ncen[gcell];
// 	    ew_cen_send[nc] = global_ew_cen[gcell];
// 	    cenrn_send[nc]  = global_cenrn[gcell];
// 	}

//       // send source number info to IMC-processors
// 	Send (ncen_send,   ncells_on_proc, p_send, 440);
// 	Send (ew_cen_send, ncells_on_proc, p_send, 443);
// 	Send (cenrn_send,  ncells_on_proc, p_send, 447);

//       // reclaim memory
// 	delete [] ncen_send;
// 	delete [] ew_cen_send;
// 	delete [] cenrn_send;
//     }

//   // assign data on master processor
//     int ncells_on_proc = mesh.get_num_cells();
//     for (int nc = 1; nc <= ncells_on_proc; nc++)
//     {
// 	int gcell  = cells_on_proc[0][nc-1];
// 	ncen(nc)   = global_ncen[gcell];
// 	ew_cen(nc) = global_ew_cen[gcell];
// 	cenrn(nc)  = global_cenrn[gcell];
//     }
// }


CSPACE

//---------------------------------------------------------------------------//
//                              end of Parallel_Source_Init.cc
//---------------------------------------------------------------------------//
