//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/DD_Source_Builder.t.hh
 * \author Todd J. Urbatsch
 * \date   Tue May  2 15:11:21 2000
 * \brief  Implementation file for DD_Source_Builder.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "DD_Source_Builder.hh"
#include "Mesh_Operations.hh"
#include "Global.hh"
#include "c4/global.hh"
#include "mc/Math.hh"
#include "mc/Parallel_Data_Operator.hh"
#include <iomanip>
#include <numeric>

namespace rtt_imc
{

using C4::nodes;
using C4::node;
using rtt_dsxx::SP;
using rtt_mc::global::soft_equiv;
using rtt_mc::Parallel_Data_Operator;

using std::cout;
using std::endl;
using std::ios;
using std::accumulate;

//---------------------------------------------------------------------------//
// constructor
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for DD_Source_Builder.
 */
template<class MT, class PT>
template<class IT>
DD_Source_Builder<MT,PT>::DD_Source_Builder(SP<IT> interface, SP_Mesh mesh, 
					    SP_Topology top)
    : Source_Builder<MT,PT>(interface, mesh, top),
    local_ncen(mesh),
    local_ncentot(0),
    local_eloss_cen(0),
    global_eloss_cen(0),
    local_nvol(mesh),
    local_nvoltot(0),
    local_nss(mesh),
    local_nsstot(0),
    local_eloss_vol(0),
    global_eloss_vol(0),
    local_eloss_ss(0),
    global_eloss_ss(0)
{ 
    // at the beginning of the timestep, random number stream ID should be
    // the same on every processor.
    Check(parallel_data_op.check_global_equiv(rtt_rng::rn_stream));
    
    // Update the persistent census energy data, unless the census does not
    // exist yet, which is the case on the first IMC cycle.
    if (census)
    {
	// a check on the total census energy
	double ecentot_check = 0.0;

	for (int cell = 1; cell <= mesh->num_cells(); cell++)
	{
	    // get local values of ecen from the interface
	    ecen(cell)     = interface->get_ecen(cell);
	    ecentot_check += ecen(cell);

	    // more updates may follow -- probably time-cumulative edits
	}

	// sum up census energy check from all processors
	C4::gsum(ecentot_check);

	// get total global census energy from the interface
	global_ecentot = interface->get_ecentot();

	// check consistency of energies and totals
	Check (rtt_mc::global::soft_equiv(global_ecentot, ecentot_check,
					  1.0e-12));
    }
}

//---------------------------------------------------------------------------//
// BUILD THE SOURCE ON THIS DD TOPOLOGY
//---------------------------------------------------------------------------//
template<class MT, class PT>
DD_Source_Builder<MT,PT>::SP_Source
DD_Source_Builder<MT,PT>::build_Source(SP_Mesh mesh,
				       SP_Mat_State state,
				       SP_Opacity opacity,
				       SP_Rnd_Control rnd_control,
				       SP_Comm_Patterns patterns)
{
    int num_cells = mesh->num_cells();
    Require(num_cells == state->num_cells());
    Require(num_cells == opacity->num_cells());

    // if we don't have a census object, build it 
    if (!census)
	calc_initial_census(mesh, state, opacity, rnd_control);
    else
    {
	// calculate the volume emission and surface source energies
	calc_source_energies(*state, *opacity);

	// globally sum up volume emission and surface source energies
	global_evoltot = evoltot;
	C4::gsum(global_evoltot);
	global_esstot = esstot;
	C4::gsum(global_esstot);

	// make sure global energies are the same on each processor
	Check(parallel_data_op.check_global_equiv(global_evoltot));
	Check(parallel_data_op.check_global_equiv(global_esstot));
	Check(parallel_data_op.check_global_equiv(global_ecentot));
    }

    // calculate the number of source particles for all source particles
    calc_source_numbers();

    // initialize data for resurrecting dead, combed-out census particles
    SP_Census dead_census(new Particle_Buffer<PT>::Census());
    ccsf_int max_dead_rand_id(mesh);

    // comb census -- even if census is empty; must update local_ncen
    double eloss_comb = 0;
    comb_census(rnd_control, local_ncen, local_ncentot, eloss_comb,
		dead_census, max_dead_rand_id);

    // calculate global values of energy loss and numbers of particles due
    // to combing the census
    double global_eloss_comb = eloss_comb;
    C4::gsum(global_eloss_comb);

    global_ncentot = local_ncentot;
    C4::gsum(global_ncentot);

    // recalculate the census particles' energy-weights
    recalc_census_ew_after_comb(mesh, max_dead_rand_id, dead_census,
				global_eloss_comb);

    // update total numbers of census particles
    global_ncentot = local_ncentot;
    C4::gsum(global_ncentot);

    // reset the census particles' energy-weights
    reset_ew_in_census(local_ncentot, global_eloss_comb, global_ecentot);

    // add energy loss from the comb to the global census energy loss
    global_eloss_cen += global_eloss_comb;
    
    // build Mesh_Operations class for source
    SP<Mesh_Operations<MT> > mesh_op
	(new Mesh_Operations<MT>(mesh, state, topology, patterns)); 

    // build the source --> the source gets a SP to the census on this
    // processor, this means that as particles are taken out of the source
    // the "master" census (in a global_state object) on processor is emptied 
    // as well
    SP_Source source(new Source<MT,PT>(volrn, local_nvol, ew_vol, ssrn,
				       local_nss, ss_face_in_cell, ew_ss, 
				       census, ss_dist, local_nvoltot,
				       local_nsstot, rnd_control, state,
				       mesh_op, topology));

    // return the source
    return source;
}




//---------------------------------------------------------------------------//
// CALCULATE INITIAL CENSUS
//---------------------------------------------------------------------------//
/*!
 * \brief Build the initial census on a full domain decomposed topology
 *
 * This function is a pure virtual function in the Source_Builder base class.
 * We make it part of the public interface to give the option of building an
 * initial census and then exiting the source building process.  This is
 * useful in rad-hydro applications where the initial census must be build on
 * cycle 0, before a hydro step, from the initial radiation temperature.  The
 * build_Source() function checks to see if a census exists and then builds
 * it if it does not.  The census is loaded into the Source_Builder base
 * class through the run-time interface.
 *
 * This function builds a local initial census on each processor.  The
 * rtt_dsxx::SP to the completed census may be accessed through the
 * Source_Builder::get_census() function.  
 *
 * The full domain decomposed (DD) version of this function differs from the
 * full replication (Rep) version in that it calculates global energy values
 * via communication.
 * 
 */
template<class MT, class PT>
void DD_Source_Builder<MT,PT>::calc_initial_census(SP_Mesh mesh, 
						   SP_Mat_State state, 
						   SP_Opacity opacity,
						   SP_Rnd_Control rcontrol)
{
    Require(!census);

    // make the census on this processor
    census = new Particle_Buffer<PT>::Census();

    // calculate local source energies for volume emission and surface
    // source.  In the DD topology, these are Data_Distributed quantities.
    // calc_source_energies calculates local values of evoltot and esstot
    calc_source_energies(*state, *opacity);

    // calculate global source energies (required to estimate the number of
    // initial census particles)
    global_evoltot = evoltot;
    C4::gsum(global_evoltot);
    global_esstot = esstot;
    C4::gsum(global_esstot);

    // weak checks on local and global energies
    Check (evoltot <= global_evoltot);
    Check (esstot <= global_esstot);

    // make sure each processor has the same global energy values
    Check (parallel_data_op.check_global_equiv(global_evoltot));
    Check (parallel_data_op.check_global_equiv(global_esstot));

    // calculate local initial census energy
    calc_initial_ecen();

    // calculate global initial census energy
    global_ecentot = ecentot;
    C4::gsum(global_ecentot);
    Check (parallel_data_op.check_global_equiv(global_ecentot));

    // make a local, local-mesh sized field for the initial census random
    // number stream IDs
    ccsf_int cenrn(mesh);

    // calculate local (Data_Distributed) number of census particles per cell
    calc_initial_ncen(cenrn);

    // write out the initial census (local, same for all topologies)
    if (local_ncentot > 0)
	write_initial_census(mesh, rcontrol, local_ncen, local_ncentot,
			     cenrn); 
}

//---------------------------------------------------------------------------//
// CALCULATE THE NUMBER OF INITIAL CENSUS PARTICLES
//---------------------------------------------------------------------------//
/*!
 * \brief Calculate the number of initial census particles for the full
 * domain decomposition topology.
 *
 * This function calculates the number of initial census particles on each
 * processor.  Local number of census particles are calculated based on the
 * global totals of energy and numbers of census particles.  
 *
 * \param cenrn temporary, mutable, local field of random number stream IDs
 * for the initial census particles. 
 *
 */
template<class MT, class PT>
void DD_Source_Builder<MT,PT>::calc_initial_ncen(ccsf_int &cenrn)
{
    // The problem must have a nonzero source
    Insist ((global_evoltot + global_esstot + global_ecentot) != 0, 
	    "You must specify some source!");

    // No random numbers should have been used at this point
    Require (rtt_rng::rn_stream == 0);

    // Check consistency on vector lengths on this processor
    Require (cenrn.size() == local_ncen.size());

    // estimate the desired number of census particles from relative global
    // values of energy for each source species
    int ncenwant = static_cast<int>(global_ecentot / 
				    (global_evoltot + global_esstot + 
				     global_ecentot) * npwant); 

    // particles per unit energy
    double part_per_e;

    // attempting to make all census particles have the same energy weight, 
    // iterate on number of initial census particles
    bool retry = true;
    int  ntry = 0;
    int  ncenguess = ncenwant;

    // calculate global number of census particles per cell
    while (retry)
    {
	// update the attempt counter
	ntry++;
	
	// initialize calculated local number of intitial census particles 
	local_ncentot = 0;

	// limit guess to nonzero and nonnegative number
	if (ncenguess < 1)
	    ncenguess = 1;

	// update part_per_e
	if (global_ecentot > 0)
	    part_per_e = ncenguess / global_ecentot;
	else
	    part_per_e = 0.0;

	// calculate local number of census particles based on estimate of
	// (global) part_per_e.  
	calc_num_src_particles(part_per_e, ecen, local_ncen, local_ncentot); 

	// calculate global total number of intital census particles
	global_ncentot = local_ncentot;
	C4::gsum(global_ncentot);
    
	// check to see we haven't exceeded total particles for this cycle
	if (global_ncentot > ncenwant  &&  ntry < 10  &&  ncenguess > 1)
	    ncenguess -= (global_ncentot - ncenwant);
	else
	    retry = false;
    }

    // temporary variable for old ecen per cell
    double old_ecen = 0.0;

    // calculate local energy weights and sampling contribution to the energy 
    // loss.  NOTE that ecentot is not modified to be the true total; the
    // true total is ecentot - global_eloss_cen.
    for (int cell = 1; cell <= local_ncen.size(); cell++)
    {
	// calculate energy-weights of initial census particles in each cell 
	if (local_ncen(cell) > 0)
	    ew_cen(cell) = ecen(cell) / local_ncen(cell);
	else
	    ew_cen(cell) = 0;

	// save old ecen
	old_ecen = ecen(cell);

	// calculate new ecen which accounts for the energy loss from the
	// initial sampling.
	ecen(cell) = local_ncen(cell) * ew_cen(cell);

	// add up this processor's energy loss due to initial sampling. 
	local_eloss_cen += old_ecen - local_ncen(cell) * ew_cen(cell);
    }

    // sum up global sampling loss due to sampling.  The way we do this
    // calculation is possible because no cell is ever replicated in full
    // domain decomposition.
    global_eloss_cen = local_eloss_cen;
    C4::gsum(global_eloss_cen);

    // calculate the local number of particles and local random number stream 
    // IDs for the initial census
    calc_fullDD_rn_fields(topology->num_cells(), local_ncen, rtt_rng::rn_stream,
			  cenrn, global_ncentot);

    // Check global consistency of the random number stream state, noting
    // that rn_stream always starts the problem at zero
    Check  (parallel_data_op.check_global_equiv(rtt_rng::rn_stream));
    Ensure (rtt_rng::rn_stream == global_ncentot);
}

//---------------------------------------------------------------------------//
// RECALCULATE ENERGY-WEIGHTS OF CENSUS PARTICLES
//---------------------------------------------------------------------------//
/*!  
 * \brief Recalculate the census particles' energy-weights on each
 * processor.
 *
 * Because of its reproducible nature, the census comb does not conserve
 * energy exactly.  Therefore, we use the actual, pre-combed census energy
 * and the number of post-combed census particles in each cell to recalculate
 * the energy-weight in each cell.  Thus, energy is conserved at the cost of
 * some degradation in variance reduction.  In a full domain decomposition
 * topology, the census energies and number of particles in each cell ARE the
 * global values since no cells are replicated across processors.  (These two
 * quantities are "Data_Distributed.")
 * 
 */
template<class MT, class PT>
void DD_Source_Builder<MT,PT>::recalc_census_ew_after_comb(SP_Mesh mesh,
							   ccsf_int
							    &max_dead_rand_id, 
							   SP_Census
							    dead_census,
							   double 
							    &global_eloss_comb)
{
    // check sizes
    Require (local_ncen.size() == mesh->num_cells());

    // initialize checks
    int check_local_ncentot = 0;

    // <<<< resurrect dead census particles, if necessary  >>>>

    // initialize number and energy of resurrected particles
    int num_resurrected  = 0;
    double e_resurrected = 0.0;

    // make a temporary cell-centered scalar field to hold the post-comb, 
    // pre-resurrected values of local_ncen
    ccsf_int old_local_ncen(mesh);
    for (int cell = 1; cell <= mesh->num_cells(); cell++)
	old_local_ncen(cell) = local_ncen(cell);

    // loop through dead census particles and resurrect if there exists
    // globally unsampled energy and if the particle has the maximum random
    // number stream id in its cell
    while(dead_census->size() > 0)
    {
	// get the dead particle
	SP<PT> dead_particle = dead_census->top();
	dead_census->pop();

	// get the dead particle's attributes
	int loc_dead_cell = topology->local_cell(dead_particle->get_cell()); 
	Sprng dead_random = dead_particle->get_random();

	// resurrect the particle if needed
        if (old_local_ncen(loc_dead_cell) == 0)
	    if (ecen(loc_dead_cell) > 0.0)
		if (dead_random.get_num() == max_dead_rand_id(loc_dead_cell))
		{
		    e_resurrected += dead_particle->get_ew();
		    census->push(dead_particle);
		    local_ncen(loc_dead_cell)++;
		    local_ncentot++;
		    num_resurrected++;
		}
    }
    // obtain global number of resurrected particles
    C4::gsum(num_resurrected);

    // obtain global value of resurrected energy and update global census
    // energy loss due to combing
    if (num_resurrected > 0)
    {
	C4::gsum(e_resurrected);
	global_eloss_comb -= e_resurrected;
    }

    // <<<< set the post-comb census energy-weight in each cell >>>>
    for (int cell = 1; cell <= mesh->num_cells(); cell++)
    {
	if (local_ncen(cell) > 0)
	    ew_cen(cell) = ecen(cell) / local_ncen(cell);
	else
	    ew_cen(cell) = 0.0;

	check_local_ncentot += local_ncen(cell);
    }
    
    // have we conserved the total number of census particles on this proc?
    Ensure (check_local_ncentot == local_ncentot);
}


//---------------------------------------------------------------------------//
// CALCULATE NUMBERS OF SOURCE PARTICLES 
//---------------------------------------------------------------------------//
/*!
 * \brief Calculate the number of source particles for all sources.
 *
 * This function uses Source_Builder::calc_num_src_particles to get local
 * values (Data_Distributed) of the source numbers on processor.  It then uses
 * DD_Source_Builder::calc_fullDD_rn_fields to calculate local values
 * of the random number stream IDs per species per cell.
 *
 * Census numbers are calculated in preparation for combing.  This occurs
 * even if the census has just been created.  The census is combed \b every
 * cycle.
 */
template<class MT, class PT>
void DD_Source_Builder<MT,PT>::calc_source_numbers()
{
    Check(parallel_data_op.check_global_equiv(rtt_rng::rn_stream));

    // calculate total, global source energy 
    double total_energy = global_evoltot + global_ecentot + global_esstot;
    Require(total_energy > 0.0);

    // iterate over the problem energies to get source numbers

    // controllers for iteration
    int np_try_for = npwant;
    bool retry     = true;
    int num_try    = 0;

    // desired particles per unit energy
    double part_per_e;

    // global, total number of source particles
    int numtot;

    // iteration
    while (retry)
    {
	// increment counter
	num_try++;
	
	// initialize global, total number of source particles
	numtot = 0;

	// check on desired number of source particles
	if (np_try_for < 1)
	    np_try_for = 1;
	
	// calculate particles per energy for this iteration
	part_per_e = static_cast<double>(np_try_for) / total_energy;

	// estimate the number of local source particles for each source
	// species --> calculations are done on each processor without
	// communication 
	calc_num_src_particles(part_per_e, evol, local_nvol,
			       local_nvoltot); 
	
	calc_num_src_particles(part_per_e, ess, local_nss, local_nsstot); 
	
	calc_num_src_particles(part_per_e, ecen, local_ncen,
			       local_ncentot);

	// sum up source particles on this processor (used only for iteration 
	// stopping criterion)
	numtot = local_nvoltot + local_nsstot + local_ncentot;

	// calculate global number of source particles
	C4::gsum(numtot);
	
	// ending condition
	if (numtot > npwant && num_try < 10 && np_try_for > 1)
	    np_try_for -= (numtot - npwant);
	else
	    retry = false;
    }
    
    // now, calculate global number totals per species
    global_nvoltot = local_nvoltot;
    C4::gsum(global_nvoltot);
    global_nsstot = local_nsstot;
    C4::gsum(global_nsstot);
    global_ncentot = local_ncentot;
    C4::gsum(global_ncentot);

    // check global equivalence of source numbers
    Check(parallel_data_op.check_global_equiv(global_nvoltot));
    Check(parallel_data_op.check_global_equiv(global_nsstot));
    Check(parallel_data_op.check_global_equiv(global_ncentot));
    
    // calculate local energy weights and global sampling loss 
    for (int cell = 1; cell <= ew_cen.size(); cell++)
    {
	// calculate energy weights of volume source particles per cell
	if (local_nvol(cell) > 0) 
	    ew_vol(cell) = evol(cell) /
		static_cast<double>(local_nvol(cell));
	else 
	    ew_vol(cell) = 0;
	
	// calculate energy weights of surface source particles per cell
	if (local_nss(cell) > 0) 
	    ew_ss(cell) = ess(cell) /
		static_cast<double>(local_nss(cell));
	else 
	    ew_ss(cell) = 0;

	// calculate energy weights of census particles per cell
	if (local_ncen(cell) > 0) 
	    ew_cen(cell) = ecen(cell) /
		static_cast<double>(local_ncen(cell));
	else 
	    ew_cen(cell) = 0;	
	
	// add up local sampling contribution to energy loss for volume
	// emission, surface source, and census
	local_eloss_vol += evol(cell) - local_nvol(cell) * ew_vol(cell);
	local_eloss_ss  += ess(cell)  - local_nss(cell)  * ew_ss(cell);
	local_eloss_cen += ecen(cell) - local_ncen(cell) * ew_cen(cell);
    }

    // calculate global energy losses for this cycle.  some global census
    // energy loss may exist from the initial census calculation if this is
    // the first IMC cycle.
    global_eloss_vol = local_eloss_vol;
    C4::gsum(global_eloss_vol);
    global_eloss_ss = local_eloss_ss;
    C4::gsum(global_eloss_ss);
    double new_eloss_cen = local_eloss_cen;
    C4::gsum(new_eloss_cen);
    global_eloss_cen += new_eloss_cen;


    // calculate local source numbers from the global values for each source
    // species
    int global_numcells = topology->num_cells();

    calc_fullDD_rn_fields(global_numcells, local_nvol, rtt_rng::rn_stream,
			  volrn, global_nvoltot);

    calc_fullDD_rn_fields(global_numcells, local_nss, rtt_rng::rn_stream,
			  ssrn, global_nsstot);

    // every processor have the same updated random number stream ID?
    Check(parallel_data_op.check_global_equiv(rtt_rng::rn_stream));
}

//---------------------------------------------------------------------------//
// CALCULATE RANDOM NUMBER IDS AND NUMBER OF SOURCE PARTICLES
//---------------------------------------------------------------------------//
/*!  
 * \brief Calculate the local starting random number stream IDs in each cell
 * for each species of source.
 *
 * calc_fullDD_rn_fields calculates the starting random number stream ID in
 * each local cell.  This function is general in that it works for each
 * species of source particle (initial census, volume emission, surface
 * source).  We assume that the local number of source particles per species
 * has already been calculated.  This function calculates the random number
 * stream ID offset based on the numbers of particles in the globally indexed
 * cells, and then it calculates the local starting random number stream IDs.
 *
 * The underlying reason for this function is to make the source calculation
 * reproducible regardless of the number of processors.
 *
 * \param global_numcells constant number of cells in the global mesh
 *
 * \param local_n_field constant numbers of source particles (for a given species)
 * on the local mesh.
 *
 * \param next_avail_rn mutable value of next available global random number
 * stream ID.
 *
 * \param rn_field mutable field of starting random number stream IDs
 *
 * \param global_ntot constant global total of source particles for the given 
 * species -- used as a check.
 *
 */
template<class MT, class PT> 
void DD_Source_Builder<MT,PT>::calc_fullDD_rn_fields(const int
						     global_numcells, 
						     ccsf_int &local_n_field,
						     int &next_avail_rn, 
						     ccsf_int &rn_field, 
						     const int global_ntot)
{
    // require that source number and random number stream IDs are equally
    // sized.  require (weak check) that local and global cell numbers are OK
    int num_cells = local_n_field.size();
    Require (rn_field.size() == topology->num_cells(C4::node()));
    Require (rn_field.size() == num_cells);
    Require ((C4::nodes() > 1) ? (num_cells < global_numcells) 
	     : (num_cells == global_numcells));

    // make a temporary vector for global source numbers
    sf_int global_n_field(global_numcells, 0);

    // map local source numbers into temporary global vector
    parallel_data_op.local_to_global(local_n_field, global_n_field,
				     Parallel_Data_Operator::Data_Distributed());

    // check the sum of the global source numbers vector
    Check (global_ntot == accumulate(global_n_field.begin(),
				     global_n_field.end(), 0));
    
    // temporary vector holding the global random number stream ID offsets.
    // Note that the extra step and storage from having the global_offset
    // vector absolves us from making assumptions about the monotonicity of
    // the local-to-global cell index mapping. The two steps and the extra
    // storage allow us to handle, e.g., global_cell(1)=45, global_cell(2)=3.
    sf_int global_offset(global_numcells, 0);
    int rn_offset = 0;
    for (int gc = 1; gc <= global_numcells; gc++)
    {
	global_offset[gc-1] = rn_offset;
	rn_offset          += global_n_field[gc-1];
    }
    
    // map local random number stream ID offset from global offsets
    for (int loc = 1; loc <= num_cells; loc++)
	rn_field(loc) = global_offset[topology->global_cell(loc)-1] +
	    next_avail_rn;  
    
    // check that final stream ID offset (not explicitly used) matches total
    // number of source particles.
    Check (rn_offset == global_ntot);
    
    // update the next available random number stream ID
    next_avail_rn += global_ntot;   
}


} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                        end of imc/DD_Source_Builder.t.hh
//---------------------------------------------------------------------------//
