//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Rep_Source_Builder.t.hh
 * \author Todd Urbatsch and Thomas M. Evans
 * \date   Thu Dec  9 10:31:16 1999
 * \brief  Implementation file for Rep_Source_Builder.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Rep_Source_Builder.hh"
#include "Mesh_Operations.hh"
#include "mc/Comm_Patterns.hh"
#include "rng/Random.hh"

namespace rtt_imc
{

//---------------------------------------------------------------------------//
// SOURCE BUILDER MAIN INTERFACE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Build the source in a full replication (Domain Replication)
 * parallel topology.
 *
 * \param mesh rtt_dsxx::SP to a IMC MT
 * \param state rtt_dsxx::SP to a rtt_imc::Mat_State
 * \param opacity rtt_dsxx::SP to an rtt_imc::Opacity
 * \param rnd_control rtt_dsxx::SP to a rtt_rng::Rnd_Control
 * \param patterns null rtt_dsxx::SP for domain replication
 */
template<class MT, class FT, class PT>
typename Rep_Source_Builder<MT,FT,PT>::SP_Source
Rep_Source_Builder<MT,FT,PT>::build_Source(SP_Mesh          mesh,
					   SP_Mat_State     state,
					   SP_Opacity       opacity,
					   SP_Rnd_Control   rnd_control,
					   SP_Comm_Patterns patterns)
{
    int num_cells = mesh->num_cells();
    Require(num_cells == state->num_cells());
    Require(num_cells == opacity->num_cells());

    // if we don't have a census object, build it
    if (!Base::census)
	calc_initial_census(mesh, state, opacity, rnd_control);
    else
    {
	// calculate the volume emission and surface source energies
	Base::calc_source_energies(*state, *opacity);

	// make sure global energies are the same on each processor
	Check(Base::parallel_data_op.check_global_equiv(Base::evoltot));
	Check(Base::parallel_data_op.check_global_equiv(Base::esstot));
	Check(Base::parallel_data_op.check_global_equiv(Base::ecentot));
    }

    // calculate the number of source particles for all source particles
    calc_source_numbers();

    // comb census -- called even if census is empty; must update local_ncen
    double    eloss_comb = 0;
    SP_Census dead_census(new Census());
    ccsf_int  max_dead_rand_id(mesh);
    Base::comb_census(rnd_control, local_ncen, local_ncentot, eloss_comb,
		      dead_census, max_dead_rand_id);

    // calculate global values of energy loss and numbers of particles due
    // to combing the census
    double global_eloss_comb = eloss_comb;
    C4::gsum(global_eloss_comb);
    
    global_ncentot = local_ncentot;
    C4::gsum(global_ncentot);

   // recalculate and reset the census particles' energy-weights
    recalc_census_ew_after_comb(mesh, max_dead_rand_id, dead_census,
				global_eloss_comb); 

    // reset the census particles' energy-weights and update (presumably
    // improve) the global census energy loss due to combing. actual_ecentot
    // is the actual, total, global census energy represented by the
    // particles; it doesn't contain the energy that was lost during
    // sampling; ecentot is what the total energy is supposed to be if there
    // was no sampling loss.
    double actual_ecentot = Base::ecentot - global_eloss_cen;
    Base::reset_ew_in_census(local_ncentot, global_eloss_comb, actual_ecentot); 

    // add energy loss from the comb to the global census energy loss
    global_eloss_cen += global_eloss_comb;

    // build Mesh_Operations class for source
    rtt_dsxx::SP<Mesh_Operations<MT> > mesh_op(
	new Mesh_Operations<MT>(mesh, state, Base::topology, patterns)); 

    // build the source --> the source gets a SP to the census on this
    // processor, this means that as particles are taken out of the source
    // the "master" census (in a global_state object) on processor is emptied 
    // as well
    SP_Source source(new Source<MT,FT,PT>(
			 Base::volrn, local_nvol, Base::ew_vol, Base::ssrn,
			 local_nss, Base::ss_face_in_cell, Base::ew_ss, 
			 Base::census, Base::ss_dist, local_nvoltot,
			 local_nsstot, rnd_control, state, mesh_op, 
			 Base::topology, opacity, Base::freq_samp_data));

    // return the source
    return source;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Build a local initial census from radiation temperature.
 *
 * This function is a pure virtual function in the Source_Builder base class.
 * We make it part of the public interface to give the option of building an
 * initial census and then exiting the source building process.  This is
 * useful in rad-hydro applications where the initial census must be built on
 * cycle 0, before a hydro step, from the initial radiation temperature.  The
 * build_Source() function checks to see if a census exists and then builds
 * it if it does not.  The census is loaded into the Source_Builder base
 * class through the run-time interface.
 *
 * This function builds a local initial census on each processor.  The
 * rtt_dsxx::SP to the completed census may be accessed through the
 * Source_Builder::get_census() function. 
 */
template<class MT, class FT, class PT>
void Rep_Source_Builder<MT,FT,PT>::calc_initial_census(SP_Mesh        mesh,
						       SP_Mat_State   state, 
						       SP_Opacity     opacity,
						       SP_Rnd_Control rcontrol)
{
    using rtt_dsxx::SP;

    Require(!Base::census);

    // make the census
    Base::census = new Census();

    // calculate global source energies for volume emission and surface
    // source -> Data_Replicated quantities -> calc_source_energies
    // calculates local values of evoltot and esstot, which for full
    // replication, are equivalent to the global values
    Base::calc_source_energies(*state, *opacity);
    Check(Base::parallel_data_op.check_global_equiv(Base::evoltot));
    Check(Base::parallel_data_op.check_global_equiv(Base::esstot));
    
    // calculate local (Data_Replicated) values of initial census energy
    Base::calc_initial_ecen(*opacity);
    Check(Base::parallel_data_op.check_global_equiv(Base::ecentot));

    // make a local, global-mesh sized, field holding the census random
    // number IDs on each processor
    ccsf_int cenrn(mesh);

    // calculate local (Data_Replicated) number of census particles per cell
    calc_initial_ncen(cenrn);

    // get the frequency
    SP<FT> frequency = opacity->get_Frequency();

    // write out the initial census
    if (local_ncentot > 0)
	Base::write_initial_census(mesh, rcontrol, *frequency, *state,
				   local_ncen, local_ncentot, cenrn);   
}

//===========================================================================//
// REP_SOURCE_BUILDER IMPLEMENTATION FUNCTIONS
// Private functions in the derived class that are used to build sources.
// These functions are unique to the derived class and usually deal with the
// way that data is parsed across processors.
//===========================================================================//

//---------------------------------------------------------------------------//
// CALCULATE INITIAL NUMBER OF CENSUS PARTICLES
//---------------------------------------------------------------------------//
/*!  
 * \brief Calculate the initial number of census particles on each
 * processor.
 *
 * This function calculates the number of initial census particles.  Both
 * local numbers and the global number of census particles are calculated.
 *
 * \param cenrn temporary, mutable, local field of initial census random
 * number IDs 
 */
template<class MT, class FT, class PT>
void Rep_Source_Builder<MT,FT,PT>::calc_initial_ncen(ccsf_int &cenrn)
{
    // require that the global source energies are not zero
    Insist ((Base::evoltot + Base::esstot + Base::ecentot) != 0,
	    "You must specify some source!");
    Require(rtt_rng::rn_stream == 0);
    Require(cenrn.size() == global_ncen.size());

    // estimate the number of census particles that are desired from global
    // values of the total energy of each source
    int ncenwant = static_cast<int>(
	Base::ecentot / (Base::evoltot + Base::esstot + Base::ecentot) * 
	Base::npwant);

    // particles per unit energy
    double part_per_e = 0.0;

    // attempt to make all census particles have the same energy weight,
    // iterate on number of initial census particles
    bool retry     = true;
    int  ntry      = 0;
    int  ncenguess = ncenwant;

    // calculate global number of census particles per cell
    while (retry)
    {
	// update the attempt counter
	ntry++;
	
	// initialize calculated quantities
	global_ncentot = 0;

	// iteration parameters
	if (ncenguess < 1)
	    ncenguess = 1;

	// update part_per_e
	if (Base::ecentot > 0)
	    part_per_e = ncenguess / Base::ecentot;
	else
	    part_per_e = 0.0;

	// calculate global number of census particles (on each processor)
	// based on estimate of part_per_e -> Note: this requires no
	// processor-processor communication
	Base::calc_num_src_particles(part_per_e, Base::ecen, global_ncen,
				     global_ncentot);
    
	// check to see we haven't exceeded total particles for this cycle
	if (global_ncentot > ncenwant  &&  ntry < 10  &&  ncenguess > 1)
	    ncenguess -= (global_ncentot - ncenwant);
	else
	    retry = false;
    }

    // calculate global energy weights and sampling contribution to energy
    // loss.  NOTE that ecentot is not modified to be the actual total
    // represented by the census particles; ecentot is what the total is
    // supposed to be.  The actual total is ecentot - global_eloss_cen =
    // sum(ncen*ewcen).
    for (int cell = 1; cell <= global_ncen.size(); cell++)
    {
	if (global_ncen(cell) > 0)
	    ew_cen(cell) = ecen(cell) / global_ncen(cell);
	else
	{
	    // add up energy loss due to unsampled census energy
	    global_eloss_cen += ecen(cell);

	    // zero out ew_cen and ecen since cell is unsampled or has no
	    // energy.
	    ew_cen(cell) = 0.0;
	    ecen(cell)   = 0.0;
	}
    }

    // calculate the local number of particles and local random number IDs
    calc_num_part_and_rn_fields(global_ncen, global_ncentot,
				rtt_rng::rn_stream, local_ncen,
				local_ncentot, cenrn);

    // make checks, noting that rn_stream always starts the problem at zero
    Check  (Base::parallel_data_op.check_global_equiv(rtt_rng::rn_stream));
    Ensure (rtt_rng::rn_stream == global_ncentot);
}

//---------------------------------------------------------------------------//
// RECALCULATE ENERGY-WEIGHTS OF CENSUS PARTICLES
//---------------------------------------------------------------------------//
/*!  
 * \brief Recalculate the census particles' energy-weights on each
 * processor using global information.
 *
 * First, if global census energy is unsampled because of the comb, this
 * function resurrects one census particle and modifies the global census
 * energy loss accordingly.  The particle to be resurrected is the one with
 * the maximum random number stream id.  We chose this criterion because it
 * is reproducible and because any bias would presumably be negligible.  If,
 * due to spawning, more than one particle has the maximum stream id, all
 * particles with the maximum stream id get resurrected.  The resurrection
 * gives up some reduction in variance in order to conserve energy globally
 * and, in a spatial sense, locally. 
 *  
 * Because of its reproducible nature, the census comb does not conserve
 * energy exactly.  Therefore, we use the actual, global, pre-combed census
 * energy and the global number of post-combed (and resurrected) census
 * particles in each cell to recalculate the energy-weight in each cell.  In
 * a fully replicated topology, the census energies and number of particles
 * in each cell must be summed over all processors.  (These two quantities
 * are "Data_Decomposed.")
 *
 * \param mesh smart pointer to the mesh.
 * \param max_dead_rand_id mutable mesh-sized array containing local maximum
 *                         values of the random number stream id's of census 
 *                         particles that were combed out.
 * \param dead_census smart pointer to the dead, combed out census particles.
 * \param global_eloss_comb mutable global value of energy loss due to
 *                          combing; this loss will be decreased by the
 *                          energy of any resurrected census particles.
 *  */
template<class MT, class FT, class PT>
void Rep_Source_Builder<MT,FT,PT>::recalc_census_ew_after_comb(
    SP_Mesh    mesh,
    ccsf_int  &max_dead_rand_id, 
    SP_Census  dead_census,
    double    &global_eloss_comb)
{
    using rtt_mc::Parallel_Data_Operator;
    using rtt_rng::Sprng;
    using rtt_dsxx::SP;

    // check sizes
    Require (global_ncen.size() == mesh->num_cells());
    Require (local_ncen.size()  == mesh->num_cells());

    // initialize checks
    int check_global_ncentot = 0;

    // <<<< map local census numbers to global census numbers >>>>
    for (int cell = 1; cell <= mesh->num_cells(); cell++)
	global_ncen(cell) = local_ncen(cell);
    Base::parallel_data_op.local_to_global(
	local_ncen, global_ncen, Parallel_Data_Operator::Data_Decomposed());

    // <<<< resurrect dead census particles, if necessary  >>>>

    // initialize number and energy of resurrected particles
    int num_resurrected  = 0;
    double e_resurrected = 0.0;

    // map, in full replication topology, the global max random number stream 
    // id to local array
    int *local_max = new int[mesh->num_cells()];
    for (int cell = 1; cell <= mesh->num_cells(); cell++)
	local_max[cell-1] = max_dead_rand_id(cell);
    C4::gmax(local_max, mesh->num_cells());
    for (int cell = 1; cell <= mesh->num_cells(); cell++)
	max_dead_rand_id(cell) = local_max[cell-1];
    delete [] local_max;

    // loop through dead census particles and resurrect if there exists
    // globally unsampled energy and if the particle has the maximum random
    // number stream id in its cell
    while(dead_census->size() > 0)
    {
	// get the dead particle and its attributes
	SP<PT> dead_particle = dead_census->top();
	dead_census->pop();
	int dead_cell        = dead_particle->get_cell();
	Sprng dead_random    = dead_particle->get_random();

	// resurrect the particle if needed
        if (global_ncen(dead_cell) == 0)
	    if (ecen(dead_cell) > 0.0)
		if (dead_random.get_num() == max_dead_rand_id(dead_cell))
		{
		    e_resurrected += dead_particle->get_ew();
		    Base::census->push(dead_particle);
		    local_ncen(dead_cell)++;
		    local_ncentot++;
		    num_resurrected++;
		}
    }
    // obtain global number of resurrected particles
    C4::gsum(num_resurrected);

    // <<<< recalculate global census numbers if any census particle were
    // resurrected >>>>
    if (num_resurrected > 0)
    {
	// map new local census numbers to global census numbers 
	for (int cell = 1; cell <= mesh->num_cells(); cell++)
	    global_ncen(cell) = local_ncen(cell);

	Base::parallel_data_op.local_to_global(
	    local_ncen, global_ncen, Parallel_Data_Operator::Data_Decomposed());

	// update the total
	global_ncentot += num_resurrected;

	// obtain global value of resurrected energy and update global census 
	// energy loss due to combing
	C4::gsum(e_resurrected);
	global_eloss_comb -= e_resurrected;
    }


    // <<<< set the post-comb census energy-weight in each cell <<<<
    for (int cell = 1; cell <= mesh->num_cells(); cell++)
    {
	if (global_ncen(cell) > 0)
	    ew_cen(cell) = ecen(cell) / global_ncen(cell);
	else
	    ew_cen(cell) = 0.0;

	check_global_ncentot += global_ncen(cell);
    }
    
    // make sure we have considered the proper total number of census
    // particles during our recalculation of the ew's
    Ensure (check_global_ncentot == global_ncentot);
}

//---------------------------------------------------------------------------//
// CALCULATE SOURCE NUMBERS IN A CYCLE
//---------------------------------------------------------------------------//
/*!
 * \brief Calculate the number of source particles for all sources.
 *
 * This function uses Source_Builder::calc_num_src_particles to get global
 * values (Data_Replicated) of the source numbers on processor.  It then uses
 * Rep_Source_Builder::calc_num_part_and_rn_fields to calculate local values
 * of the number of source particles per species per cell.
 *
 * Census numbers are calculated in preparation for combing.  This occurs
 * even if the census has just been created.  The census is combed \b every
 * cycle.
 */
template<class MT, class FT, class PT>
void Rep_Source_Builder<MT,FT,PT>::calc_source_numbers()
{
    Check(Base::parallel_data_op.check_global_equiv(rtt_rng::rn_stream));

    // calculate total, global source energy 
    double total_energy = Base::evoltot + Base::ecentot + Base::esstot;
    Require(total_energy > 0.0);

    // iterate over the problem energies to get source numbers

    // controllers for iteration
    int np_try_for = Base::npwant;
    bool retry     = true;
    int num_try    = 0;

    // desired particles per unit energy
    double part_per_e = 0.0;

    // global, total number of source particles
    int numtot = 0;

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

	// estimate the number of global source particles for each source
	// species --> calculations are done on each processor without
	// communication 
	Base::calc_num_src_particles(part_per_e, Base::evol, global_nvol,
				     global_nvoltot); 
	
	Base::calc_num_src_particles(part_per_e, Base::ess, global_nss,
				     global_nsstot); 
	
	Base::calc_num_src_particles(part_per_e, Base::ecen, global_ncen,
				     global_ncentot);

	// sum up source particles
	numtot = global_nvoltot + global_nsstot + global_ncentot;
	
	// ending condition
	if (numtot > Base::npwant && num_try < 10 && np_try_for > 1)
	    np_try_for -= (numtot - Base::npwant);
	else
	    retry = false;
    }
    
    // check global equivalence of source numbers
    Check(Base::parallel_data_op.check_global_equiv(global_nvoltot));
    Check(Base::parallel_data_op.check_global_equiv(global_nsstot));
    Check(Base::parallel_data_op.check_global_equiv(global_ncentot));
    
    // calculate energy weights sampling loss constribution to energy loss
    // based on global source numbers and energies
    for (int cell = 1; cell <= Base::ew_cen.size(); cell++)
    {
	// calculate energy weights of volume source particles per cell
	if (global_nvol(cell) > 0) 
	    ew_vol(cell) = evol(cell) /
		static_cast<double>(global_nvol(cell));
	else 
	    ew_vol(cell) = 0;
	
	// calculate energy weights of surface source particles per cell
	if (global_nss(cell) > 0) 
	    ew_ss(cell) = ess(cell) / static_cast<double>(global_nss(cell));
	else 
	    ew_ss(cell) = 0;

	// calculate energy weights of census particles per cell
	if (global_ncen(cell) > 0) 
	    ew_cen(cell) = ecen(cell) / static_cast<double>(global_ncen(cell));
	else 
	    ew_cen(cell) = 0;	
	
	// add up sampling contribution to energy loss for volume, ss, and
	// census.  Some global census energy loss may exist from the initial 
	// census calculation if this is the first IMC cycle.
	global_eloss_vol += evol(cell) - (global_nvol(cell) * ew_vol(cell));
	global_eloss_ss  += ess(cell)  - (global_nss(cell)  * ew_ss(cell));
	global_eloss_cen += ecen(cell) - (global_ncen(cell) * ew_cen(cell));
    }

    // calculate local source numbers from the global values for each source
    // species
    calc_num_part_and_rn_fields(global_nvol, global_nvoltot,
				rtt_rng::rn_stream, local_nvol,
				local_nvoltot, Base::volrn);
    
    calc_num_part_and_rn_fields(global_nss, global_nsstot,
				rtt_rng::rn_stream, local_nss,
				local_nsstot, Base::ssrn);

    Check(Base::parallel_data_op.check_global_equiv(rtt_rng::rn_stream));
}

//---------------------------------------------------------------------------//
// CALCULATE RANDOM NUMBER IDS AND NUMBER OF SOURCE PARTICLES
//---------------------------------------------------------------------------//
/*!  
 * \brief Calculate the local random number stream IDs and the local number
 * of particles per cell per species on the local processor.
 *
 * This function calculates the local number of particles per cell per
 * species on each processor.  We assume that the global number of particles
 * per cell per species has been calculated previously. Additionally, the
 * starting random number stream IDs are calculated for each cell for the
 * given source species.  Data calculated in this function is written
 * directly into mutable function arguments described below.
 *
 * \param global_n_field constant field of global number of source particles
 *                       for a particular species
 * \param global_numtot  constant global total number of source particles for
 *                       a particular species
 * \param next_avail_rn  mutable value of next available global random number
 *                       stream ID
 * \param local_n_field  mutable field of local number of source particles for
 *                       a particular species on processor
 * \param local_numtot   mutable value of the local total number of source
 *                       particles for a particular species on processor
 * \param rn_field       mutable field of starting random number stream IDs
 */
template<class MT, class FT, class PT> 
void Rep_Source_Builder<MT,FT,PT>::calc_num_part_and_rn_fields(
    const ccsf_int &global_n_field,
    const int       global_numtot,
    int            &next_avail_rn,
    ccsf_int       &local_n_field,
    int            &local_numtot, 
    ccsf_int       &rn_field)
{
    // requirements
    int num_cells = global_n_field.size();
    Require(local_n_field.size() == num_cells);
    Require(rn_field.size() == num_cells);

    // set total number to zero
    local_numtot = 0;

    // initial starting random number stream id
    int initial_avail_rn = next_avail_rn;

    // number of global particles leftover after doing integer math on the
    // global number of particles per cell
    int leftover;

    // integer number of particles per cell
    int even_spread;

    // processor on which to place the first of the next leftover particles;
    // incremented by the value of "leftover" after each cell with leftovers.
    int first_proc_for_leftovers = 0;

    // sweep through cells and determine the number of source particles per
    // cell for this species
    for (int cell = 1; cell <= num_cells; cell++)
    {
	if (global_n_field(cell) > 0)
	{
	    // calculate the integer number of particles per cell
	    even_spread = global_n_field(cell) / C4::nodes();
	    
	    // calculate the remainder left from integer division
	    leftover = global_n_field(cell) - C4::nodes() * even_spread;
	    Check(leftover < C4::nodes());

	    // for this cell, shifted_node is a local renumbering of
	    // processors beginning with the first processor that gets any
	    // leftover particles
	    int shifted_node = (C4::node() + C4::nodes() - 
				first_proc_for_leftovers) % C4::nodes();  
	    Check (shifted_node >= 0);
	    Check (shifted_node < C4::nodes());
	    
	    // set the starting random number stream ID for this cell on this
	    // processor. Recall that, when processor nodes are numbered
	    // beginning with zero, there are node() processors "below"
	    // processor = node().
	    rn_field(cell) = next_avail_rn + shifted_node * even_spread +
		rtt_mc::global::min(leftover, shifted_node);

	    // set local field; from the remainder (leftover) put one on each
	    // of the first leftover processors
	    local_n_field(cell) = even_spread;
	    if (shifted_node < leftover)
		local_n_field(cell)++;	

	    // increment next_avail_rn number by global number of particles
	    // in this cell
	    next_avail_rn += global_n_field(cell);
	    
	    // update check of total number of source particles
	    local_numtot += local_n_field(cell);

	    // increment the ID of the next processor to get any leftovers
	    first_proc_for_leftovers = (first_proc_for_leftovers + leftover)
		% C4::nodes();
	    Check (first_proc_for_leftovers >= 0);
	    Check (first_proc_for_leftovers < C4::nodes());
	}
	else
	{
	    local_n_field(cell) = 0;
	    rn_field(cell)      = next_avail_rn;
	}
    }
    
    // do a check on random number id
    Ensure(next_avail_rn == initial_avail_rn + global_numtot);
}

} // end of rtt_imc

//---------------------------------------------------------------------------//
//                        end of imc/Rep_Source_Builder.t.hh
//---------------------------------------------------------------------------//
