//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Rep_Source_Builder.t.hh
 * \author Thomas M. Evans
 * \date   Thu Dec  9 10:31:16 1999
 * \brief  Implementation file for Rep_Source_Builder.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Rep_Source_Builder.hh"
#include "Mesh_Operations.hh"
#include "Global.hh"
#include "c4/global.hh"
#include <iomanip>

namespace rtt_imc
{

using C4::nodes;
using C4::node;
using rtt_dsxx::SP;

using std::cout;
using std::endl;
using std::ios;

//---------------------------------------------------------------------------//
// constructor
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for Rep_Source_Builder.
 */
template<class MT, class PT>
template<class IT>
Rep_Source_Builder<MT,PT>::Rep_Source_Builder(SP<IT> interface, SP_Mesh mesh, 
					      SP_Topology top)
    : Source_Builder<MT,PT>(interface, mesh, top),
      global_ncen(mesh),
      global_ncentot(0),
      local_ncen(mesh),
      local_ncentot(0),
      global_eloss_cen(0),
      global_nvol(mesh),
      global_nvoltot(0),
      local_nvol(mesh),
      local_nvoltot(0),
      global_nss(mesh),
      global_nsstot(0),
      local_nss(mesh),
      local_nsstot(0),
      global_eloss_vol(0),
      global_eloss_ss(0)
{ 
    Check(parallel_data_op.check_global_equiv(rtt_rng::rn_stream));

    // if the census does not exist yet then we build it --> if the census
    // does not exist it means that this is the initial IMC cycle
    if (census)
    {
	// sweep through cells and fill in persistent data

	// checks of persistent data
	double ecentot_check = 0.0;
	for (int cell = 1; cell <= mesh->num_cells(); cell++)
	{
	    // get global values of ecen from the interface
	    ecen(cell)     = interface->get_ecen(cell);
	    ecentot_check += ecen(cell);
	    
	    // more may follow --> probably cumulative edits
	}
	
	// get ecentot
	ecentot = interface->get_ecentot();
    
	// checks
	Check(rtt_mc::global::soft_equiv(ecentot, ecentot_check, 1.e-12));
    }
}

//---------------------------------------------------------------------------//
// SOURCE BUILDER MAIN INTERFACE FUNCTIONS
//---------------------------------------------------------------------------//

template<class MT, class PT>
Rep_Source_Builder<MT,PT>::SP_Source
Rep_Source_Builder<MT,PT>::build_Source(SP_Mesh mesh,
					SP_Mat_State state,
					SP_Opacity opacity,
					SP_Rnd_Control rnd_control)
{
    int num_cells = mesh->num_cells();
    Require(num_cells == state->num_cells());
    Require(num_cells == opacity->num_cells());

    // if we don't have a census build it
    if (!census)
	calc_initial_census(mesh, state, opacity, rnd_control);
    else
    {
	// calculate the volume emission and surface source energies
	calc_source_energies(*state, *opacity);

	// make sure global energies are the same on each processor
	Check(parallel_data_op.check_global_equiv(evoltot));
	Check(parallel_data_op.check_global_equiv(esstot));
	Check(parallel_data_op.check_global_equiv(ecentot));
    }

    // calculate the number of source particles for all source particles
    calc_source_numbers();

    // comb census
    double eloss_comb = 0;
    if (census->size() > 0)
	comb_census(rnd_control, local_ncentot, eloss_comb);

    // update the global values of eloss_cen and the post comb total numbers
    // of census particles
    C4::gsum(eloss_comb);
    global_eloss_cen += eloss_comb;
    
    global_ncentot = local_ncentot;
    C4::gsum(global_ncentot);

    // build Mesh_Operations class for source
    SP<Mesh_Operations<MT> > mesh_op
	(new Mesh_Operations<MT>(mesh, state, topology)); 

    // build the source --> the source gets a SP to the census on this
    // processor, this means that as particles are taken out of the source
    // the "master" census (in a global_state object) on processor is emptied 
    // as well
    SP_Source source(new Source<MT,PT>(volrn, local_nvol, ew_vol, ssrn,
				       local_nss, ss_face_in_cell, ew_ss, 
				       census, ss_dist, local_nvoltot,
				       local_nsstot, rnd_control, state,
				       mesh_op));

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
 * useful in rad-hydro applications where the initial census must be build on
 * cycle 0, before a hydro step, from the initial radiation temperature.  The
 * build_Source() function checks to see if a census exists and then builds
 * it if it does not.  The census is loaded into the Source_Builder base
 * class through the run-time interface.
 *
 * This function builds a local initial census on each processor.  The
 * rtt_dsxx::SP to the completed census may be accessed through the
 * Source_Builder::get_census() function. 
 */
template<class MT, class PT>
void Rep_Source_Builder<MT,PT>::calc_initial_census(SP_Mesh mesh,
						    SP_Mat_State state, 
						    SP_Opacity opacity,
						    SP_Rnd_Control rcontrol)
{
    Require(!census);

    // make the census
    census = new Particle_Buffer<PT>::Census();

    // calculate global source energies for volume emission and surface
    // source -> Data_Replicated quantities -> calc_source_energies
    // calculates local values of evoltot and esstot, which for full
    // replication, are equivalent to the global values
    calc_source_energies(*state, *opacity);
    Check(parallel_data_op.check_global_equiv(evoltot));
    Check(parallel_data_op.check_global_equiv(esstot));
    
    // calculate local (Data_Replicated) values of initial census energy
    calc_initial_ecen();
    Check(parallel_data_op.check_global_equiv(ecentot));

    // make a local, global-mesh sized, field holding the census random
    // number IDs on each processor
    ccsf_int cenrn(mesh);

    // calculate local (Data_Replicated) number of census particles per cell
    calc_initial_ncen(cenrn);

    // write out the initial census
    if (local_ncentot > 0)
	write_initial_census(mesh, rcontrol, local_ncen, local_ncentot,
			     cenrn);   
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
template<class MT, class PT>
void Rep_Source_Builder<MT,PT>::calc_initial_ncen(ccsf_int &cenrn)
{
    // require that the global source energies are not zero
    Insist ((evoltot+esstot+ecentot) != 0, "You must specify some source!");
    Require(rtt_rng::rn_stream == 0);
    Require(cenrn.size() == global_ncen.size());

    // estimate the number of census particles that are desired from global
    // values of the total energy of each source
    int ncenwant = static_cast<int>((ecentot) / (evoltot + esstot + ecentot) 
				    * npwant);

    // particles per unit energy
    double part_per_e;

    // attempt to make all census particles have the same energy weight,
    // iterate on number of initial census particles
    bool   retry = true;
    int    ntry = 0;
    int    ncenguess = ncenwant;

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
	if (ecentot > 0)
	    part_per_e = ncenguess / ecentot;
	else
	    part_per_e = 0.0;

	// calculate global number of census particles (on each processor)
	// based on estimate of part_per_e -> Note: this requires no
	// processor-processor communication
	calc_num_src_particles(part_per_e, ecen, global_ncen,
			       global_ncentot);
    
	// check to see we haven't exceeded total particles for this cycle
	if (global_ncentot > ncenwant  &&  ntry < 10  &&  ncenguess > 1)
	    ncenguess -= (global_ncentot - ncenwant);
	else
	    retry = false;
    }

    // temporary variable for the old global ecen per cell
    double old_ecen = 0.0;

    // calculate global energy weights and sampling contribution to
    // energy loss.  NOTE that ecentot is not modified to be the true
    // total; the true total is ecentot - global_eloss_cen.
    for (int cell = 1; cell <= global_ncen.size(); cell++)
    {
	if (global_ncen(cell) > 0)
	    ew_cen(cell) = ecen(cell) / global_ncen(cell);
	else
	    ew_cen(cell) = 0.0;

        // save old ecen
        old_ecen = ecen(cell);

        // calculate new global ecen, which won't contain energy lost
        // from initial sampling
        ecen(cell) = global_ncen(cell) * ew_cen(cell);

        // add up energy loss due to initial census sampling
        global_eloss_cen += old_ecen - global_ncen(cell) * ew_cen(cell);
    }

    // calculate the local number of particles and local random number IDs
    calc_num_part_and_rn_fields(global_ncen, global_ncentot,
				rtt_rng::rn_stream, local_ncen,
				local_ncentot, cenrn);

    Check(parallel_data_op.check_global_equiv(rtt_rng::rn_stream));
    Ensure(rtt_rng::rn_stream == global_ncentot);
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
template<class MT, class PT>
void Rep_Source_Builder<MT,PT>::calc_source_numbers()
{
    Check(parallel_data_op.check_global_equiv(rtt_rng::rn_stream));

    // calculate total, global source energy 
    double total_energy = evoltot + ecentot + esstot;
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

	// estimate the number of global source particles for each source
	// species --> calculations are done on each processor without
	// communication 
	calc_num_src_particles(part_per_e, evol, global_nvol,
			       global_nvoltot); 
	
	calc_num_src_particles(part_per_e, ess, global_nss, global_nsstot); 
	
	calc_num_src_particles(part_per_e, ecen, global_ncen,
			       global_ncentot);

	// sum up source particles
	numtot = global_nvoltot + global_nsstot + global_ncentot;
	
	// ending condition
	if (numtot > npwant && num_try < 10 && np_try_for > 1)
	    np_try_for -= (numtot - npwant);
	else
	    retry = false;
    }
    
    // check global equivalence of source numbers
    Check(parallel_data_op.check_global_equiv(global_nvoltot));
    Check(parallel_data_op.check_global_equiv(global_nsstot));
    Check(parallel_data_op.check_global_equiv(global_ncentot));
    
    // calculate energy weights sampling loss constribution to energy loss
    // based on global source numbers and energies
    for (int cell = 1; cell <= ew_cen.size(); cell++)
    {
	// calculate energy weights of volume source particles per cell
	if (global_nvol(cell) > 0) 
	    ew_vol(cell) = evol(cell) /
		static_cast<double>(global_nvol(cell));
	else 
	    ew_vol(cell) = 0;
	
	// calculate energy weights of surface source particles per cell
	if (global_nss(cell) > 0) 
	    ew_ss(cell) = ess(cell) /
		static_cast<double>(global_nss(cell));
	else 
	    ew_ss(cell) = 0;

	// calculate energy weights of census particles per cell
	if (global_ncen(cell) > 0) 
	    ew_cen(cell) = ecen(cell) /
		static_cast<double>(global_ncen(cell));
	else 
	    ew_cen(cell) = 0;	
	
	// add up sampling contribution to energy loss for volume, ss, and
	// census
	global_eloss_vol += evol(cell) - global_nvol(cell) * ew_vol(cell);
	global_eloss_ss  += ess(cell)  - global_nss(cell)  * ew_ss(cell);
	global_eloss_cen += ecen(cell) - global_ncen(cell) * ew_cen(cell);
    }

    // calculate local source numbers from the global values for each source
    // species
    calc_num_part_and_rn_fields(global_nvol, global_nvoltot,
				rtt_rng::rn_stream, local_nvol,
				local_nvoltot, volrn);
    
    calc_num_part_and_rn_fields(global_nss, global_nsstot,
				rtt_rng::rn_stream, local_nss,
				local_nsstot, ssrn);

    Check(parallel_data_op.check_global_equiv(rtt_rng::rn_stream));
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
 * for a particular species
 * \param global_numtot constant global total number of source particles for
 * a particular species
 * \param next_avail_rn mutable value of next available global random number
 * stream ID
 * \param local_n_field mutable field of local number of source particles for
 * a particular species on processor
 * \param local_numtot mutable value of the local total number of source
 * particles for a particular species on processor
 * \param rn_field mutable field of starting random number stream IDs
 */
template<class MT, class PT> 
void Rep_Source_Builder<MT,PT>::calc_num_part_and_rn_fields
(const ccsf_int &global_n_field, const int global_numtot,
 int &next_avail_rn, ccsf_int &local_n_field, int &local_numtot, 
 ccsf_int &rn_field)
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

    // sweep through cells and determine the number of source particles per
    // cell for this species
    for (int cell = 1; cell <= num_cells; cell++)
    {
	if (global_n_field(cell) > 0)
	{
	    // calculate the integer number of particles per cell
	    even_spread = global_n_field(cell) / nodes();
	    
	    // calculate the remainder left from integer division
	    leftover = global_n_field(cell) - nodes() * even_spread;
	    Check(leftover < nodes());
	    
	    // set the starting random number stream ID for this cell on this
	    // processor. Recall that there are node() processors "below" me
	    rn_field(cell) = next_avail_rn + node() * even_spread +
		rtt_mc::global::min(leftover, node());

	    // set local field; from the remainder (leftover) put one on each
	    // of the first leftover processors
	    local_n_field(cell) = even_spread;
	    if (node() < leftover)
		local_n_field(cell)++;	

	    // increment next_avail_rn number by global number of particles
	    // in this cell
	    next_avail_rn += global_n_field(cell);
	    
	    // update check of total number of source particles
	    local_numtot += local_n_field(cell);
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
