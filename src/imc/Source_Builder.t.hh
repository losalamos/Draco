//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Source_Builder.t.hh
 * \author Thomas M. Evans and Todd Urbatsch
 * \date   Wed Dec  8 14:35:43 1999
 * \brief  Source_Builder implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Source_Builder.hh"
#include "Gray_Particle.hh"
#include "Multigroup_Particle.hh"
#include "Frequency.hh"
#include "rng/Random.hh"
#include "mc/Sampler.hh"

namespace rtt_imc
{

//===========================================================================//
// IMPLEMENTATION INHERITANCE FOR SOURCE BUILDERS
// The functions in the Source_Builder implementation inheritance are
// protected.  They access protected data.  Thus, in the derived classes
// they become protected (private) functionality.  Because these functions
// are used by all Source_Builder derived classes, the data they modify is
// local to a processor.  For example, calc_source_energies calculates the
// on-processor source energies for "DD" and "REPLICATION" topologies.
//
// Associate functions that are called inside of the protected interface are
// often private.  We do this because this is a pure implementation detail.
// The derived classes never need to see these functions.  We include private 
// implementation functions beneath the protected interface functions from
// which they were called.
//===========================================================================//

//---------------------------------------------------------------------------//
// FUNCTIONS THAT BUILD SOURCE ENERGIES
//---------------------------------------------------------------------------//
/*!  
 * \brief Calculate source energies for volume emission and surface sources
 * on local processor.
 *
 * This function fills up the evol, evol_net, ess, and mat_vol_src fields.
 * The data calculated by this function is local.  That is, no communication
 * is done yet to determine the total source energies across all processors.
 *
 * \param state local (on-processor) Mat_State object
 * \param opacity local (on-processor) Opacity object */
template<class MT, class FT, class PT>
void Source_Builder<MT,FT,PT>::calc_source_energies(
    const Mat_State<MT>  &state,
    const Opacity<MT,FT> &opacity)
{
    // make sure that the number of cells on this mesh is consistent (it
    // should be a submesh)
    Require(topology->num_cells(C4::node()) == state.num_cells());
    Require(topology->num_cells(C4::node()) == opacity.num_cells());

    // calc volume emission energy per cell, total
    calc_evol(state, opacity);

    // calc surface source energy per cell, total
    calc_ess(opacity);

    // set the surface source temperature in the Source object so, in the
    // case of multigroup frequency treatment, the Source can sample
    // frequency group from a Planckian.

}

//---------------------------------------------------------------------------//
/*!  
 * \brief Calculate source energies for volume emission sources on local
 * processor.
 *
 * This function is used by the Source_Builder base class inside of
 * calc_source_energies to calculate volume emission, net volume emission,
 * and external material volume emission source energies.  As such, it is a
 * private function because no derived classes will call this function
 * directly.
 *
 * \param state local (on-processor) Mat_State object
 * \param opacity local (on-processor) Opacity object
 */
template<class MT, class FT, class PT>
void Source_Builder<MT,FT,PT>::calc_evol(const Mat_State<MT>  &state, 
					 const Opacity<MT,FT> &opacity)
{
    // draco necessities
    using rtt_imc::global::Type_Switch;

    // reset evoltot
    evoltot        = 0.0;
    mat_vol_srctot = 0.0;

    // calc evol net depending on frequency type
    calc_evol_net(Type_Switch<FT>(), state, opacity);
    
    // calc volume source and tot volume source
    // evol_net and mat_vol_src needed for temperature update
    for (int cell = 1; cell <= evol.size(); cell++)
    {
	// calc cell centered radiation volume source

	// add time-implicit portion of material energy source
	evol(cell) = evol_net(cell) + 
	    evol_ext[cell-1] * (1.0 - opacity.get_fleck(cell)) *  
	    evol.get_Mesh().volume(cell) * delta_t;

	// accumulate evoltot
	evoltot += evol(cell);

	// calc cell centered material volume source
	mat_vol_src(cell) = opacity.get_fleck(cell) * evol_ext[cell-1] *
	    evol.get_Mesh().volume(cell) * delta_t; 

	// accumulate mat_vol_srctot
	mat_vol_srctot += mat_vol_src(cell);
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
	    double evol_add;
	    for (int cell = 1; cell <= evol.size(); cell++)
	    {
		// calculate radiation source
		evol_add = rad_source[cell-1] * evol.get_Mesh().volume(cell)
		    * duration; 

		// add the radiation source to evol
		evol(cell) += evol_add;
		evoltot    += evol_add;

		// calculate fraction of evol_net to evol
		calc_prob_Planck_emission<Dummy_Type>(Type_Switch<FT>(), 
						      evol_add, cell);
	    }
	}
    }
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Calculate surface source energies on local processor.
 *
 * This function calculates local surface source energies in cells that
 * contain a surface source.  It uses the defined_surcells data from a
 * run-time interface.  This data gives a list of local cells that contain a
 * surface source.  Thus, this function is domain independent because the
 * surface source cells are listed individually.  
 */
template<class MT, class FT, class PT>
void Source_Builder<MT,FT,PT>::calc_ess(const Opacity<MT,FT> &opacity)
{
    // draco niceties
    using rtt_mc::global::a;
    using rtt_mc::global::c;
    using std::vector;
    using rtt_imc::global::Type_Switch;

    // reset esstot
    esstot = 0.0;

    // loop over surface sources in problem
    for (int ss = 0; ss < ss_pos.size(); ss++)
    {
        // surface src cells (local id) must be defined by the host
	vector<int> surcells = defined_surcells[ss];
	
	// do the surface source calculation only if there are actually cells 
	// on this processor that contain the surface source
	if (!surcells.empty())
	{
	    Check (ss_desc[ss] == "allow_refl_bc" ? true : 
		   ess.get_Mesh().check_defined_surcells(ss_pos[ss],
							 surcells));
	    
	    int local_cell;
	    for (int sc = 0; sc < surcells.size(); sc++)
	    {      
		// local cell index for the ss'th surface source
		local_cell = surcells[sc];
		Check (local_cell > 0 && local_cell <=
		       ss_face_in_cell.get_Mesh().num_cells()); 
		
		// make sure this cell doesn't already have a surface source
		Check (ss_face_in_cell(local_cell) == 0);
		
		// assign source face to surface source cell
		ss_face_in_cell(local_cell) = ss_face_in_cell.get_Mesh().
		    get_bndface(ss_pos[ss], local_cell);
		
		// assign energy to surface source cell
		ess(local_cell) = a * c * 0.25 *
		    ess.get_Mesh().face_area
		    (local_cell, ss_face_in_cell(local_cell)) * 
		    std::pow(ss_temp[ss],4) * 
		    opacity.get_integrated_norm_Planck(ss_temp[ss]) * delta_t;

		// accumulate esstot
		esstot += ess(local_cell);

		// set the surface source temperature in the source object's
		// Frequency_Sampling_Data struct
		set_ss_temperature_in_source<Dummy_Type>(
		    Type_Switch<FT>(), ss_temp[ss], local_cell);
	    }
	}
    }
}

//---------------------------------------------------------------------------//
/*!  
 * \brief Calculate source energies for the initial census on the local
 * processor.
 *
 * This function is used to calculate the \b initial census energies on the
 * local processor.  
 */
template<class MT, class FT, class PT>
void Source_Builder<MT,FT,PT>::calc_initial_ecen(
    const Opacity<MT,FT> &opacity)
{
    // draco stuff
    using rtt_mc::global::a;

    // reset ecentot
    ecentot = 0.0;

    // calc census radiation energy in each cell and accumulate
    for (int cell = 1; cell <= ecen.size(); cell++)
    {
	// calc cell centered census radiation energy
	ecen(cell) = a * ecen.get_Mesh().volume(cell) *
	    std::pow(rad_temp[cell-1], 4) * 
	    opacity.get_integrated_norm_Planck(rad_temp[cell-1]);

	// accumulate ecentot on-processor
	ecentot += ecen(cell);
    }
}

//---------------------------------------------------------------------------//
// CALCULATE NUMBER OF SOURCE PARTICLES FOR A GIVEN SOURCE FIELD
//---------------------------------------------------------------------------//
/*!
 * \brief Calculate the number of source particles for a given source field
 * on the local processor.
 *
 * This functions takes a field of source energies and the desired number of
 * particles per unit energy, and it fills in the source numbers field and
 * returns (as a mutable argument) the total number of source particles for
 * the given field type.  These calculations are performed on the local
 * processor.  Typically, one will iterate in a derived source class on the
 * number of source particles across processor space.
 *
 * \param part_per_e desired number of particles per unit energy
 * \param e_field ccsf_double field of energies
 * \param n_field mutable ccsf_int field of number of particles
 * \param ntot mutable total number of particles for this species 
 */
template<class MT, class FT, class PT> 
void Source_Builder<MT,FT,PT>::calc_num_src_particles(
    const double       part_per_e,
    const ccsf_double &e_field, 
    ccsf_int          &n_field,
    int               &ntot)
{
    Require(e_field.size() == n_field.size());
    Require(part_per_e >= 0);

    // return value of total number of particles on this processor
    int num_particles = 0;

    // (double) valued number of particles per cell
    double d_num;

    // sweep through cells and calculate number of particles per cell
    for (int cell = 1; cell <= e_field.size(); cell++)
    {
	// if the cell has any energy try to put some particles in it
	if (e_field(cell) > 0.0)
	{
	    // get estimate of number of particles per cell to nearest
	    // integer per species
	    d_num = e_field(cell) * part_per_e;
	    n_field(cell) = static_cast<int>(d_num + 0.5);

	    // try to get at least one particle per cell per species
	    if (n_field(cell) == 0)
 		n_field(cell) = static_cast<int>(d_num + 0.9999);
	    
	    // increment particle counter
	    num_particles += n_field(cell);
	}
	else
	    n_field(cell) = 0;
    }

    // conditions on the number of particles
    Ensure(num_particles >= 0);
    
    // return the total number of particles
    ntot = num_particles;
}

//---------------------------------------------------------------------------//
// WRITE INITIAL CENSUS
//---------------------------------------------------------------------------//
/*!
 * \brief Calculate the initial census on a local processor.
 *
 * This function takes local initial census data calculated in the derived
 * source builder classes and writes an intial census.  It is used in the
 * derived class implementations of calc_initial_census().  All source
 * builder derived classes utilize this function; however, each derived class
 * must provide the correct local census fields.
 *
 * \param mesh rtt_dsxx::SP to the local mesh
 * \param rcon random number controller
 * \param ncen local field of census particles on processor
 * \param ncentot total number of local census particles
 * \param cenrn local field of census random number IDs
 */
template<class MT, class FT, class PT>
void Source_Builder<MT,FT,PT>::write_initial_census(
    SP_Mesh              mesh, 
    SP_Rnd_Control       rcon,
    const FT            &frequency,
    const Mat_State<MT> &mat_state,
    const ccsf_int      &ncen,
    const int           &ncentot,
    const ccsf_int      &cenrn)
{
    using rtt_imc::global::Type_Switch;
    using rtt_rng::Sprng;
    using rtt_dsxx::SP;
    using std::vector;

    Require(mesh);
    Require(census);

    Check(mesh->num_cells() == ncen.size());
    Check(mesh->num_cells() == cenrn.size());
    Check(mesh->num_cells() == ew_cen.size());

    // loop over all cells
    for (int cell = 1; cell <= mesh->num_cells(); cell++)
    {
	// calculate the global index for this cell
	int global_cell = topology->global_cell(cell);

	// loop over the locally desired number of census particles 
	for (int i = 1; i <= ncen(cell); i++)
	{
	    // make a new random number for delivery to Particle
	    int rn_str_id = INTEGER_MODULO_1E9(cenrn(cell) + i-1); 
	    Sprng random  = rcon->get_rn(rn_str_id);
	    
	    // sample particle location
	    vector<double> r = mesh->sample_pos(cell, random);

	    // sample particle direction
	    vector<double> omega = mesh->get_Coord().
		sample_dir("isotropic", random);

	    // create Particle
	    SP<PT> particle = make_particle<Dummy_Type>(
		Type_Switch<FT>(), Type_Switch<PT>(), frequency,
		mat_state.get_T(cell), r, omega, ew_cen(cell),
		global_cell, random);

	    // write particle to census
	    census->push(particle);
	}
    }

    // a final assertion
    Ensure (census->size() == ncentot);
}

//---------------------------------------------------------------------------//
// Specialization of write census for Gray_Particle type.

template<class MT, class FT, class PT>
template<class Stop_Explicit_Instantiation>
rtt_dsxx::SP<Gray_Particle<MT> >
Source_Builder<MT,FT,PT>::make_particle(
    rtt_imc::global::Type_Switch<Gray_Frequency>,
    rtt_imc::global::Type_Switch<Gray_Particle<MT> >,
    Gray_Frequency        frequency,
    const double          T,
    const sf_double      &r,
    const sf_double      &omega,
    const double          ew,
    int                   global_cell,
    const rtt_rng::Sprng &random)
{
    using rtt_dsxx::SP;

    SP<Gray_Particle<MT> > particle(
	new Gray_Particle<MT>(r, omega, ew, global_cell, random));

    Check (particle);
    return particle;
}

//---------------------------------------------------------------------------//
// Specialization of write census for Multigroup_Particle type.

template<class MT, class FT, class PT>
template<class Stop_Explicit_Instantiation>
rtt_dsxx::SP<Multigroup_Particle<MT> > 
Source_Builder<MT,FT,PT>::make_particle(
    rtt_imc::global::Type_Switch<Multigroup_Frequency>,
    rtt_imc::global::Type_Switch<Multigroup_Particle<MT> >,
    Multigroup_Frequency  frequency,
    const double          T,
    const sf_double      &r,
    const sf_double      &omega,
    const double          ew,
    int                   global_cell,
    const rtt_rng::Sprng &random)
{
    using rtt_mc::sampler::sample_planckian_frequency;
    using rtt_dsxx::SP;

    // sample group from Planckian
    int group   = 0;
    int counter = 0;
    double freq = 0.0;

    while (group == 0)
    {
	freq   = sample_planckian_frequency(random, T);
	group  = frequency.find_group_given_a_freq(freq);

	// advance counter
	counter++;

	Insist (counter < 100, "Unable to sample group.");
    }

    SP<Multigroup_Particle<MT> > particle(
	new Multigroup_Particle<MT>(r, omega, ew, global_cell, random, group));

    Check (particle);
    return particle;
}

//---------------------------------------------------------------------------//
// COMB CENSUS
//---------------------------------------------------------------------------//
/*!
 * \brief Comb a local census object on each processor.
 *
 * This functions performs a comb on the on-processor, local census object.
 * No communication is done in the comb.  The client must ensure that the
 * (global) value of the desired energy weight per cell has been previously
 * calculated.  The function returns the on-processor, local energy loss
 * resulting from the comb.
 *
 * \param rcon rtt_dsxx::SP to a Rnd_Control object
 * \param local_ncen mutable cell values of post-comb total number of census
 *                   particles on this processor
 * \param local_ncentot mutable value of post-comb total number of census
 *                      particles on this processor.  On entry, after the
 *                      initial cycle, it is the total value of desired 
 *                      census particles on this processor
 * \param eloss_comb mutable value of the energy loss incurred through
 *                   combing
 */
template<class MT, class FT, class PT>
void Source_Builder<MT,FT,PT>::comb_census(SP_Rnd_Control  rcon,
					   ccsf_int       &local_ncen,
					   int            &local_ncentot,
					   double         &eloss_comb,
					   SP_Census       dead_census,
					   ccsf_int       &max_dead_rand_id)
{
    using rtt_rng::Sprng;
    using rtt_dsxx::SP;

    // silly checks
    Require(census);
    Require(local_ncentot >= 0);
    
    // reset eloss_comb to zero
    eloss_comb = 0;

    // initialize number of (new, combed) census particles
    std::fill(local_ncen.begin(), local_ncen.end(), 0);
    local_ncentot = 0;

    // in-function data required for combing
    double dbl_cen_part      = 0.0;
    int    numcomb           = 0;
    double ecencheck         = 0.0;
    int    cencell           = 0;
    double cenew             = 0.0;
    int    local_cell        = 0;
    double local_eloss_check = 0.0;

    // local accumulation of census energy from the census particles on this
    // processor 
    double local_precombed_ecentot = 0.0;

    // make new census bank to hold combed census particles
    SP_Census comb_census(new Census());

    // comb census if there are particles in it
    if (census->size() > 0)
    {
	while (census->size())
	{
	    // read census particle
	    rtt_dsxx::SP<PT> particle = census->top();
	    census->pop();
	    
	    // get pertinent census particle attributes
	    cencell      = particle->get_cell();
	    cenew        = particle->get_ew();
	    Sprng random = particle->get_random();
	    
	    // get local cell corresponding to census particle's global cell
	    Check (cencell > 0);
	    local_cell = topology->local_cell(cencell);
	    Check (local_cell > 0 && local_cell <= ew_cen.size());
	    
	    // add up census energy
	    local_precombed_ecentot += cenew;
	    
	    // comb
	    if (ew_cen(local_cell) > 0)
	    {
		dbl_cen_part = (cenew / ew_cen(local_cell)) + random.ran();
		numcomb = static_cast<int>(dbl_cen_part);
		
		// create newly combed census particles
		if (numcomb > 0)
		{
		    particle->set_ew(ew_cen(local_cell));
		    comb_census->push(particle);
		    
		    if (numcomb > 1)
			for (int nc = 1; nc <= numcomb-1; nc++)
			{
			    // COPY a new particle and spawn a new RN state
			    SP<PT> another(new PT(*particle));
			    Sprng nran = rcon->spawn(particle->get_random());
			    another->set_random(nran);
			    comb_census->push(another);
			}
		    
		    // add up newly combed census particles
		    local_ncen(local_cell) += numcomb;
		    local_ncentot          += numcomb;
		    
		    // check census energy total
		    ecencheck  += numcomb * ew_cen(local_cell);
	    
		}
		// put the combed-out particle into the dead_census
		else
		{
		    if (random.get_num() > max_dead_rand_id(local_cell))
			max_dead_rand_id(local_cell) = random.get_num();
		    dead_census->push(particle);
		}

		// add energy loss to eloss_comb (returned argument)
		eloss_comb += cenew - numcomb * ew_cen(local_cell);
	    }
	    else
	    {
		// if there is no census energy weight in the cell we sampled
		// zero census particles previously
		numcomb = 0;
		
		// if ewcen == 0 and a census particle exists in this cell
		// then the energy lost was already tabulated in the energy
		// loss due to sampling in function calc_source_numbers
		// and/or calc_initial_ncen
	    }

	    // perform a local sum of the total, local energy loss.  This sum 
	    // will include loss due sampling, which is accounted for outside 
	    // this function in global_eloss_cen.
	    local_eloss_check += cenew - numcomb * ew_cen(local_cell);
	}
	
	Check(census->size() == 0);
	Check(comb_census->size() == local_ncentot);
	
	// assign newly combed census to census
	census = comb_census;
    }

    // check consistency of census size and external count (should be the
    // same after combing); energy checks and balances.
    Ensure(census->size() == local_ncentot);
    Ensure(rtt_mc::global::soft_equiv(ecencheck + local_eloss_check,
				      local_precombed_ecentot,  
				      (local_ncentot+1) * 1.0e-12));
}

//---------------------------------------------------------------------------//
// RESET THE PARTICLE EW'S IN THE CENSUS
//---------------------------------------------------------------------------//
/*!
 * \brief Go through the census particles on this processor and reset the
 * energy-weights. 
 *
 * Once we have combed the census particles, we recalculate the census
 * particle energy-weights using recalc_census_ew_after_comb, and, in this
 * function, we go through each particle in this processor's census and reset
 * its energy-weight, ew, according to ew_cen.  This post-comb readjustment
 * of the energy-weights better conserves energy at the cost of increasing
 * the variance somewhat. This readjustment is necessary because of the
 * nature of the reproducible comb.  This function also takes the global
 * census energy loss due to combing, global_eloss_comb, as an argument and
 * updates it.  Since the census is always local, this function applies to
 * any topology.
 *
 * \param local_ncentot         const value of post-comb total number of
 *                              census particles on this processor--used as a
 *                              check 
 * \param global_eloss_comb     mutable value of the energy loss in the
 *                              census due only to combing
 * \param global_ecentot_actual const total value of global, actual, 
 *                              pre-combed census energy, used as a check
 *
 */
template<class MT, class FT, class PT>
void Source_Builder<MT,FT,PT>::reset_ew_in_census(
    const int    &local_ncentot,
    double       &global_eloss_comb,
    const double &global_ecentot_actual) 
{
    using rtt_mc::global::soft_equiv;

    // check consistency of the census size and local_ncentot
    Require(census->size() == local_ncentot);

    // calculate global ncentot inside the scope of this function
    int global_ncentot = local_ncentot;
    C4::gsum(global_ncentot);

    // make a new census bank to hold census particle with updated ew's
    SP_Census updated_census(new Census());

    // initialize local, global energy loss due to the ew readjustment.
    // hopefully the global_incremental_eloss takes global_eloss_cen closer
    // to zero!
    double incremental_eloss        = 0.0;
    double global_incremental_eloss = 0.0;
    double check_ecentot            = 0.0;
    double check_global_ecentot     = 0.0;

    // update census particles' energy-weights
    while (census->size())
    {
	// get particle off census list
	rtt_dsxx::SP<PT> particle = census->top();
	census->pop();

	// get census particle's cell (always global); convert to local
	int cencell = particle->get_cell();
	Check (cencell > 0);
	int local_cell = topology->local_cell(cencell);
	Check (local_cell > 0);
	Check (local_cell <= num_cells());

	// sum up the energy change due to ew readjustment
	incremental_eloss += (particle->get_ew() - ew_cen(local_cell));

	// readjust the particle's energy-weight
	particle->set_ew(ew_cen(local_cell));

	// accumulate on-processor ecentot check
	check_ecentot += ew_cen(local_cell);

	// push particle into the updated census 
	updated_census->push(particle);
    }

    // check for appropriate sizes
    Check (census->size() == 0);
    Check (updated_census->size() == local_ncentot);

    // reassign the census pointer
    census = updated_census;

    // sum up all processors' energy loss due to ew readjustment
    global_incremental_eloss = incremental_eloss;
    C4::gsum(global_incremental_eloss);

    // update (presumably decrease) the global census energy loss due to
    // combing 
    global_eloss_comb += global_incremental_eloss;

    // sum up all processors' check_ecentot and test 
    check_global_ecentot = check_ecentot;
    C4::gsum(check_global_ecentot);

    // the "reference" and "value to be compared" have been reversed in this
    // check on equivalence.  Here, the "value to be compared"
    // (check_global_ecentot+global_eloss_cen) can be identically zero if the
    // census is entirely unsampled.  The soft_equiv function can rigorously
    // treat the case when the "reference"==0; the converse is not
    // necessarily true.  Thus, we reverse the arguments to soft_equiv.
    Check ((global_ecentot_actual > 1.0e-14 &&
	    check_global_ecentot+global_eloss_comb > 1.0e-14) ?
	   soft_equiv(global_ecentot_actual,
		      check_global_ecentot + global_eloss_comb, 
		      static_cast<double>(global_ncentot+1)*1.0e-14) :
	   true);
}

} // end of rtt_imc

//---------------------------------------------------------------------------//
//                        end of imc/Source_Builder.t.hh
//---------------------------------------------------------------------------//
