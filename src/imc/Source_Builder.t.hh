//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Source_Builder.t.hh
 * \author Thomas M. Evans
 * \date   Wed Dec  8 14:35:43 1999
 * \brief  Source_Builder implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Source_Builder.hh"
#include "Particle.hh"
#include "Global.hh"
#include "ds++/Assert.hh"
#include <cmath>

namespace rtt_imc
{

using rtt_rng::Sprng;
using C4::nodes;
using C4::node;
using dsxx::SP;

using std::pow;
using std::fabs;
using std::vector;
using std::string;
using std::cout;
using std::endl;

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Source_Builder base class constructor.
 *
 * The Source_Builder base class constructor gets data from a run-time
 * interface and a valid Topology to construct its data.  The data in
 * Source_Builder is divided into two basic types: \arg interface data comes
 * from the interface and is used to calculate source fields, \arg source
 * data fields are the products of the Source_Builder process.  The interface
 * data must be given to Source_Builder locally.  That is, the data should be
 * dimensioned to the local (on-processor) mesh size.  All fields in
 * Source_Builder are dimensioned to the local mesh.
 *
 * When the constructor is called the following requirements must be met:
 * \arg a mesh must exist and be properly sized, \arg the rtt_rng::rn_stream
 * must be the same on all processors.
 *
 * \param interface dsxx::SP to a valid run-time interface
 * \param mesh dsxx::SP to a local mesh object
 * \param top dsxx::SP to a topology object 
 */
template<class MT, class PT>
template<class IT>
Source_Builder<MT,PT>::Source_Builder(SP<IT> interface, SP_Mesh mesh, 
				      SP_Topology top)
    : elapsed_t(interface->get_elapsed_t()), 
      evol_ext(interface->get_evol_ext()),
      rad_s_tend(interface->get_rad_s_tend()),
      rad_source(interface->get_rad_source()),
      rad_temp(interface->get_rad_temp()),
      ss_pos(interface->get_ss_pos()),
      ss_temp(interface->get_ss_temp()),
      defined_surcells(interface->get_defined_surcells()),
      ss_desc(interface->get_ss_desc()),
      npnom(interface->get_npnom()),
      npmax(interface->get_npmax()),
      dnpdt(interface->get_dnpdt()),
      cycle(interface->get_cycle()), 
      delta_t(interface->get_delta_t()), 
      ss_dist(interface->get_ss_dist()),
      topology(top), 
      parallel_data_op(top),
      census(interface->get_census()),
      npwant(0),
      ecen(mesh),
      ew_cen(mesh),
      ecentot(0),
      evol(mesh),
      ew_vol(mesh),
      evol_net(mesh),
      mat_vol_src(mesh),
      evoltot(0),
      mat_vol_srctot(0),
      ess(mesh),
      ew_ss(mesh),
      ss_face_in_cell(mesh),
      esstot(0),
      volrn(mesh),
      ssrn(mesh)
{
    using rtt_mc::global::min;

    Require(mesh);
    Require(mesh->num_cells() == topology->num_cells(node()));
    Check(parallel_data_op.check_global_equiv(rtt_rng::rn_stream));
    
    int num_cells = mesh->num_cells();

    // calculate the desired number of source particles
    npwant = min(npmax, static_cast<int>(npnom + dnpdt * elapsed_t)); 

    Ensure(evol_ext.size() == num_cells);
    Ensure(rad_source.size() == num_cells);
    Ensure(rad_temp.size() == num_cells || rad_temp.size() == 0);
    Ensure(ss_pos.size() == ss_temp.size());
    Ensure(ss_pos.size() == defined_surcells.size());
    Ensure(npwant > 0);
}

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
template<class MT, class PT>
void Source_Builder<MT,PT>::calc_source_energies(const Mat_State<MT> &state,
						 const Opacity<MT> &opacity)
{
    // make sure that the number of cells on this mesh is consistent (it
    // should be a submesh)
    Require(topology->num_cells(node()) == state.num_cells());
    Require(topology->num_cells(node()) == opacity.num_cells());

    // calc volume emission energy per cell, total
    calc_evol(state, opacity);

    // calc surface source energy per cell, total
    calc_ess();
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
template<class MT, class PT>
void Source_Builder<MT,PT>::calc_evol(const Mat_State<MT> &state, 
				      const Opacity<MT> &opacity)
{
    // draco necessities
    using rtt_mc::global::a;
    using rtt_mc::global::c;

    // reset evoltot
    evoltot        = 0.0;
    mat_vol_srctot = 0.0;

    // calc volume source and tot volume source
    // evol_net and mat_vol_src needed for temperature update
    for (int cell = 1; cell <= evol.size(); cell++)
    {
	// calc cell centered radiation volume source
	evol_net(cell) = opacity.fplanck(cell) * a * c *
	    pow(state.get_T(cell), 4) * evol.get_Mesh().volume(cell) * 
	    delta_t;
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
		evol_add = rad_source[cell-1] * evol.get_Mesh().volume(cell)
		    * duration; 
		evol(cell) += evol_add;
		evoltot    += evol_add;
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
template<class MT, class PT>
void Source_Builder<MT,PT>::calc_ess()
{
    // draco niceties
    using rtt_mc::global::a;
    using rtt_mc::global::c;

    // reset esstot
    esstot = 0.0;

    // loop over surface sources in problem
    for (int ss = 0; ss < ss_pos.size(); ss++)
    {
        // surface src cells (local id) must be defined by the host
	Check (defined_surcells[ss].size() > 0);
	vector<int> surcells = defined_surcells[ss];
	Check (ss_desc[ss] == "allow_refl_bc" ? true : 
	       ess.get_Mesh().check_defined_surcells(ss_pos[ss], surcells));

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
		pow(ss_temp[ss],4) * delta_t;

	    // accumulate esstot
	    esstot += ess(local_cell);
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
template<class MT, class PT>
void Source_Builder<MT,PT>::calc_initial_ecen()
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
	    pow(rad_temp[cell-1], 4);

	// accumulate evoltot
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
 * the given field type.  These calculations are preformed on the local
 * processor.  Typically, one will iterate in a derived source class on the
 * number of source particles across processor space.
 *
 * \param part_per_e desired number of particles per unit energy
 * \param e_field ccsf_double field of energies
 * \param n_field mutable ccsf_int field of number of particles
 * \param ntot mutable total number of particles for this species 
 */
template<class MT, class PT> void
Source_Builder<MT,PT>::calc_num_src_particles(const double part_per_e,
					      const ccsf_double &e_field, 
					      ccsf_int &n_field,
					      int &ntot)
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
 * \param mesh dsxx::SP to the local mesh
 * \param rcon random number controller
 * \param ncen local field of census particles on processor
 * \param ncentot total number of local census particles
 * \param cenrn local field of census random number IDs
 */
template<class MT, class PT>
void Source_Builder<MT,PT>::write_initial_census(SP_Mesh mesh, 
						 SP_Rnd_Control rcon,
						 const ccsf_int &ncen,
						 const int &ncentot,
						 const ccsf_int &cenrn)
{
    Require(mesh);
    Require(census);
    Require(mesh->num_cells() == ncen.size());
    Require(mesh->num_cells() == cenrn.size());
    Require(mesh->num_cells() == ew_cen.size());

    // loop over all cells
    for (int cell = 1; cell <= mesh->num_cells(); cell++)
	for (int i = 1; i <= ncen(cell); i++)
	{
	    // make a new random number for delivery to Particle
	    Sprng random = rcon->get_rn(cenrn(cell) + i - 1);
	    
	    // sample particle location
	    vector<double> r = mesh->sample_pos(cell, random);

	    // sample particle direction
	    vector<double> omega = mesh->get_Coord().
		sample_dir("isotropic", random);
	    
	    // sample frequency (not now; 1 group)

	    // create Particle
	    SP<PT> particle(new PT(r, omega, ew_cen(cell), cell, random));

	    // write particle to census
	    census->push(particle);
	}

    // a final assertion
    Ensure (census->size() == ncentot);
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
 * \param rcon dsxx::SP to a Rnd_Control object
 * \param local_ncentot mutable value of post-comb total number of census
 * particles on this processor
 * \param eloss_comb mutable value of the energy loss incurred through
 * combing
 */
template<class MT, class PT>
void Source_Builder<MT,PT>::comb_census(SP_Rnd_Control rcon,
					int &local_ncentot,
					double &eloss_comb)
{
    // make double sure we have census particles to comb
    Require(census->size() > 0);
    
    // reset eloss_comb to zero
    eloss_comb = 0;

    // initialize number of (new, combed) census particles
    local_ncentot = 0;

    // in-function data required for combing
    double dbl_cen_part = 0.0;
    int    numcomb      = 0;
    double ecencheck    = 0.0;
    int    cencell      = 0;
    double cenew        = 0.0;

    // local accumulation of census energy from the census particles on this
    // processor 
    double local_ecentot = 0.0;

    // make new census bank to hold combed census particles
    SP_Census comb_census(new Particle_Buffer<PT>::Census());

    // comb census
    while (census->size())
    {
	// read census particle and get cencell, ew
	SP<PT> particle = census->top();
	census->pop();
	cencell      = particle->get_cell();
	cenew        = particle->get_ew();
	Sprng random = particle->get_random();
	Check(cencell > 0 && cencell <= ew_cen.size());

	// add up census energy
	local_ecentot += cenew;
	
	// comb
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
		      	SP<PT> another(new PT(*particle));
			Sprng nran     = rcon->spawn(particle->get_random());
			another->set_random(nran);
			comb_census->push(another);
		    }
		   
		// add up newly combed census particles
		local_ncentot += numcomb;
	    
		// check census energy
		ecencheck  += numcomb * ew_cen(cencell);
	    }
	}
	else
	{
	    // if there is no census energy weight in the cell we sampled
	    // zero census particles previously
	    numcomb = 0;
	
	    // if ewcen == 0 and a census particle exists in this cell then
	    // the energy lost was already tabulated in the energy loss due
	    // to sampling
	}

	// add energy loss to eloss_comb
	eloss_comb += cenew - numcomb * ew_cen(cencell);
    }

    Check(census->size() == 0);
    Check(comb_census->size() == local_ncentot);

    // assign newly combed census to census
    census = comb_census;

    Ensure(census->size() == local_ncentot);
    Ensure(rtt_mc::global::soft_equiv(ecencheck+eloss_comb, local_ecentot, 
				      1.e-6));
}

} // end of rtt_imc

//---------------------------------------------------------------------------//
//                        end of imc/Source_Builder.t.hh
//---------------------------------------------------------------------------//
