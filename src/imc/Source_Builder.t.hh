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
#include "Global.hh"
#include "c4/global.hh"
#include "ds++/Assert.hh"
#include <cmath>

namespace rtt_imc
{

using C4::nodes;
using C4::node;
using dsxx::SP;

using std::pow;
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
 * \param interface dsxx::SP to a valid run-time interface
 * \param mesh dsxx::SP to a local mesh object
 * \param top dsxx::SP to a topology object
 */
template<class MT, class PT>
template<class IT>
Source_Builder<MT,PT>::Source_Builder(SP<IT> interface, SP_Mesh mesh, 
				      SP_Topology top)
    : evol_ext(interface->get_evol_ext()),
      rad_s_tend(interface->get_rad_s_tend()),
      rad_source(interface->get_rad_source()),
      rad_temp(interface->get_rad_temp()),
      ss_pos(interface->get_ss_pos()),
      ss_temp(interface->get_ss_temp()),
      defined_surcells(interface->get_defined_surcells()),
      cycle(interface->get_cycle()), delta_t(interface->get_delta_t()), 
      topology(top), 
      census(interface->get_census()),
      ecen(mesh),
      ncen(mesh),
      evol(mesh),
      evol_net(mesh),
      mat_vol_src(mesh),
      evoltot(0),
      mat_vol_srctot(0),
      ess(mesh),
      ss_face_in_cell(mesh),
      esstot(0),
      volrn(mesh),
      ssrn(mesh),
      cenrn(0)
{
    Require(mesh);
    Require(mesh->num_cells() == topology->num_cells(node()));
    
    int num_cells = mesh->num_cells();

    Ensure(evol_ext.size() == num_cells);
    Ensure(rad_source.size() == num_cells);
    Ensure(rad_temp.size() == num_cells);
    Ensure(ss_pos.size() == ss_temp.size());
    Ensure(ss_pos.size() == defined_surcells.size());   
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
	ess.get_Mesh().check_defined_surcells(ss_pos[ss], surcells);

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
void Source_Builder<MT,PT>::calc_ecen_init()
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
 * returns the total number of source particles for the given field type.
 * These calculations are preformed on the local processor.  Typically, one
 * will iterate in a derived source class on the number of source particles
 * across processor space.
 *
 * \param part_per_e desired number of particles per unit energy
 * \param e_field ccsf_double field of energies
 * \param n_field ccsf_int field of number of particles
 * \return total number of particles for this source species
 */
template<class MT, class PT> int
Source_Builder<MT,PT>::calc_num_src_particles(const double part_per_e,
					      const ccsf_double &e_field, 
					      ccsf_int &n_field)
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
    return num_particles;
}

} // end of rtt_imc

//---------------------------------------------------------------------------//
//                        end of imc/Source_Builder.t.hh
//---------------------------------------------------------------------------//
