//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Rep_Transporter.t.hh
 * \author Thomas M. Evans
 * \date   Thu Apr 13 11:41:37 2000
 * \brief  Rep_Transporter template definitions.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_imc_Rep_Transporter_t_hh
#define rtt_imc_Rep_Transporter_t_hh

#include "Rep_Transporter.hh"
#include "mc/Particle_Stack.hh"
#include "c4/global.hh"
#include "ds++/Assert.hh"
#include <iostream>
#include <iomanip>

namespace rtt_imc
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Rep_Transporter constructor.
 *
 * Builds an IMC Transporter for full replication topologies.  The
 * constructor checks to make sure that the proper Topology exists for this
 * construction.  Additionally, the communicator must be null.
 */
template<class MT, class FT, class PT>
Rep_Transporter<MT,FT,PT>::Rep_Transporter(SP_Topology top)
    : Transporter<MT,FT,PT>(),
      topology(top)
{
    Require (!mesh);
    Require (!opacity);
    Require (!mat_state);
    Require (!source);
    Require (!tally);
    Require (!random_walk);
    Require (!communicator);
    Require (topology);

    // check the topology
    Insist (topology->get_parallel_scheme() == "replication", 
	    "Invalid topology assigned to Transporter!");
}

//---------------------------------------------------------------------------//
// IMC FULL REPLICATION TRANSPORT
//---------------------------------------------------------------------------//
/*!
 * \brief Do IMC, full replication transport.
 
 * This does F&C IMC transport for one timestep.  The new census created
 * during the timestep is returned.  All fundamental objects are unset at the
 * end of the transport step.

 * \param dt timestep size
 * \param cycle current cycle
 * \param print_f particle notification frequency
 * \param num_to_run total number of particles to run across all processors
 * \param verbose do verbose messaging

 * \return census for this timestep

 */
template<class MT, class FT, class PT>
typename Rep_Transporter<MT,FT,PT>::SP_Census 
Rep_Transporter<MT,FT,PT>::transport(double dt,
				     int    cycle, 
				     int    print_f, 
				     int    num_to_run,
				     bool   verbose) 
{
    using std::cerr;
    using std::cout;
    using std::endl;
    using std::setw;

    Require (ready());    
    
    cerr << ">> Doing transport for cycle " << cycle
	 << " on proc " << C4::node() << " using full replication." << endl;

    // make a new census Comm_Buffer on this node
    SP_Census new_census_bank(new Census());

    // particle history diagnostic
    rtt_dsxx::SP<typename PT::Diagnostic> check;
    if (verbose)
	check = new typename PT::Diagnostic(cout, true); 

    // begin timing the transport on this processor
    double trans_begin = C4::Wtime();

    // get source particles and run them to completion
    int counter = 0;
    while (*source)
    {
	// get a particle from the source
	rtt_dsxx::SP<PT> particle = source->get_Source_Particle(dt); 
	Check (particle->status());

	// transport the particle
	particle->transport(*mesh, *opacity, *tally, random_walk, check);
	counter++;

	// after the particle is no longer active take appropriate action
	Check (!particle->status());
	
	// if census write to file
	if (particle->get_descriptor() == PT::CENSUS)
	    new_census_bank->push(particle);

	// message particle counter
	if (!(counter % print_f)) 
	    cerr << setw(10) << counter << " particles run on proc " 
		 << C4::node() << endl;
    }

    // stop timing the transport on this processor
    double trans_end = C4::Wtime();

    // finished with this timestep
    cerr << ">> Finished transporting " << counter 
         << " particles for cycle "
	 << cycle << " on proc " << C4::node() 
	 << " in " << trans_end - trans_begin << " seconds." << endl;
    
    // unset the transport stuff
    unset();
    
    Ensure (!ready());
    Ensure (new_census_bank);

    return new_census_bank;
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Set fundamental objects for a transport step.
 
 * This function sets the fundamental transport objects: MT, Opacity,
 * Mat_State, Source, Tally, and Communicator for a current timestep.  After
 * the transport function is executed these objects are automatically unset
 * using the unset function.  Thus, the set is valid for \b one execution of
 * transport.

 * The Communicator should be null for Rep_Transporter as a communicator is
 * unnecessary for full replication transport.
 
 * \param mesh_in rtt_dsxx::SP to a fully replicated mesh.
 * \param mat_state_in rtt_dsxx::SP to a valid Mat_State object
 * \param opacity_in rtt_dsxx::SP to a valid Opacity object
 * \param source_in rtt_dsxx::SP to a valid Source object
 * \param tally_in rtt_dsxx::SP to a valid Tally object
 * \param random_walk_in rtt_dsxx::SP to a random walk object (can be null)
 * \param communicator_in null rtt_dsxx::SP to a Communicator (must be null)

 */
template<class MT, class FT, class PT>
void Rep_Transporter<MT,FT,PT>::set(SP_Mesh         mesh_in,
				    SP_Mat_State    mat_state_in,
				    SP_Opacity      opacity_in,
				    SP_Source       source_in,
				    SP_Tally        tally_in,
				    SP_Random_Walk  random_walk_in,
				    SP_Communicator communicator_in)
{
    Require (mesh_in);
    Require (opacity_in);
    Require (source_in);
    Require (mat_state_in);
    Require (tally_in);
    Require (!communicator_in);

    // assign objects (no need to assign communicator as it should be null)
    mesh        = mesh_in;
    opacity     = opacity_in;
    source      = source_in;
    mat_state   = mat_state_in;
    tally       = tally_in;
    random_walk = random_walk_in;
    
    // number of global cells is the same number of cells on processor
    int num_cells = topology->num_cells();

    Ensure (num_cells == mesh->num_cells());
    Ensure (num_cells == opacity->num_cells());
    Ensure (num_cells == source->num_cells());
    Ensure (num_cells == mat_state->num_cells());
    Ensure (num_cells == tally->num_cells());
    Ensure (random_walk ? num_cells == random_walk->num_cells() : true);
    Ensure (topology);
    Ensure (!communicator);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Unset transport objects.

 * This function is used to unset the fundamental transport objects assigned
 * in set.  It is provided in the public interface as a memory conservation
 * feature.  It is automatically called at the end of transport.  We include
 * it in the public interface so that users may optionally unset object.

 * The function does \b not check to see if the fundamental objects are
 * assigned prior to unsetting them.

 */
template<class MT, class FT, class PT>
void Rep_Transporter<MT,FT,PT>::unset()
{
    Require (topology);

    // assign the fundamental objects to null pointers
    mesh         = SP_Mesh();
    opacity      = SP_Opacity();
    mat_state    = SP_Mat_State();
    tally        = SP_Tally();
    source       = SP_Source();
    communicator = SP_Communicator();
    random_walk  = SP_Random_Walk();

    Ensure (!mesh);
    Ensure (!opacity);
    Ensure (!mat_state);
    Ensure (!tally);
    Ensure (!source);
    Ensure (!random_walk);
    Ensure (!communicator);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Query to see if all objects are ready for transport.
 
 * This function checks to make sure that all fundamental IMC transport
 * objects are set in the transporter.  If everything is ready a value of
 * true is returned; otherwise, ready returns false.
 
 */
template<class MT, class FT, class PT>
bool Rep_Transporter<MT,FT,PT>::ready() const
{
    Require (!communicator);
    Require (topology);

    bool indicator = true;

    // if any of the fundamental objects does not exist we are not ready
    if (!mesh) 
	indicator = false;
    else if (!opacity)
	indicator = false;
    else if (!mat_state)
	indicator = false;
    else if (!source) 
	indicator = false;
    else if (!tally)
	indicator = false;

    return indicator;
}

} // end namespace rtt_imc

#endif                          // rtt_imc_Rep_Transporter_t_hh

//---------------------------------------------------------------------------//
//                        end of imc/Rep_Transporter.t.hh
//---------------------------------------------------------------------------//
