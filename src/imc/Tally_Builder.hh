//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Tally_Builder.hh
 * \author Thomas M. Evans
 * \date   Thu Aug 21 09:29:35 2003
 * \brief  Tally_Builder class definition.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_imc_Tally_Builder_hh
#define rtt_imc_Tally_Builder_hh

#include <vector>
#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include "Azimuthal_Mesh.hh"
#include "Surface_Sub_Tally.hh"
#include "Random_Walk_Sub_Tally.hh"
#include "Hybrid_Diffusion.hh"

namespace rtt_imc
{

template<class MT> class Tally;

//===========================================================================//
/*!
 * \class Tally_Builder
 * \brief Build a Tally object.
 *
 * \sa Tally_Builder.t.hh for member definitions.
 *
 * Code Sample:
 * \code
 *     Tally_Builder<MT> builder(interface);
 *     SP<Tally<MT> > tally = builder.build_Tally(mesh);
 * \endcode
 */
// revision history:
// -----------------
// 0) (Thu Aug 21 09:29:35 2003) Thomas M. Evans: original
// 
//===========================================================================//

template<class MT>
class Tally_Builder 
{
  public:
    typedef rtt_dsxx::SP<MT>                    SP_Mesh;
    typedef rtt_dsxx::SP<Tally<MT> >            SP_Tally;
    typedef rtt_dsxx::SP<Random_Walk_Sub_Tally> SP_RW_Sub_Tally;
    typedef rtt_dsxx::SP<Surface_Sub_Tally>     SP_Surface_Sub_Tally;
    typedef rtt_dsxx::SP<Azimuthal_Mesh>        SP_Azimuthal_Mesh;
    
  private:
    // >>> DATA

    // Smart pointer to the surface sub tally.
    std::vector<SP_Surface_Sub_Tally> surface_sub_tallies;

    // Smart pointer to randomw walk sub tally.
    SP_RW_Sub_Tally rw_sub_tally;

    // Number of energy groups.
    int energy_groups;

  public:
    // Constructor.
    template<class IT>
    Tally_Builder(rtt_dsxx::SP<IT>, int groups = 1);

    // >>> PUBLIC INTERFACE

    // Build the tally.
    SP_Tally build_Tally(SP_Mesh);
};

//---------------------------------------------------------------------------//
// MEMBER TEMPLATE DEFINITIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for the Tally_Builder
 *
 * The interface must be a model of rtt_imc::Interface and a derived class of
 * rtt_imc::Surface_Tracking_Interface.
 *
 * We build a surface sub tally if there are global tally surfaces defined.
 * There may be no local tally surfaces; the Surface_Sub_Tally will still be
 * defined. 
 *
 * We define a Random_Walk_Sub_Tally if random walk is on.
 *
 * \param interface rtt_dsxx::SP to an interface that implements (is a
 * realization) of rtt_imc::Interface and is a derived class of
 * rtt_imc::Surface_Tracking_Interface.
 */
template<class MT>
template<class IT>
Tally_Builder<MT>::Tally_Builder(rtt_dsxx::SP<IT> interface, int groups)
    : energy_groups(groups)
{
    Require (interface);
    Require (energy_groups > 0);

    // build the surface sub tally
    if (interface->number_of_surfaces())
    {
	// build an Azimuthal_Mesh
	SP_Azimuthal_Mesh az_mesh(new Azimuthal_Mesh(*interface));

	// build the surface sub tallies
	surface_sub_tallies.resize(energy_groups);

	for (int group = 0; group < groups; ++group)
	{
	    surface_sub_tallies[group] = new Surface_Sub_Tally(
		az_mesh, *interface);
	    Check (surface_sub_tallies[group]);
	}
    }

    // build the random walk sub tally
    if (interface->get_hybrid_diffusion_method() == 
	Hybrid_Diffusion::RANDOM_WALK)
    {
	rw_sub_tally = new Random_Walk_Sub_Tally();

	Check (rw_sub_tally);
    }	
}

} // end namespace rtt_imc

#endif // rtt_imc_Tally_Builder_hh

//---------------------------------------------------------------------------//
//              end of imc/Tally_Builder.hh
//---------------------------------------------------------------------------//
