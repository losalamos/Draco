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

namespace rtt_imc
{

using dsxx::SP;

using std::cout;
using std::endl;

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
    : Source_Builder<MT,PT>(interface, mesh, top)
{
    // Nothing Yet
}

//---------------------------------------------------------------------------//
// SOURCE BUILDER MAIN INTERFACE FUNCTION
//---------------------------------------------------------------------------//

template<class MT, class PT>
Rep_Source_Builder<MT,PT>::SP_Source
Rep_Source_Builder<MT,PT>::build_Source(SP_Mesh mesh,
					SP_Mat_State state,
					SP_Opacity opacity,
					SP_Rnd_Control rnd_control)
{
    return SP_Source();
}

//===========================================================================//
// REP_SOURCE_BUILDER IMPLEMENTATION FUNCTIONS
// Private functions in the derived class that are used to build sources.
// These functions are unique to the derived class and usually deal with the
// way that data is parsed across processors.
//===========================================================================//

//---------------------------------------------------------------------------//
// CALCULATE RANDOM NUMBER IDS AND NUMBER OF SOURCE PARTICLES
//---------------------------------------------------------------------------//
/*!  
 * \brief Calculate the local random number stream IDs and the local number
 * of particles per cell per species on the local processor.
 *
 * This function calculates the local number of particles per cell per
 * species on each processor.  We assume that the global number of particles
 * per cell per species has been calculated previously. 
 */
template<class MT, class PT>
void Rep_Source_Builder<MT,PT>::calc_num_part_and_rn_fields()
{
    
}

} // end of rtt_imc

//---------------------------------------------------------------------------//
//                        end of imc/Rep_Source_Builder.t.hh
//---------------------------------------------------------------------------//
