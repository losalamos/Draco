//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Tally_Builder.t.hh
 * \author Thomas M. Evans
 * \date   Thu Aug 21 09:29:35 2003
 * \brief  Tally_Builder member definitions.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_imc_Tally_Builder_t_hh
#define rtt_imc_Tally_Builder_t_hh

#include "Tally_Builder.hh"
#include "Tally.hh"

namespace rtt_imc
{

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Build a tally.
 *
 * \return an instantiated tally with the appropriate sub tallies as
 * determined by the interface
 */
template<class MT>
typename Tally_Builder<MT>::SP_Tally 
Tally_Builder<MT>::build_Tally(SP_Mesh mesh)
{
    Require (mesh);

    // build the tally
    SP_Tally tally(new Tally<MT>(mesh));

    // assign sub tallies (these will be null if not required by the
    // interface)
    tally->assign_RW_Sub_Tally(rw_sub_tally);
    tally->assign_Surface_Sub_Tally(sur_sub_tally);

    Ensure (tally);
    return tally;
}

} // end namespace rtt_imc

#endif // rtt_imc_Tally_Builder_t_hh

//---------------------------------------------------------------------------//
//                   end of imc/Tally_Builder.t.hh
//---------------------------------------------------------------------------//
