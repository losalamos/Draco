//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Global_Mesh_Data.i.hh
 * \author Thomas M. Evans
 * \date   Thu Dec  4 17:36:07 2003
 * \brief  Global_Mesh_Data member definitions.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

namespace rtt_mc
{

//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * \param t rtt_dsxx::SP to a rtt_mc::Topology class
 * \param mesh mesh type
 */
template<class MT>
Global_Mesh_Data<MT>::Global_Mesh_Data(SP_Topology  t,
				       const MT    &mesh)
    : topology(t)
{
    Require (topology);
    Require (topology->get_parallel_scheme() == "replication" ||
	     topology->get_parallel_scheme() == "DD");

    // calculate the global mesh data
    calc_global_mesh_data(mesh);
}

//---------------------------------------------------------------------------//
// PRIVATE IMPLEMENTATION
//---------------------------------------------------------------------------//
/*!
 * \brief General mesh defintion of global mesh data calculation.
 *
 * This function must be specialized for each mesh type, so it will always
 * throw an assertion if called.
 */
template<class MT>
void Global_Mesh_Data<MT>::calc_global_mesh_data(const MT &mesh)
{
    throw rtt_dsxx::assertion(
	"No specialization provided for MT in Global_Mesh_Data.");
}

} // end namespace rtt_mc

//---------------------------------------------------------------------------//
//                   end of mc/Global_Mesh_Data.i.hh
//---------------------------------------------------------------------------//
