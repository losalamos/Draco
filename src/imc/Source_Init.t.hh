//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Source_Init.t.hh
 * \author Thomas M. Evans
 * \date   Tue Jan 11 09:50:00 2000
 * \brief  Source_Init implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Source_Init.hh"

namespace rtt_imc
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Source_Init constructor.
 *
 * \param intrface dsxx::SP to a realization of an interface defined by
 * Interface.hh 
 */
template<class IT, class MT>
Source_Init<IT,MT>::Source_Init(SP_Interface intrface)
    : interface(intrface)
{
    Ensure (interface);
}

//---------------------------------------------------------------------------//
// CALCULATE DEFINED SURFACE CELLS
//---------------------------------------------------------------------------//
/*!
 * \brief Use a global mesh to calculate a list of defined surface cells from
 * a global problem boundary.
 *
 * This function uses a global mesh to convert a global boundary surface
 * source into a list of explicit, global cells along the boundary.  For
 * example, a "lox" boundary could be defined as a surface source.  This
 * description will be converted into a list of all global cells along the
 * low-x boundary.  The list of cells is assigned back into the interface
 * using the Interface::set_defined_surcells function. 
 *
 * Surface source cell lists are only calculated for surface sources that do
 * not already have user-defined lists of cells.
 *
 * In general, this function is run in the initial problem setup on the host
 * node.  However, we do not require that this function be run on the host
 * node.  We only require that the mesh be global, ie. MT::full_Mesh() ==
 * true.
 *
 * \param mesh dsxx::SP to a global mesh
 */
template<class IT, class MT>
void Source_Init<IT,MT>::calc_defined_surcells(SP_Mesh mesh)
{
    Require (mesh->full_Mesh());

    // get surface source positions
    sf_string ss_pos = interface->get_ss_pos();

    // define number of surface sources
    int num_ss = ss_pos.size();

    // exit if no surface sources
    if (num_ss == 0)
	return;

    // get the defined surface source list
    vf_int defined_ss = interface->get_defined_surcells();
    Check (defined_ss.size() == num_ss);

    // check to see if we need to do any work
    for (int i = 0; i < num_ss; i++)
    {
	if (defined_ss[i].size() == 0)
	{
	    // get a vector of surface cells along the boundary
	    defined_ss[i] = mesh->get_surcells(ss_pos[i]);
	    Check (defined_ss[i].size() > 0);

	    // set the new defined surface cell list in the interface
	    interface->set_defined_surcells(i+1, defined_ss[i]);
	}
	else if (defined_ss[i].size() > 0)
	{
	    Check (mesh->check_defined_surcells(ss_pos[i], defined_ss[i]));
	    continue;
	}
	else
	{
	    Insist (0, "Number of defined surface source cells < 0!");
	}
    }
}

} // end of rtt_imc

//---------------------------------------------------------------------------//
//                        end of imc/Source_Init.t.hh
//---------------------------------------------------------------------------//
