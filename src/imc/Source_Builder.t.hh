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
#include "c4/global.hh"
#include "ds++/Assert.hh"
#include <cmath>

namespace rtt_imc
{

using C4::nodes;
using C4::node;
using C4::Send;
using C4::Recv;
using C4::gsync;
using C4::C4_Req;
using C4::SendAsync;

using std::fabs;

//---------------------------------------------------------------------------//
// SPECIAL GLOBAL EQUIVALENCE FUNCTIONALITY
//---------------------------------------------------------------------------//
/*!
 * \brief Function to check the equivalence of an int across all processors.
 *
 * This function is (hopefully) a temporary parallel check function that more
 * properly belongs in C4.  It is used to check the equivalence of a given
 * integer across all processors.  This is used for Design By Contract
 * analysis in the Source_Builder codes.
 *
 * \param local_value integer value to check against
 * \return true if equivalent across all processors; false if not
 */
template<class MT, class PT>
bool Source_Builder<MT,PT>::check_global_equiv(int local_value) const
{
    // passing condition
    bool pass = false;
    
    // return true if serial, if not then do check on all processors
    if (nodes() == 1)
	pass = true;
    else
    {
	// value from processor above local processor
	int neighbors_value;

	if (node() > 0 && node() < nodes() - 1)
	{
	    Send(local_value, node()-1, 600);
	    Recv(neighbors_value, node()+1, 600);
	    if (local_value == neighbors_value) pass = true;
	}
	else if (node() == nodes() - 1)
	{
	    Send(local_value, node()-1, 600);
	    pass = true;
	}
	else if (node() == 0)
	{
	    Recv(neighbors_value, node()+1, 600);
	    if (local_value == neighbors_value) pass = true;
	}
	else
	{
	    Insist(0, "Something is wrong with nodes!");
	}
    }

    // sync everything so we don't leave before all processors are finished
    gsync();

    // return result
    return pass;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Function to check the equivalence of a double across all processors.
 *
 * This function is the same as check_global_equiv(int) except that doubles
 * are compared to precision eps.
 *
 * \param local_value integer value to check against
 * \param eps precision of double, default 1e-8
 * \return true if equivalent across all processors; false if not 
 */
template<class MT, class PT>
bool Source_Builder<MT,PT>::check_global_equiv(double local_value, 
					       double eps) const 
{
    // passing condition
    bool pass = false;
    
    // return true if serial, if not then do check on all processors
    if (nodes() == 1)
	pass = true;
    else
    {
	// value from processor above local processor
	double neighbors_value;

	if (node() > 0 && node() < nodes() - 1)
	{
	    Send(local_value, node()-1, 600);
	    Recv(neighbors_value, node()+1, 600);
	    if (fabs(local_value - neighbors_value) < eps * local_value) 
		pass = true;
	}
	else if (node() == nodes() - 1)
	{
	    Send(local_value, node()-1, 600);
	    pass = true;
	}
	else if (node() == 0)
	{
	    Recv(neighbors_value, node()+1, 600);
	    if (fabs(local_value - neighbors_value) < eps * local_value) 
		pass = true;
	}
	else
	{
	    Insist(0, "Something is wrong with nodes!");
	}
    }

    // sync everything so we don't leave before all processors are finished
    gsync();

    // return result
    return pass;
}

} // end of rtt_imc

//---------------------------------------------------------------------------//
//                        end of imc/Source_Builder.t.hh
//---------------------------------------------------------------------------//
