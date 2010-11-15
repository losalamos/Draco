//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file   ds++/dim.hh
 * \author Kent G. Budge
 * \date   Wed Jan 22 15:18:23 MST 2003
 * \brief  Replacement for the Fortran dim function
 */
//---------------------------------------------------------------------------//
// $Id$ 
//---------------------------------------------------------------------------//

#ifndef rtt_dsxx_dim_hh
#define rtt_dsxx_dim_hh

#include <algorithm>
#include <iterator>
#include <functional>
#include "ds++/Assert.hh"

namespace rtt_dsxx
{

//! Return the positive difference of the arguments.
template <class Ordered_Group_Element>
inline Ordered_Group_Element dim(Ordered_Group_Element a, 
				 Ordered_Group_Element b)
{
    if (a<b)
    {
        return Ordered_Group_Element(0);
    }
    else
    {
        return a-b;
    }
}

} // ane of namespace rtt_dsxx

#endif // rtt_dsxx_dim_hh
//---------------------------------------------------------------------------//
//                           end of dim.hh
//---------------------------------------------------------------------------//



