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
 * \param intrface rtt_dsxx::SP to a realization of an interface defined by
 * Interface.hh 
 */
template<class IT, class MT>
Source_Init<IT,MT>::Source_Init(SP_Interface intrface)
    : interface(intrface)
{
    Ensure (interface);
}

} // end of rtt_imc

//---------------------------------------------------------------------------//
//                        end of imc/Source_Init.t.hh
//---------------------------------------------------------------------------//
