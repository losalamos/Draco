//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/Release.cc
 * \author Thomas M. Evans
 * \date   Fri Aug 24 12:28:30 2001
 * \brief  Release function implementation for cdi_analytic library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_cdi_analytic
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form cdi_analytic-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "cdi_analytic(draco-6_1_0)";
    return pkg_release;
}

}  // end of rtt_cdi_analytic

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
