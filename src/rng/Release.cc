//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file   Release.cc
 * \author Thomas M. Evans
 * \date   Thu May 27 15:24:02 1999
 * \brief  Release function implementation for rng library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_rng
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form rng-#_#_# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "@(#)rng-1_0_0x";
    return pkg_release;
}

}  // end of rtt_rng

//---------------------------------------------------------------------------//
//                              end of Release.cc
//---------------------------------------------------------------------------//
