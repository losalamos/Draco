//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Release.cc
 * \author Thomas M. Evans
 * \date   Wed Apr 14 15:36:52 1999
 * \brief  Release function implementation for imc library.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_imc
{

using std::string;

/*!
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form imc-\#_\#_\# in pkg_release variable
 */
const string release()
{
    string pkg_release = "@(#)imc-2_2_0b";
    return pkg_release;
}

}  // end of rtt_imc

//---------------------------------------------------------------------------//
//                              end of Release.cc
//---------------------------------------------------------------------------//
