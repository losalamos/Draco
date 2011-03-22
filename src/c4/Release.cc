//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/Release.cc
 * \author Thomas M. Evans
 * \date   Thu Jul 15 09:41:09 1999
 * \brief  Release function implementation for c4 library.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_c4
{

using std::string;

/*!
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form c4-\#_\#_\# in pkg_release variable
 */
const string release()
{
    string pkg_release = "c4(draco-6_1_0)";
    return pkg_release;
}

}  // end of rtt_c4

//---------------------------------------------------------------------------//
//                              end of Release.cc
//---------------------------------------------------------------------------//
