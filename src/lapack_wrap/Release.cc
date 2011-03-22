//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   lapack_wrap/Release.cc
 * \author Thomas M. Evans
 * \date   Thu Aug 29 11:06:46 2002
 * \brief  Release function implementation for lapack_wrap library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_lapack_wrap
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form lapack_wrap-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "lapack_wrap(draco-6_1_0)";
    return pkg_release;
}

}  // end of rtt_lapack_wrap

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
