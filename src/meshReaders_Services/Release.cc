//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file   meshReaders_Services/Release.cc
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Implementation file for meshReaders_Services library release 
 *         function.
 */
//---------------------------------------------------------------------------//
// @> Release function implementation for meshReaders_Servicess library
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace rtt_meshReaders_Services
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
    string pkg_release = "meshReaders_Services(draco-5_4_0)";
    return pkg_release;
}

}  // end of rtt_meshReaders_Servicess namespace

//---------------------------------------------------------------------------//
//                              end of Release.cc
//---------------------------------------------------------------------------//
