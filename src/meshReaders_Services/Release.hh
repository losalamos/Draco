//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file   meshReaders_Services/Release.hh
 * \author B.T. Adams
 * \date   Wed June 7 10:33:26 2000
 * \brief  Header file for meshReaders_Services library release function.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __meshReaders_Services_Release_hh__
#define __meshReaders_Services_Release_hh__

//===========================================================================//
// namespace version - 
//
// Purpose : Return the version of meshReaders_Services; this can be used to 
//           get exact version information in codes that use 
//           meshReaders_Services
// 
//===========================================================================//

#include <string>

/*!
 * \brief Namespace to contain the meshReaders_Services utilities.
 *
 * Provides namespace protection for the Draco meshReaders_Services utilities.
 *
 *\sa The meshReaders_Services class constructor automatically instantiates 
 *    and executes the readMesh member function used to parse the mesh data. 
 *    Accessor functions are provided for all of the remaining member classes 
 *    to allow data retrieval. The \ref rtt_mesh_reader_overview page presents
 *    a summary of the capabilities provided by the namespace.
 */
namespace rtt_meshReaders_Services 
{
/*!
 * \brief  Gets the release number for the meshReaders_Services package. 
 * \return release number as a string in the form 
 *         "meshReaders_Services-\#_\#_\#"
 */
const std::string release();

}  // end of rtt_meshReaders_Services namespace

#endif                          // __meshReaders_Services_Release_hh__

//---------------------------------------------------------------------------//
//                              end of meshReaders_Services/Release.hh
//---------------------------------------------------------------------------//
