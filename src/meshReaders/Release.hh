//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file   meshReaders/Release.hh
 * \author B.T. Adams
 * \date   Fri Aug 27 10:33:26 1999
 * \brief  Header file for meshReaders library release function.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_meshReaders_Release_hh
#define rtt_meshReaders_Release_hh

//===========================================================================//
/*!

 * \namespace rtt_meshReaders

 * \brief Namespace to contain the RTT mesh Reader utilities.
 *
 * Provides namespace protection for the Draco (RTT) mesh reader
 * utilities.

 * The Element_Definition class describes the mesh elements supported. The
 * RTT_Format class provides a capability to read RTT format mesh files. The
 * Hex_Mesh_Reader class provides a capability to read CIC-19 Hex Mesh format
 * mesh files.
 
 * \sa The \ref index presents a summary of the capabilities provided within
 * the namespace.

 */
//===========================================================================//

#include <string>

namespace rtt_meshReaders 
{

/*!
 * \brief  Gets the release number for the meshReaders package. 
 * \return release number as a string in the form meshReaders-\#_\#_\#
 */
const std::string release();

}  // end of rtt_meshReaders namespace

#endif                          // rtt_meshReaders_Release_hh
//---------------------------------------------------------------------------//
//                              end of Release.hh
//---------------------------------------------------------------------------//
