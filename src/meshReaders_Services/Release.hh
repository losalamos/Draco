//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file   meshReaders_Services/Release.hh
 * \author B.T. Adams
 * \date   Wed June 7 10:33:26 2000
 * \brief  Header file for meshReaders_Services library release function.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_meshReaders_Services_Release_hh
#define rtt_meshReaders_Services_Release_hh

#include <string>

//===========================================================================//
/*!
 * \brief Namespace to contain the meshReaders_Services utilities.
 *
 * Provides namespace protection for the Draco meshReaders_Services utilities.
 *
 * \sa The meshReaders_Services class constructor automatically instantiates
 * and executes the readMesh member function used to parse the mesh
 * data. Accessor functions are provided for all of the remaining member
 * classes to allow data retrieval. The <a href="./index.html"> Main Page
 * </a> presents a summary of the *capabilities provided by the namespace.
 */
//===========================================================================//

namespace rtt_meshReaders_Services 
{

const std::string release();

}  // end of rtt_meshReaders_Services namespace

#endif                          // rtt_meshReaders_Services_Release_hh

//---------------------------------------------------------------------------//
//                              end of meshReaders_Services/Release.hh
//---------------------------------------------------------------------------//
