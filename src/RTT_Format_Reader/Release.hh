//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   RTT_Format_Reader/Release.hh
 * \author Thomas M. Evans
 * \date   Mon Apr 19 22:08:44 2004
 * \brief  Release function for the RTT_Format_Reader library
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef RTT_Format_Reader_Release_hh
#define RTT_Format_Reader_Release_hh

#include <string>

//===========================================================================//
/*!
 * \namespace rtt_RTT_Format_Reader
 * 
 * \brief Namespace to contain the RTT_Format_Reader utilities.
 *
 * Provides namespace protection for the Draco RTT_Format_Reader
 * utilities. The RTT_Format_Reader class constructor automatically
 * instantiates and executes the readMesh member function used to parse the
 * mesh data.  Accessor functions are provided for all of the remaining
 * member classes to allow data retrieval.
 *
 * The rtt_RTT_Format_Reader::RTT_Format_Reader class is also wrapped into
 * the rtt_meshReader namespace.  Inside of the meshReaders package there is
 * a file, RTT_Mesh_Reader.hh, that one can include.  This file puts the
 * RTT_Mesh_Reader class inside of the rtt_meshReaders namespace for
 * convenience.  It does not violate levelization because this file is only a
 * header.  It is only included in an upper level component when both the
 * meshReaders and RTT_Format_Reader packages are already included.
 *
 * \sa The <a href="./index.html">Main Page</a> presents a summary of the
 * capabilities provided by the namespace.
 */
//===========================================================================//

namespace rtt_RTT_Format_Reader
{
    const std::string release();
}

#endif // RTT_Format_Reader_Release_hh

//---------------------------------------------------------------------------//
//                        end of RTT_Format_Reader/Release.hh
//---------------------------------------------------------------------------//
