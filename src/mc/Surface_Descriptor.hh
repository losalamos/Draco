//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Surface_Descriptor.hh
 * \author Mike Buksas
 * \date   Mon Aug 11 13:17:49 2003
 * \brief  Header file for Surface_Descriptor
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_mc_Surface_Descriptor_hh
#define rtt_mc_Surface_Descriptor_hh

#include <vector> 

namespace rtt_mc
{

//===========================================================================//
/*!
 * \class Surface_Descriptor
 * \brief A Struct for storing a surface description.
 *
 * This class enumerates the available surface types (exactly one, so far)
 * and their data storage requirements.
 *
 */
/*! 
 * \example mc/test/mc_test.cc 
 * 
 * description of example
 */
// revision history:
// -----------------
// 0) (Mon Aug 11 13:17:49 2003) Mike Buksas: original
// 
//===========================================================================//

struct Surface_Descriptor
{

    enum Surface_Type { SPHERE = 0 };

    static const int kinds = 1;
    static const int sizes[1];

    Surface_Type type;

    std::vector<double> data;

};

} // end namespace rtt_mc

#endif // rtt_mc_Surface_Descriptor_hh

//---------------------------------------------------------------------------//
//              end of mc/Surface_Descriptor.hh
//---------------------------------------------------------------------------//
