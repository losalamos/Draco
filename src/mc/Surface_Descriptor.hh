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
 * \example mc/test/tstSurface_Descriptor.cc
 * 
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

    // DATA

    Surface_Type type;
    std::vector<double> data;

    //! Default constructor
    Surface_Descriptor() { /* ... */ }

    //! Destructor
    ~Surface_Descriptor() { /* ... */ }

    //! Constructor
    Surface_Descriptor(int type, const std::vector<double>& data);

    //! Unpacking constructor
    Surface_Descriptor(const std::vector<char> &data);

    //! Packing operator
    std::vector<char> pack() const;

    bool operator==(const Surface_Descriptor& rhs) const;
    bool operator!=(const Surface_Descriptor& rhs) const;

};

} // end namespace rtt_mc

#endif // rtt_mc_Surface_Descriptor_hh

//---------------------------------------------------------------------------//
//              end of mc/Surface_Descriptor.hh
//---------------------------------------------------------------------------//
