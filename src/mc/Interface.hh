//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Interface.hh
 * \author Thomas M. Evans
 * \date   Mon Jul 10 16:38:49 2000
 * \brief  Interface abstract base class definition for the mc package.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __mc_Interface_hh__
#define __mc_Interface_hh__

#include <string>

namespace rtt_mc
{
 
//===========================================================================//
/*!
 * \class Interface

 * \brief Interface base class to rtt_mc package components.

 * The interface base class defines the interface functionality required by
 * the mesh builders in the rtt_mc package.  At this point the interface has
 * not been formally defined; however, this is a goal in the immediate
 * future.

 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class Interface 
{
  public:
    // Useful typedefs.
    typedef std::string std_string;

  public:
    //! Constructor.
    Interface() { /* no data to construct */ }
   
    //! Virtual constructor to make life happy down the inheritance chain.
    virtual ~Interface() { /* need a destructor for inhertance chain */ }

    // >>> FUNCTIONS REQUIRED BY ALL MESH BUILDERS

    //! Get the name of the mesh description file.
    std_string get_mesh_file() const = 0;
};

} // end namespace rtt_mc

#endif                          // __mc_Interface_hh__

//---------------------------------------------------------------------------//
//                              end of mc/Interface.hh
//---------------------------------------------------------------------------//
