//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Surface_Tracking_Interface.hh
 * \author Mike Buksas
 * \date   Wed Jul 30 12:26:37 2003
 * \brief  Header file for Surface_Tracking_Interface
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_imc_Surface_Tracking_Interface_hh
#define rtt_imc_Surface_Tracking_Interface_hh

#include "mc/Surface_Descriptor.hh"

namespace rtt_imc
{

//===========================================================================//
/*!
 * \class Surface_Tracking_Interface
 * \brief An abstract interface representing the operations necessary to
 * build surface tracking applications.
 *
 *
 */
// revision history:
// -----------------
// 0) (Wed Jul 30 12:26:37 2003) Mike Buksas: original
// 
//===========================================================================//

class Surface_Tracking_Interface 
{
  protected:

    //! protected constructor
    Surface_Tracking_Interface() { /* No data */ }

  public:

    //! destructor
    virtual ~Surface_Tracking_Interface() { /* ... */ }

    // ACCESSORS

    //! Get the number of surfaces
    virtual int number_of_surfaces() const = 0;
    
    //! Get the entire vector of descriptors:
    virtual const std::vector<rtt_mc::Surface_Descriptor>& get_surface_data() 
	const = 0;

    //! Get the vector of angular bin data
    virtual const std::vector<double>& get_bin_cosines() const = 0;
    
};

} // end namespace rtt_imc

#endif // rtt_imc_Surface_Tracking_Interface_hh

//---------------------------------------------------------------------------//
//              end of imc/Surface_Tracking_Interface.hh
//---------------------------------------------------------------------------//
