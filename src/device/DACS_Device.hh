//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   device/DACS_Device.hh
 * \author Gabriel M. Rockefeller
 * \date   Thu Jun 16 14:55:48 2011
 * \brief  DACS_Device header file.
 * \note   Copyright (C) 2011 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef device_DACS_Device_hh
#define device_DACS_Device_hh

#include "ds++/Assert.hh"

#include <dacs.h>
#include <string>

namespace rtt_device
{

//===========================================================================//
/*!
 * \class DACS_Device
 * \brief A singleton for DACS host-side resource and process management.
 *
 * DACS_Device provides a minimal representation of the host side of a hybrid
 * DACS process.  It initializes DACS, reserves a Cell, launches the
 * accel-side process, and provides accessors for the DACS de_id and child
 * pid.
 *
 * Because it's implemented as a (Meyers) singleton, DACS_Device ensures that
 * DACS initialization, reservation, and process-launch only happen once.  It
 * also waits for the accel-side process to complete and calls dacs_exit at
 * program termination.
 *
 * \sa DACS_Device.cc for detailed descriptions.
 */
/*! 
 * \example device/test/tstDACS_Device.cc
 *
 * Test of DACS_Device.
 */
//===========================================================================//

class DACS_Device
{
  public:

    //! Return a reference to the DACS_Device.
    static DACS_Device& instance()
    {
        static DACS_Device the_instance;
        return the_instance;
    }

    //! Set the filename of the accel-side binary.
    static void init(const std::string & fname);

    //! Return the element id of the reserved DACS child.
    de_id_t get_de_id() { return de_id; }

    //! Return the process id of the DE process.
    dacs_process_id_t get_pid() { return pid; }
    
  private:

    //! Private default constructor; DACS_Device constructs itself.
    DACS_Device();

    //! Private destructor; the DACS_Device instance is destroyed by atexit.
    ~DACS_Device();

    //! Private copy constructor; DACS_Device is unique.
    DACS_Device(const DACS_Device&);

    //! Private assignment operator; any assignment would be self-assignment.
    DACS_Device& operator=(const DACS_Device&);

    // DATA

    static std::string accel_filename;
    de_id_t de_id;
    dacs_process_id_t pid;
};

} // end namespace rtt_device

#endif // device_DACS_Device_hh

//---------------------------------------------------------------------------//
//              end of device/DACS_Device.hh
//---------------------------------------------------------------------------//
