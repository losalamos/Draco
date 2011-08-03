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

#include "DACS_Process.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"

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
 * DACS process.  It initializes DACS, reserves a Cell, launches an accel-side
 * process, and provides accessors for the DACS de_id and child pid.
 *
 * Because it's implemented as a (Meyers) singleton, DACS_Device ensures that
 * DACS initialization, reservation, and process-launch only happen once.  It
 * also tells the accel-side process when to stop, waits for the process to
 * complete, and calls dacs_exit at program termination.
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

    // TYPEDEFS

    typedef rtt_dsxx::SP<DACS_Process> SP_Proc;

    //! Return a reference to the DACS_Device.
    static DACS_Device& instance()
    {
        static DACS_Device the_instance;
        return the_instance;
    }

    //! Specify a DACS_Process to be launched on the device.
    static void init(SP_Proc proc);

    //! Return the element id of the reserved DACS child.
    de_id_t get_de_id() const { return de_id; }

    //! Return the process id of the DE process.
    dacs_process_id_t get_pid() const { return process->get_pid(); }

    //! Return the handle to the combined host-accel DACS group.
    dacs_group_t get_group() const { return group; }
    
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

    static SP_Proc process;
    de_id_t de_id;
    dacs_group_t group;
};

} // end namespace rtt_device

#endif // device_DACS_Device_hh

//---------------------------------------------------------------------------//
//              end of device/DACS_Device.hh
//---------------------------------------------------------------------------//
