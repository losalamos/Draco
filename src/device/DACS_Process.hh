//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   device/DACS_Process.hh
 * \author Gabriel M. Rockefeller
 * \date   Thu Jul 14 16:15:29 2011
 * \brief  Define class DACS_Process
 * \note   Copyright (C) 2011-2014 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef device_DACS_Process_hh
#define device_DACS_Process_hh

#include <dacs.h>
#include <string>

namespace rtt_device
{

//===========================================================================//
/*!
 * \class DACS_Process
 * \brief A abstract host-side representation of an accel-side process.
 *
 * DACS_Process represents an accel-side process.  It provides methods to
 * start and stop the accel-side process, and accessors that return the
 * filename of the accel-side binary and the pid of the running process.
 *
 * Clients should inherit from DACS_Process and implement the stop method.
 *
 * \sa DACS_Process.cc for detailed descriptions.
 * \sa DACS_External_Process.hh for a concrete implementation.
 */
/*!
 * \example device/test/tstDACS_Device_Process.cc
 *
 * Test of DACS_Process.
 */
//===========================================================================//

class DACS_Process 
{
  public:

    //! Constructor
    DACS_Process(const std::string & filename);

    //! Destructor.
    virtual ~DACS_Process() { }

    //! Comparison operator.
    bool operator==(const DACS_Process &other) const
    {
        return (accel_filename == other.accel_filename);
    }

    //! Comparison operator.
    bool operator!=(const DACS_Process &other) const
    {
        return (accel_filename != other.accel_filename);
    }
    
    // SERVICES

    //! Launch the accel-side process represented by DACS_Process.
    virtual void start(const de_id_t & de_id);

    //! Terminate the accel-side process.
    virtual void stop() = 0;

    // ACCESSORS

    //! Return the filename of the binary associated with this process.
    std::string const & get_filename() const { return accel_filename; }

    //! Return the pid of the running accel-side process.
    dacs_process_id_t const & get_pid() const { return pid; }

  private:

    // DATA

    std::string accel_filename;
    dacs_process_id_t pid;
};

} // end namespace rtt_device

#endif // device_DACS_Process_hh

//---------------------------------------------------------------------------//
// end of device/DACS_Process.hh
//---------------------------------------------------------------------------//
