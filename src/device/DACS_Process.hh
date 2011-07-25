//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   device/DACS_Process.hh
 * \author Gabriel M. Rockefeller
 * \date   Thu Jul 14 16:15:29 2011
 * \brief  Define class DACS_Process
 * \note   Copyright (C) 2007 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef device_DACS_Process_hh
#define device_DACS_Process_hh

#include <string>

#include <dacs.h>

namespace rtt_device
{

//===========================================================================//
/*!
 * \class DACS_Process
 * \brief
 *
 * Long description or discussion goes here.  Information about Doxygen
 * commands can be found at http://www.doxygen.org.
 *
 * \sa DACS_Process.cc for detailed descriptions.
 *
 * \par Code Sample:
 * \code
 *     cout << "Hello, world." << endl;
 * \endcode
 */
/*! 
 * \example device/test/tstDACS_Process.cc
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

    std::string const & get_filename() const { return accel_filename; }
    dacs_process_id_t const & get_pid() const { return pid; }

  private:

    // DATA

    std::string accel_filename;
    dacs_process_id_t pid;
};

} // end namespace rtt_device

#endif // device_DACS_Process_hh

//---------------------------------------------------------------------------//
//              end of device/DACS_Process.hh
//---------------------------------------------------------------------------//
