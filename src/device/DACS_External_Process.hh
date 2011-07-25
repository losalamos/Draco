//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   device/DACS_External_Process.hh
 * \author Gabriel M. Rockefeller
 * \date   Thu Jul 14 16:03:44 2011
 * \brief  Define class DACS_External_Process
 * \note   Copyright (C) 2007 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef device_DACS_External_Process_hh
#define device_DACS_External_Process_hh

#include "DACS_Process.hh"

namespace rtt_device
{

//===========================================================================//
/*!
 * \class DACS_External_Process
 * \brief A DACS_Process that accepts an external "stop" function pointer.
 *
 * DACS_External_Process derives from DACS_Process and provides an
 * implementation of the stop method that invokes a function passed to the
 * constructor to stop the accel-side process.
 */
/*! 
 * \example device/test/tstDACS_External_Process.cc
 *
 * Test of DACS_External_Process.
 *
 * \example device/test/tstDACS_Device_Process.cc
 *
 * A combined test of DACS_External_Process and DACS_Device, through the flat
 * interface.
 */
//===========================================================================//

class DACS_External_Process : public DACS_Process
{
  public:

    // TYPEDEFS

    typedef void (*voidfunc)();
    
    //! Default constructor.
    DACS_External_Process(const std::string & filename,
                          voidfunc final_func) : DACS_Process(filename),
                                                 f(final_func) { }

    // SERVICES

    //! Terminate the accel-side process.
    virtual void stop()
    {
        if (f)
            f();
    }

  private:

    // DATA

    voidfunc f;
};

} // end namespace rtt_device

#endif // device_DACS_External_Process_hh

//---------------------------------------------------------------------------//
//              end of device/DACS_External_Process.hh
//---------------------------------------------------------------------------//
