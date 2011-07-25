//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   device/test/dacs_wait_for_cmd_ppe.cc
 * \author Gabriel M. Rockefeller
 * \date   Wed Jul 20 09:16:12 2011
 * \brief  A tiny DACS accel-side binary, for testing DACS_Process.
 * \note   Copyright (C) 2011 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "ds++/Assert.hh"

#include <dacs.h>
#include <unistd.h>

int main()
{
    DACS_ERR_T err = dacs_init(DACS_INIT_FLAGS_NONE);
    Insist(err == DACS_SUCCESS, dacs_strerror(err));

    // wait for command
    dacs_wid_t wid;
    err = dacs_wid_reserve(&wid);
    Insist(err == DACS_SUCCESS, dacs_strerror(err));

    uint32_t command;
    err = dacs_recv(&command,
                    sizeof(uint32_t),
                    DACS_DE_PARENT,
                    DACS_PID_PARENT,
                    1,
                    wid,
                    DACS_BYTE_SWAP_DISABLE);
    Insist(err == DACS_SUCCESS, dacs_strerror(err));

    err = dacs_wait(wid);
    Insist(err == DACS_WID_READY, dacs_strerror(err));

    Insist(command == 2121256449, "Failed to receive expected command");

    err = dacs_wid_release(&wid);
    Insist(err == DACS_SUCCESS, dacs_strerror(err));

    err = dacs_exit();
    Insist(err == DACS_SUCCESS, dacs_strerror(err));

    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of dacs_wait_for_cmd_ppe.cc
//---------------------------------------------------------------------------//
