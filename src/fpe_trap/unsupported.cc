//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   fpe_trap/unsupported.cc
 * \author Rob Lowrie
 * \date   Thu Oct 13 16:52:05 2005
 * \brief  Default implementation of fpe_trap functions.
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <fpe_trap/config.h>

#ifdef FPETRAP_UNSUPPORTED

namespace rtt_fpe_trap
{

bool enable_fpe()
{
    return false;
}

} // end namespace rtt_shared_lib

#endif // FPETRAP_UNSUPPORTED

//---------------------------------------------------------------------------//
//                 end of unsupported.cc
//---------------------------------------------------------------------------//
