//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   UnitTestFrame/ScalarProcessPolicy.hh
 * \author Randy M. Roberts
 * \date   Wed Nov 29 10:16:22 2000
 * \brief  Process Policy class for initing and finalizing processes w/out C4.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __UnitTestFrame_ScalarProcessPolicy_hh__
#define __UnitTestFrame_ScalarProcessPolicy_hh__

#include "ProcessPolicy.hh"

namespace rtt_UnitTestFrame
{

/*!
 * Define ProcessInit in the namespace.
 * This function will call C4::Init
 * to initialize the process.
 */

void ProcessInit(int &argc, char *argv[])
{
    /* empty */
}

/*!
 * Define ProcessFinalize in the namespace.
 * This function will call C4::Finalize
 * to finalize the process.
 */

void ProcessFinalize()
{
    /* empty */
}

} // end namespace rtt_UnitTestFrame

#endif                          // __UnitTestFrame_ScalarProcessPolicy_hh__

//---------------------------------------------------------------------------//
//                              end of UnitTestFrame/ScalarProcessPolicy.hh
//---------------------------------------------------------------------------//
