//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   UnitTestFrame/C4ProcessPolicy.hh
 * \author Randy M. Roberts
 * \date   Wed Nov 29 10:16:22 2000
 * \brief  Process Policy class for initing and finalizing processes w/ C4.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __UnitTestFrame_C4ProcessPolicy_hh__
#define __UnitTestFrame_C4ProcessPolicy_hh__

#include "ProcessPolicy.hh"
#include "c4/global.hh"

namespace rtt_UnitTestFrame
{

/*!
 * Define ProcessInit in the namespace.
 * This function will call C4::Init
 * to initialize the process.
 */

void ProcessInit(int &argc, char *argv[])
{
    C4::Init(argc, argv);
}

/*!
 * Define ProcessFinalize in the namespace.
 * This function will call C4::Finalize
 * to finalize the process.
 */

void ProcessFinalize()
{
    C4::Finalize();
}

} // end namespace rtt_UnitTestFrame

#endif                          // __UnitTestFrame_C4ProcessPolicy_hh__

//---------------------------------------------------------------------------//
//                              end of UnitTestFrame/C4ProcessPolicy.hh
//---------------------------------------------------------------------------//
