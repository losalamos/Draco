//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   UnitTestFrame/ProcessPolicy.hh
 * \author Randy M. Roberts
 * \date   Wed Nov 29 11:13:12 2000
 * \brief  Declarations for ProcessPolicies.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __UnitTestFrame_ProcessPolicy_hh__
#define __UnitTestFrame_ProcessPolicy_hh__

namespace rtt_UnitTestFrame
{
 
/*!
 * Declare ProcessInit in the namespace.
 * This function will call C4::Init
 * to initialize the process.
 */

void ProcessInit(int &argc, char *argv[]);

/*!
 * Declare ProcessFinalize in the namespace.
 * This function will call C4::Finalize
 * to finalize the process.
 */

void ProcessFinalize();

} // end namespace rtt_UnitTestFrame

#endif                          // __UnitTestFrame_ProcessPolicy_hh__

//---------------------------------------------------------------------------//
//                              end of UnitTestFrame/ProcessPolicy.hh
//---------------------------------------------------------------------------//
