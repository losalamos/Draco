//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/test/tGandolfFile.hh
 * \author Kelly Thompson
 * \date   Tue Aug 22 15:48:55 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_gandolf_tGandolfFile_hh__
#define __cdi_gandolf_tGandolfFile_hh__

#include "UnitTestFrame/TestApp.hh"

namespace rtt_cdi_gandolf_test
{
 
//===========================================================================//
/*!
 * \class tGandolfFile
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class tGandolfFile : public rtt_UnitTestFrame::TestApp
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA
    
  public:

    // CREATORS
    
    tGandolfFile( int argc, char *argv[], std::ostream &os_in );
    // (defaulted) tGandolfFile(const tGandolfFile &rhs);
    // (defaulted) ~tGandolfFile();

    // MANIPULATORS
    
    // (defaulted) tGandolfFile& operator=(const tGandolfFile &rhs);

    // ACCESSORS

    std::string name() const { return "tGandolfFile"; }
    std::string version() const;

  protected:

    std::string runTest();
    
  private:

    // IMPLEMENTATION
};

} // end namespace rtt_cdi_gandolf_test

#endif // __cdi_gandolf_tGandolfFile_hh__

//---------------------------------------------------------------------------//
//                 end of cdi_gandolf/tGandolfFile.hh
//---------------------------------------------------------------------------//
