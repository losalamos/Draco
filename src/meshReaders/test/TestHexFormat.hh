//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshReaders/test/TestHexFormat.hh
 * \author John McGhee
 * \date   Thu Mar  9 08:54:59 2000
 * \brief  Header file for the Hex_Format class unit test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __meshReaders_test_TestHexFormat_hh__
#define __meshReaders_test_TestHexFormat_hh__

#include "UnitTestFrame/TestApp.hh"

namespace rtt_meshReaders_test
{
 
//===========================================================================//
/*!
 * \class TestHexFormat
 *
 * \brief Tests the Hex_Format mesh reader class.
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class TestHexFormat : public rtt_UnitTestFrame::TestApp  
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA
    
  public:

    // CREATORS
    
    TestHexFormat(int argc, char *argv[], std::ostream &os_in);
    // Defaulted TestHexFormat(const TestHexFormat &rhs);
    // Defaulted ~TestHexFormat();

    // MANIPULATORS
    
    //Defaulted TestHexFormat& operator=(const TestHexFormat &rhs);

    // ACCESSORS

    std::string name() const { return "TestHexFormat"; }

    std:: string version() const;

  protected:

    std::string runTest();
    
  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_meshReaders_test

#endif                // __meshReaders_test_TestHexFormat_hh__

//---------------------------------------------------------------------------//
//                    end of meshReaders/test/TestHexFormat.hh
//---------------------------------------------------------------------------//
