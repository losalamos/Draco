//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshReaders/test/TestElementDefinition.hh
 * \author John McGhee
 * \date   Fri Mar  3 08:41:46 2000
 * \brief  Header file for the Element_Definition class unit test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __meshReaders_test_TestElementDefinition_hh__
#define __meshReaders_test_TestElementDefinition_hh__

#include "UnitTestFrame/TestApp.hh"

namespace rtt_meshReaders_test
{
 
//===========================================================================//
/*!
 * \class TestElementDefinition
 *
 * \brief Tests the Element_Definition class.
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class TestElementDefinition : public rtt_UnitTestFrame::TestApp 
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA
    
  public:

    // CREATORS
    
    TestElementDefinition(int argc, char *argv[], std::ostream &os_in);

    // Defaulted TestElementDefinition(const TestElementDefinition &rhs);
    // Defaulted ~TestElementDefinition();

    // MANIPULATORS
    
    // Defaulted TestElementDefinition& operator=(
    //            const TestElementDefinition &rhs);

    // ACCESSORS
 
    std::string name() const { return "TestElementDefinition"; }

    std:: string version() const;

  protected:

    std::string runTest();

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_meshReaders_test

#endif           // __meshReaders_test_TestElementDefinition_hh__

//---------------------------------------------------------------------------//
//               end of meshReaders/test/TestElementDefinition.hh
//---------------------------------------------------------------------------//
