//----------------------------------*-C++-*----------------------------------//
// TestMatPropReader.hh
// Randy M. Roberts
// Mon Apr 20 15:55:22 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __matprops_test_TestMatPropReader_hh__
#define __matprops_test_TestMatPropReader_hh__

#include <string>
#include <vector>

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

BEGIN_NS_XTM
    
// Forward Declaration

class Units;

//===========================================================================//
// class TestMatPropReader - 
//
// Date created :
// Purpose      :
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class TestMatPropReader
{
    // NESTED CLASSES AND TYPEDEFS

    // DATA

    // DISALLOWED DEFAULT METHODS
    
  private:
    
    TestMatPropReader(const TestMatPropReader &rhs);
    TestMatPropReader& operator=(const TestMatPropReader &rhs);

  public:

    // CREATORS
    
    TestMatPropReader(const Units &units, const std::string &filename,
		      const std::vector<int> &materialIds);
    
    // MANIPULATORS
    
    // ACCESSORS

  private:
    
    // IMPLEMENTATION
};


END_NS_XTM  // namespace XTM

#endif                          // __matprops_test_TestMatPropReader_hh__

//---------------------------------------------------------------------------//
//                              end of matprops/test/TestMatPropReader.hh
//---------------------------------------------------------------------------//
