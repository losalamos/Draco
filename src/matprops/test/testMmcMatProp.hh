//----------------------------------*-C++-*----------------------------------//
// testMmcMatProp.hh
// John McGhee
// Mon Sep 21 09:53:46 1998
//---------------------------------------------------------------------------//
// @> Test facility for the Multi-Material Cell Mat Props class.
//---------------------------------------------------------------------------//

#ifndef __matprops_test_testMmcMatProp_hh__
#define __matprops_test_testMmcMatProp_hh__

#include <vector>

//===========================================================================//
// class testMmcMatProp - Tests the Multi-Material Cell Material Properties.

// Note that this test code only verifies the proper operation of 
// the Multi-Material Cell class. It does no checking of the imbeded
// Uni-Material Cell mat props class.

//===========================================================================//

class testMmcMatProp {

  public:
    testMmcMatProp();
    ~testMmcMatProp();

    void execute_test();

    template<class T>
    void CellAvgResultsOK(const T &results, 
			  const std::vector<double> &answer, 
			  const double eps, bool &pass) const;

    template <class T2, class T3>
    void ByMatResultsOK(const T2 &results, 
			const std::vector<std::vector<T3> > &answer, 
			const double eps, bool &pass) const;

    template <class T2>
    void ByMatResultsOK(const T2 &results, 
			const std::vector<std::vector<int> > &answer, 
			const double eps, bool &pass) const;

};

#endif                          // __matprops_test_testMmcMatProp_hh__

//---------------------------------------------------------------------------//
//                              end of matprops/test/testMmcMatProp.hh
//---------------------------------------------------------------------------//
