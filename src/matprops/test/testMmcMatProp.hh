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

// 
//===========================================================================//

class testMmcMatProp {

  public:
    testMmcMatProp();
    ~testMmcMatProp();

    void execute_test();

    template<class T>
    void CellAvgResultsOK(T &results, 
			  std::vector<double> &answer, 
			  double eps, bool &pass);

    template <class T2>
    void ByMatResultsOK(T2 &results, 
			std::vector<std::vector<double> > &answer, 
			double eps, bool &pass);

};

#endif                          // __matprops_test_testMmcMatProp_hh__

//---------------------------------------------------------------------------//
//                              end of matprops/test/testMmcMatProp.hh
//---------------------------------------------------------------------------//
