//----------------------------------*-C++-*----------------------------------//
// testMmcMatProp.t.cc
// John McGhee
// Fri Sep 25 07:58:49 1998
//---------------------------------------------------------------------------//
// @> Template code for the Multi-Material Cell mat props tester.
//---------------------------------------------------------------------------//

#include "matprops/test/testMmcMatProp.hh"
#include <string>

template <class T>
void  testMmcMatProp::CellAvgResultsOK(T &results, 
				       std::vector<double> &answer, 
				       double eps, bool &pass)
{
    std::string error_flag;
    int ncell = results.size();
    Require(ncell == answer.size());
    std::vector<double> relErr(ncell);
    double relErrMax = 0.;
    T::iterator resit = results.begin();
    for( int icell = 0; icell<ncell; icell++)
    {
	relErr[icell] = std::abs(*resit-answer[icell])/answer[icell];
	if (relErr[icell] > relErrMax) relErrMax = relErr[icell];
	if (relErr[icell] > eps)
	{
	    error_flag = " <=== FAIL";
	}
	else
	{
	    error_flag = " ";
	}
	std::cout << " [" << icell << "]   " << 
	    *resit++ << " " << relErr[icell] << 
	    error_flag << std::endl;
    }
    std::cout << std::endl;
    pass =  relErrMax <= eps;
}



template <class T2>
void testMmcMatProp::ByMatResultsOK(T2 &results, 
				    std::vector<std::vector<double> > &answer, 
				    double eps, bool &pass)
{
    std::string error_flag;
    int ncell = results.size();
    Require(ncell == answer.size());
    std::vector<std::vector<double> > relErr(ncell);
    double relErrMax = 0.; 
    T2::iterator resit = results.begin();
    for( int icell = 0; icell<ncell; icell++)
    {
	int nmat = (*resit).size();
	relErr[icell].resize(nmat);
	Require(nmat == answer[icell].size());
	T2::value_type::iterator resitit= (*resit++).begin();
	for (int imat = 0; imat<nmat; imat++)
	{
	    relErr[icell][imat] = 
		std::abs(*resitit-answer[icell][imat])/
		answer[icell][imat];
	    if (relErr[icell][imat] > relErrMax)
		relErrMax = relErr[icell][imat];
	    if (relErr[icell][imat] > eps)
	    {
		error_flag = " <=== FAIL";
	    }
	    else
	    {
		error_flag = " ";
	    }
	    std::cout << " [" << icell << 
		"]  [" << imat << "]   " << 
		*resitit++ << " " << relErr[icell][imat] << 
		error_flag << std::endl;
	}
	std::cout << std::endl;
    }
    std::cout << std::endl;
    pass =  relErrMax <= eps;
}


//---------------------------------------------------------------------------//
//                              end of testMmcMatProp.t.cc
//---------------------------------------------------------------------------//


